###################################################################
# electrostat_analysis.R
#
# Analysis relating to pairwise electrostatic interactions
###################################################################

###################################################################
# Common path utilities
###################################################################

# Cached directory to script
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
SCRIPT_DIR = getSrcDirectory(function(x) {
  x
})

# @return full path to directory containing this script
get_script_dir <- function() {
  return(SCRIPT_DIR)
}

###################################################################
# Globals
###################################################################

# Load globals
source(file.path(get_script_dir(), "electrostat_globals.R"))

###################################################################
# Random forest feature analysis
###################################################################

library(randomForest)

# Number of trees in the forest
NUM_RF_TREES = 10000

# If true, observed values are scaled by a factor to 'help' distinguish
# from NA values.
RF_SHOULD_SCALE_OBSERVED = TRUE

# Version number for random forest data
# 1 - no in-bag matrix
RF_VERSION_NO_IN_BAG = 1
# 2 - includes in-bag matrix
RF_VERSION_IN_BAG = 2
# 3 - includes in-bag matrix and does not filter intra-loop interactions
# only relevant for interaction random forests
RF_VERSION_INTRA_LOOP = 3
# Default interaction random forest version
# DOES NOT APPLY TO SRS (single residue scores)
RF_INT_VERSION = RF_VERSION_IN_BAG

# Default interaction keep setting for random forest
# DOES NOT APPLY TO SRS (single residue scores)
# Possible values are:
#   "ALL.ION" - keep all loop ion to any other ion interactions
#   "LOOP.ION" - keep only loop-to-loop ion interactions
#   "XORLOOP.ION" - keep only loop-to-non-loop ion interaction
#   "ALL.ION.POLAR" - keep all ion-ion, ion-polar, polar-polar interactions
# Can be vector like so: c("ALL.ION", "LOOP.ION", "XORLOOP.ION") and methods
# will correctly either use all args or only take the first arg.
RF_INT_KEEP = "ALL.ION"

# @return Path to random forest rdata
get_rf_path <- function(base_output_dir = get_relax_c_beta_dir(),
                        pH = 5,
                        sim_id = "wt",
                        keep = RF_INT_KEEP,
                        num_trees = NUM_RF_TREES,
                        should_scale_observed = RF_SHOULD_SCALE_OBSERVED,
                        rf_version = RF_INT_VERSION,
                        max_energy_rank = MAX_ENERGY_RANK) {
  rdata_path = file.path(
    base_output_dir,
    "capt",
    paste0("pH", pH),
    sim_id,
    paste0(
      "rf.",
      keep[1],
      ".ntree",
      num_trees,
      ".r",
      MAX_ENERGY_RANK,
      ".rdata"
    )
  )
  # @HACK - patch rdata_path with scale status
  if (!should_scale_observed) {
    rdata_path = insert_attrib_in_rdata(rdata_path, "no_scale")
  }
  
  # Check version
  if (rf_version > 1) {
    rdata_path = insert_attrib_in_rdata(rdata_path, paste0("v", rf_version))
  }
  return(rdata_path)
}

# WT RCB
WT_PH5_RF_PATH_ALL_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 5, "wt", "ALL.ION")
WT_PH7_RF_PATH_ALL_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 7, "wt", "ALL.ION")
WT_PH5_RF_PATH_LOOP_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 5, "wt", "LOOP.ION")
WT_PH7_RF_PATH_LOOP_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 7, "wt", "LOOP.ION")
WT_PH5_RF_PATH_XORLOOP_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 5, "wt", "XORLOOP.ION")
WT_PH7_RF_PATH_XORLOOP_ION_RCB = get_rf_path(get_relax_c_beta_dir(), 7, "wt", "XORLOOP.ION")

# Default random forest interaction filter
# @param elec - List with members $design and $meta, see load_data_utils.R:
#   elec_csv_2_rdata() for details
# @param keep - Specifies the interaction filter settings (only first element is used)
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
#     ALL.ION.POLAR - all ion-ion, ion-polar, polar-polar interactions are kept
# @param rf_version - Version of forest data
# @return Logical with elements to keep
rf_int_filt <- function(elec,
                        keep = RF_INT_KEEP,
                        rf_version = RF_INT_VERSION) {
  keep = keep[1]
  # Make sure keep option is recognized
  stopifnot(keep %in% c("ALL.ION", "LOOP.ION", "XORLOOP.ION", "ALL.ION.POLAR"))
  
  # Filter covariates so only ion-ion & loop-loop interactions are kept
  # Default to "ALL.ION" - all ion-ion interactions are kept and at least one residue is part of a loop
  filt = elec$meta$ION.ION &
    ((elec$meta$RESA.LOOP != "NONE") |
       (elec$meta$RESB.LOOP != "NONE"))
  
  keep = keep[1]
  if (keep == "LOOP.ION") {
    # Only ion-ion interactions between loop regions are kept
    print("...keeping ion-ion interactions between loop regions")
    filt = elec$meta$LOOP.ION
  } else if (keep == "XORLOOP.ION") {
    print("...keeping ion-ion interactions from loop to non-loop regions")
    # Only ion-ion interactions from loop to non-loop regions are kept
    filt = elec$meta$ION.ION &
      xor((elec$meta$RESA.LOOP != "NONE"),
          (elec$meta$RESB.LOOP != "NONE"))
  } else if (keep == "ALL.ION.POLAR") {
    # All ion-ion, ion-polar, polar-polar interactions kept
    print("...keeping ion-ion, ion-polar, polar-polar interactions")
    filt.resa.ion.polar = elec$meta$RESA.ION | elec$meta$RESA.POL
    filt.resb.ion.polar = elec$meta$RESB.ION | elec$meta$RESB.POL
    filt = filt.resa.ion.polar | filt.resb.ion.polar
  }
  
  # Also, only keep inter-loop interactions
  if (rf_version < 3) {
    print("...removing intra-loop interactions")
    filt = filt & (elec$meta$RESA.LOOP != elec$meta$RESB.LOOP)
  }
  
  return(filt)
}

# Default random forest SRS filter
# @param elec - List with members $design and $meta, see load_data_utils.R:
#   get_srs() for details (members are similar to those in elec_csv_2_rdata())
# @param keep - Specifies the interaction filter settings (only first element is used)
#     SRS.ALL - all residues kept
#     SRS.ION - only ionic residues kept
#     SRS.POLAR - only polar residues kept
#     SRS.ION.POLAR - only ionic or polar residues kept
#     SRS.LOOP.ION - only ionic residues within loops are kept
#     SRS.LOOP.ION.POLAR - only ionic or polar residues within loops are kept
# @param rf_version - Version of forest data
# @return Logical with elements to keep
rf_srs_filt <- function(elec,
                        keep,
                        rf_version) {
  keep = keep[1]
  # Make sure keep option is recognized
  stopifnot(
    keep %in% c(
      "SRS.ALL",
      "SRS.ION",
      "SRS.POLAR",
      "SRS.ION.POLAR",
      "SRS.LOOP.ION",
      "SRS.LOOP.ION.POLAR"
    )
  )
  # Default all residues considered
  # should be called "SRS.ALL"
  filt = rep(TRUE, nrow(elec$meta))
  
  # Only ionic residues allowed
  if (keep == "SRS.ION") {
    print("... keeping ionic residues")
    filt = elec$meta$RES.ION
  } else if (keep == "SRS.POLAR") {
    print("... keeping polar residues")
    filt = elec$meta$RES.POLAR
  } else if (keep == "SRS.ION.POLAR") {
    print("... keeping ionic and polar residues")
    # Only ionic or polar residues allowed
    filt = elec$meta$RES.ION | elec$meta$RES.POL
  } else if (keep == "SRS.LOOP.ION") {
    print("... keeping ionic loop residues")
    filt = elec$meta$RES.ION & (elec$meta$RES.LOOP != "NONE")
  } else if (keep == "SRS.LOOP.ION.POLAR") {
    print(".. keeping ionic and polar loop residues")
    filt = elec$meta$RES.ION | elec$meta$RES.POL
    filt = filt & (elec$meta$RES.LOOP != "NONE")
  }
  
  return(filt)
}

# @param elec - List with members $design and $meta, see load_data_utils.R:
#   elec_csv_2_rdata() for details
# @param keep - Specifies feature filter settings (only first element is used)
# - for interaction data:
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
#     ALL.ION.POLAR - all ion-ion, ion-polar, polar-polar interactions are kept
# - for single residue data:
#     SRS.ALL - all residues kept
#     SRS.ION - only ionic residues kept
#     SRS.POLAR - only polar residues kept
#     SRS.ION.POLAR - only ionic or polar residues kept
# @param rf_version - Version of forest data
# @param should_scale_observed - If true, observed values will be scaled by a
#   'large' factor to possibly help differentiate from non-observed (NA) values
# @return data.frame that can serve to train a random forest model
get_rf_design <- function(elec,
                          keep = RF_INT_KEEP,
                          rf_version = RF_INT_VERSION,
                          should_scale_observed = RF_SHOULD_SCALE_OBSERVED) {
  # Determine which covariates to keep
  keep = keep[1]
  filt = c()
  if (grepl(pattern = "SRS",
            x = keep,
            fixed = TRUE)) {
    # single residue scoring filter
    filt = rf_srs_filt(elec, keep, rf_version)
  } else {
    # interaction filter
    filt = rf_int_filt(elec, keep, rf_version)
  }
  
  keep_names = rownames(elec$meta)[filt]
  x = elec$design[, keep_names]
  
  # Random forest does not work well with missing feature values. There is
  # rfImpute but we know that the feature value is missing because the
  # interaction distance is too far to be recorded, in which case it's
  # treated as having 0 energy contribution. So, let's just set all
  # missing values to 0!
  select_na = is.na(x)
  x[select_na] = 0.0
  
  # Just to increase contrast between imputed and observed values, let's
  # scale the observed values by a factor (not sure if this is really needed)
  if (should_scale_observed) {
    x[!select_na] = 1000.0 * x[!select_na]
  }
  return(x)
}

# Generate random forest model
# @param oc - Logical vector with TRUE = open, FALSE = close.
# @param elec - List with members:
#   $design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples
#   $meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSE - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#     $MED.OPEN - Median among observed open samples at each interaction pair, NA if no observed samples
#     $MED.CLOSE - Median among observed closed samples at each interaction pair, NA if no observed samples
#     $MIN.OPEN - Minimum value among observed open samples at each interaction pair, NA if no observed samples
#     $Q1.OPEN - 25% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MEA.OPEN - Mean value among observed open samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.OPEN - 75% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MAX.OPEN - Maximum value among observed open samples at each interaction pair, NA if no observed samples
#     $MIN.CLOSE - Minimum value among observed closed samples at each interaction pair, NA if no observed samples
#     $Q1.CLOSE - 25% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MEA.CLOSE - Mean value among observed closed samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.CLOSE - 75% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MAX.CLOSE - Maximum value among observed closed samples at each interaction pair, NA if no observed samples
#     $FAVORS - At each interaction pair, values are: (note: NA -> not observed)
#       "OPEN" if (MED.OPEN < MED.CLOSE) or (CLOSE is NA and MED.OPEN < 0) or (OPEN is NA and MED.CLOSE > 0)
#       "CLOSE" if (MED.CLOSE < MED.OPEN) or (OPEN is NA and MED.CLOSE < 0) or (CLOSE is NA and MED.OPEN > 0)
#       "NONE" o/w
# @param rdata_path - Path to write rdata cache, will load from this file if it exists unless
#   overwrite is TRUE
# @param keep - Specifies the interaction filter settings (only first element is used)
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
#     ALL.ION.POLAR - all ion-ion, ion-polar, polar-polar interactions are kept
# @param num_trees - The number of trees in the forest
# @param rf_version - Version of forest data
# @param should_scale_observed - If true, observed values will be scaled by a 'large' factor
#     to possibly help differentiate from non-observed (NA) values
# @param overwrite - If true, will always recompute random forest model and overwrite
#   anything already stored to disk. Else, will check for previously computed rdata and
#   return the cached result.
# @return random forest object
get_rf <- function(oc,
                   elec,
                   rdata_path,
                   keep = RF_INT_KEEP,
                   num_trees = NUM_RF_TREES,
                   rf_version = RF_INT_VERSION,
                   should_scale_observed = RF_SHOULD_SCALE_OBSERVED,
                   overwrite = SHOULD_OVERWRITE) {
  # Make sure labels and samples match
  stopifnot(all(names(oc) == rownames(elec$design)))
  # Make sure number of trees is reasonable
  stopifnot(num_trees > 50)
  
  # Check if cached model exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping model build as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  print(paste0("Building random forest model for: ", rdata_path))
  
  # Construct input data.frame
  X = get_rf_design(
    elec = elec,
    keep = keep,
    rf_version = RF_INT_VERSION,
    should_scale_observed = should_scale_observed
  )
  
  # Convert label to factor
  oc = as.factor(ifelse(oc, "open", "close"))
  
  # Fit a random forest
  rf.fit = randomForest(
    x = X,
    y = oc,
    ntree = num_trees,
    importance = TRUE,
    proximity = FALSE,
    keep.inbag = TRUE
  )
  
  save_rdata(rf.fit, rdata_path)
  return (rf.fit)
}

##############################
# Random forest - probability
#   correlations

library(kknn)

# @return Corresponding rf probability correlation path from random forest path
get_rf_prob_cor_path_from_rf <- function(rf_rdata_path) {
  return(insert_attrib_in_rdata(rf_rdata_path, "pocor"))
}

# @return Corresponding knn model fits path for rf correlation data
get_rf_knn_from_cor <- function(rf_cor_rdata_path) {
  return(insert_attrib_in_rdata(rf_cor_rdata_path, "knn"))
}

# @param rf.fit - Trained random forest model
# @param elec - List with members $design and $meta, see load_data_utils.R:
#   elec_csv_2_rdata() for details
# @param rdata_path - Output data path
# @param keep - Specifies the interaction filter settings (only first element is used)
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
# @param overwrite - If true, will always overwrite anything already stored to disk.
#   Else, will check for previously computed rdata and return the cached result.
get_rf_prob_cor <- function(rf.fit,
                            elec,
                            rdata_path,
                            keep = RF_INT_KEEP,
                            rf_version = RF_INT_VERSION,
                            overwrite = SHOULD_OVERWRITE) {
  # Check if cached model exists
  if (!overwrite && file.exists(rdata_path)) {
    print(
      paste0(
        "Skipping random forest probability correlation as rdata already exists: ",
        rdata_path
      )
    )
    return(load_rdata(rdata_path))
  }
  
  print(paste0("Analysing random forest probabilities for: ", rdata_path))
  
  # Capture out-of-bag probabilities for training input
  prob = predict(object = rf.fit, type = "prob")
  
  # Extract open probabilities
  prob_open = prob[, "open"]
  
  # Construct input data.frame
  X = get_rf_design(
    elec = elec,
    keep = keep,
    rf_version = rf_version,
    should_scale_observed = FALSE
  )
  
  # Make sure design and predictions are in same order
  # (They should be already, but just being extra careful)
  X = X[rownames(prob),]
  
  # Extract meta-information for each feature
  meta = elec$meta[colnames(X),]
  
  # Compute correlations with the open probability
  # NOTE: The correlation with the close probability will simply be the negative
  # of the open probability correlation!
  # @TODO - consider (-) log transform of probability since energies are
  #   additive in log probability space. Will need to set a min probability > 0,
  #   else log(p=0) = -Inf. However, a Spearman correlation should be able to
  #   detect a non-linear (exponential) trend, so this may not be necessary.
  pears = cor(X, prob_open, method = "pearson")[, 1] # coerce to vector!
  spear = cor(X, prob_open, method = "spearman")[, 1]
  stopifnot(all(names(pears) == names(spear)))
  
  # Compute GOF-statistic from forest floor, the idea is to also be
  # be able to catch interactions that are non-monotic
  n_features = ncol(X)
  gof = rep(NA, n_features)
  names(gof) = colnames(X)
  
  # Setup KNN params
  n_obs = nrow(X)
  knn_kmax = round(sqrt(n_obs) / 2)
  knn_kernel = "gaussian"
  knn_data = data.frame(p = prob_open, x = rep(NA, n_obs))
  
  # Iterate over each feature
  knn_out = list() # stores (x,y) vectors for future plotting
  for (i_feat in 1:n_features) {
    print(paste0("GOF: Fitting KNN model for feature ", i_feat, " of ", n_features))
    knn_data$x = X[, i_feat]
    # Train KNN using leave-one-out cross-validation
    knn_model = train.kknn(
      formula = p ~ .,
      data = knn_data,
      kmax = knn_kmax,
      kernel = knn_kernel
    )
    knn_fitted = knn_model$fitted.values
    # Extract the cross-validated predictions
    knn_fitted = knn_fitted[[length(knn_fitted)]]
    # Compute R^2
    gof[i_feat] = cor(x = prob_open, y = knn_fitted)
    gof[i_feat] = gof[i_feat] * gof[i_feat]
    # Store results for easier future plotting
    feat_name = colnames(X)[i_feat]
    knn_out[[feat_name]] = list(
      fit = data.frame(
        x.feat = knn_data$x,
        y.fit = knn_fitted,
        # This is inefficient but will make code to plot
        # much more similar to feature contributions code
        y.obs = prob_open
      ),
      model = knn_model,
      gof = gof[i_feat],
      spear = spear[feat_name],
      pears = pears[feat_name]
    )
  }
  
  # Extra paranoid, make sure everything is in same order
  meta = meta[names(pears),]
  gof = gof[names(pears)]
  stopifnot(all(names(pears) == rownames(meta)))
  results = data.frame(pears = pears,
                       spear = spear,
                       gof = gof,
                       meta)
  # Save correlations
  save_rdata(results, rdata_path)
  write.csv(x = results, file = get_csv_path_from_rdata(rdata_path))
  
  # Save knn data for possible future plotting.
  knn_path = get_rf_knn_from_cor(rdata_path)
  save_rdata(knn_out, knn_path)
  
  return(results)
}

##############################
# Forest floor - Visualize
#   feature contributions

library(forestFloor)

# @return Corresponding forest floor path from random forest path
get_ff_path_from_rf <- function(rf_rdata_path) {
  return(insert_attrib_in_rdata(rf_rdata_path, "ff"))
}

# @return Corresponding forrest floor correlation path from forest floor path
get_ff_cor_path_from_ff <- function(ff_rdata_path) {
  return(insert_attrib_in_rdata(ff_rdata_path, "cor"))
}

# Generate forest floor object which can be used for
# feature contribution visualization
# http://stats.stackexchange.com/questions/21152/obtaining-knowledge-from-a-random-forest
# http://forestfloor.dk/
# @param rf.fit - A trained random forest model
# @param rdata_path - Output path for forest floor object
# @param elec - List with members $design and $meta, see load_data_utils.R:
#   elec_csv_2_rdata() for details
# @param keep - specifies the interaction filter settings (only first element is used)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param should_scale_observed - If true, observed values will be scaled by a
#   'large' factor to possibly help differentiate from non-observed (NA) values
# @param overwrite - If true, will always recompute random forest model and overwrite
#   anything already stored to disk. Else, will check for previously computed rdata and
#   return the cached result.
get_ff <- function(rf.fit,
                   rdata_path,
                   elec,
                   keep = RF_INT_KEEP,
                   rf_version = RF_INT_VERSION,
                   should_scale_observed = RF_SHOULD_SCALE_OBSERVED,
                   overwrite = SHOULD_OVERWRITE) {
  # Check if cached object exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping forest floor build as rdata already exists: ",
      rdata_path
    ))
    return(load_rdata(rdata_path))
  }
  
  print(paste0("Generating forest floor for: ", rdata_path))
  
  # Get design data.frame
  X = get_rf_design(
    elec = elec,
    keep = keep,
    rf_version = rf_version,
    should_scale_observed = should_scale_observed
  )
  
  # Generate forest floor object
  ff = forestFloor(rf.fit = rf.fit,
                   X = X,
                   # Treat 2-class as regression
                   binary_reg = TRUE)
  
  # Cache for reuse
  save_rdata(ff, rdata_path)
  return(ff)
}

# @param ff - Forest floor object
# @param elec - List with members $design and $meta, see load_data_utils.R:
#   elec_csv_2_rdata() for details
# @param oc - logical with TRUE for OPEN state, FALSE for CLOSE state
# @param rdata_path - Output data path
# @param keep - Specifies the interaction filter settings (only first element is used)
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
#     ALL.ION.POL - all ion-ion, ion-polar, and polar-polar interactions are kept
# @param overwrite - If true, will always overwrite anything already stored to disk.
#   Else, will check for previously computed rdata and return the cached result.
get_ff_cor <- function(ff,
                       elec,
                       oc,
                       rdata_path,
                       keep = RF_INT_KEEP,
                       rf_version = RF_INT_VERSION,
                       overwrite = SHOULD_OVERWRITE) {
  # Check if cached model exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping forest floor correlation as rdata already exists: ",
      rdata_path
    ))
    return(load_rdata(rdata_path))
  }
  
  print(paste0("Analysing forest floor correlations for: ", rdata_path))
  
  # Construct input data.frame
  X = get_rf_design(
    elec = elec,
    keep = keep,
    rf_version = rf_version,
    should_scale_observed = FALSE
  )
  
  # Extract meta-information for each feature
  meta = elec$meta[colnames(X), ]
  
  # Alloc structures for storing correlations
  n_features = ncol(X)
  na_vec = rep(NA, n_features)
  names(na_vec) = colnames(X)
  pears = na_vec
  spear = na_vec
  gof = na_vec
  
  sum_all = data.frame(
    fc.med = na_vec,
    fc.mea = na_vec,
    fc.min = na_vec,
    fc.q1 = na_vec,
    fc.q3 = na_vec,
    fc.max = na_vec,
    fc.var = na_vec
  )
  
  sum_open = sum_all
  sum_close = sum_all
  
  med_diff = na_vec
  mea_diff = na_vec
  
  stopifnot(colnames(X) == rownames(pears))
  knn_out = list() # stores (x,y) vectors for future plotting
  
  # Setup KNN params
  n_obs = nrow(X)
  knn_kmax = round(sqrt(n_obs) / 2)
  knn_kernel = "gaussian"
  knn_data = data.frame(fc = rep(NA, n_obs), x = rep(NA, n_obs))
  
  # Summary utility
  do_summary <- function(sum_df, fc, i_feat) {
    # Summary offsets
    IX_MIN = 1
    IX_Q1 = 2
    IX_MED = 3
    IX_MEA = 4
    IX_Q3 = 5
    IX_MAX = 6
    # Parse summary
    s = summary(fc)
    sum_df$fc.min[i_feat] = suppressWarnings(as.numeric(s[IX_MIN]))
    sum_df$fc.q1[i_feat] = suppressWarnings(as.numeric(s[IX_Q1]))
    sum_df$fc.med[i_feat] = suppressWarnings(as.numeric(s[IX_MED]))
    sum_df$fc.mea[i_feat] = suppressWarnings(as.numeric(s[IX_MEA]))
    sum_df$fc.q3[i_feat] = suppressWarnings(as.numeric(s[IX_Q3]))
    sum_df$fc.max[i_feat] = suppressWarnings(as.numeric(s[IX_MAX]))
    sum_df$fc.var[i_feat] = var(fc)
    return(sum_df)
  }
  
  # Compute correlations for each feature
  stopifnot(rownames(ff$FCmatrix) == rownames(X))
  stopifnot(rownames(oc) == rownames(X))
  for (i_feat in 1:n_features) {
    x.feat = X[, i_feat]
    fc = ff$FCmatrix[, i_feat]
    pears[i_feat] = cor(x = x.feat,
                        y = fc,
                        method = c("pearson"))
    spear[i_feat] = cor(x = x.feat,
                        y = fc,
                        method = c("spearman"))
    
    print(paste0("GOF: Fitting KNN model for feature ", i_feat, " of ", n_features))
    knn_data$x = x.feat
    knn_data$fc = fc
    # Train KNN using leave-one-out cross-validation
    knn_model = train.kknn(
      formula = fc ~ .,
      data = knn_data,
      kmax = knn_kmax,
      kernel = knn_kernel
    )
    knn_fitted = knn_model$fitted.values
    # Extract the cross-validated predictions
    knn_fitted = knn_fitted[[length(knn_fitted)]]
    # Compute R^2
    gof[i_feat] = cor(x = fc, y = knn_fitted)
    gof[i_feat] = gof[i_feat] * gof[i_feat]
    # Store results for easier future plotting
    feat_name = colnames(X)[i_feat]
    knn_out[[feat_name]] = list(
      fit = data.frame(
        x.feat = knn_data$x,
        y.fit = knn_fitted,
        # This is inefficient but will make code to plot
        # much more similar to feature contributions code
        y.obs = fc
      ),
      model = knn_model,
      gof = gof[i_feat],
      spear = spear[i_feat],
      pears = pears[i_feat]
    )
    
    # Compute summaries
    sum_all = do_summary(sum_all, fc, i_feat)
    sum_open = do_summary(sum_open, fc[oc], i_feat)
    sum_close = do_summary(sum_close, fc[!oc], i_feat)
    # Note: a good feature contribution probably has
    # a large median difference and a high GOF
    med_diff[i_feat] = sum_open$fc.med[i_feat] - sum_close$fc.med[i_feat]
    mea_diff[i_feat] = sum_open$fc.mea[i_feat] - sum_close$fc.med[i_feat]
  }
  
  # Extra paranoid, make sure everything is in same order
  meta = meta[names(pears), ]
  gof = gof[names(pears)]
  stopifnot(all(names(pears) == rownames(meta)))
  
  colnames(sum_open) = paste0(colnames(sum_open), ".open")
  colnames(sum_close) = paste0(colnames(sum_close), ".close")
  
  results = data.frame(
    pears = pears,
    pears2 = pears * pears,
    spear = spear,
    spear2 = spear * spear,
    gof = gof,
    med_diff = med_diff,
    abs_med_diff = abs(med_diff),
    mea_diff = mea_diff,
    abs_mea_diff = abs(mea_diff),
    sum_all,
    sum_open,
    sum_close,
    meta
  )
  
  # Save correlations
  save_rdata(results, rdata_path)
  write.csv(x = results, file = get_csv_path_from_rdata(rdata_path))
  
  # Save knn data for possible future plotting.
  knn_path = get_rf_knn_from_cor(rdata_path)
  save_rdata(knn_out, knn_path)
  
  return(results)
}

##############################
# Random forest plot utils

# Utility to obtain a human readable interaction name
# Converts <resno>.<resno> to <res1LetterCode><resno>.<res1LetterCode><resno>
# Eg. "68.228" becomes "R68.R228"
get_human_int_feature_name <- function(df) {
  resa.code = AA3to1[df$RESA.NAME]
  resa.no = df$RESA.NO
  resb.code = AA3to1[df$RESB.NAME]
  resb.no = df$RESB.NO
  feat.names = paste0(resa.code, resa.no, ".", resb.code, resb.no)
  names(feat.names) = rownames(df)
  return(feat.names)
}

# Generate a ggplot2 scatter plot with fit line and color coded by oc-status
# @param feat.ls - List with attributes for a single feature
#   $fit = data.frame with columns x.feat, y.fit, and y.obs where:
#     $x.feat - Observed input feature value at each sample (independent var)
#     $y.fit - KNN fitted (predicted) output feature mapping at each sample (dependent var)
#     $y.obs - Observed output feature mapping at each sample (dependent var)
#   $model - The actual trained KNN model used to generate fit$y.fit from input fit$x.feat
#   $gof - KNN goodness of fit statistic defined as squared correlation(fit$y.fit, fit$y.obs)
#     (See forest floor paper)
#   $spear - Spearman correlation(fit$x.feat, fit$y.obs)
#   $pears - Pearson correlation(fit$x.feat, fit$y.obs)
# @param oc - Logical vector with TRUE for open samples, FALSE for close samples
# @param main - The title of the plot
# @param xlab - Label for x-axis
# @param ylab - Label for y-axis
plot_ff_scatter <- function(feat.ls,
                            oc,
                            main,
                            xlab,
                            ylab) {
  # http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html
  # http://stackoverflow.com/questions/19921842/plotting-multiple-time-series-on-the-same-plot-using-ggplot
  library(ggplot2)
  
  fit.df = feat.ls$fit
  oc = oc[rownames(fit.df)]
  
  # Creat plotting data.frame
  fit.df$State = as.factor(ifelse(oc, "Open", "Close"))
  
  p <- ggplot() +
    geom_jitter(data = fit.df, aes(x = x.feat, y = y.obs, color = State)) +
    geom_line(data = fit.df, aes(x = x.feat, y = y.fit), size = 0.1) +
    labs(x = xlab, y = ylab, title = main)
  return(p)
}

# http://stackoverflow.com/questions/4559229/drawing-pyramid-plot-using-r-and-ggplot2
# https://learnr.wordpress.com/2009/09/24/ggplot2-back-to-back-bar-charts/
# Generates bar plots of the median feature contribution differences
# @param ff_cor - Forest floor correlation data.frame
# @param main - Title of plot
# @param xlab - X-axis label
# @param ylab - Y-axis label
# @param topk - Only top-k features are plotted
plot_ff_bar <- function(ff_cor,
                        feat_name_fnc = get_human_int_feature_name,
                        main = "Random Forest Most Influential Interactions",
                        xlab = "OmpG Residue Interactions",
                        ylab = "Influence Score",
                        topk = min(35, nrow(ff_cor))) {
  stopifnot(topk > 1)
  library(ggplot2)
  # Order by decreasing median feature contribution
  ix_order = order(ff_cor$abs_med_diff, decreasing = TRUE)
  ff_cor = ff_cor[ix_order,]
  ff_cor = ff_cor[1:topk,]
  feat.names = feat_name_fnc(ff_cor)
  ff_cor$feature = feat.names # as.factor(feat.names)
  p = ggplot(ff_cor) +
    geom_bar(aes(
      x = reorder(feature, abs_med_diff, function(x) {
        -x
      }),
      y = abs_med_diff,
      fill = abs_med_diff
    ),
    stat = "identity") +
    labs(x = xlab, y = ylab, title = main) +
    theme(axis.text.x = element_text(angle = -90),
          legend.position = "none")
  return(p)
}

# Plots random forest features
plot_rf_importance <- function(rf.fit,
                               meta,
                               feat_name_fnc = get_human_int_feature_name,
                               main = "Random Forest Variable Importance",
                               topk = min(30, nrow(rf.fit$importance))) {
  feat.names = feat_name_fnc(meta)
  rownames(rf.fit$importance) = feat.names[rownames(rf.fit$importance)]
  varImpPlot(
    x = rf.fit,
    scale = FALSE,
    main = main,
    n.var = topk
  )
  return(NA)
}

##############################
# Random forest - Notes

# Generate importance plots
# See https://dinsdalelab.sdsu.edu/metag.stats/code/randomforest.html
# In particular, see method varImpPlot in randomForest package:
#
# varImpPlot(WT_PH5_RF_ALL_ION_RCB)
# varImpPlot(WT_PH7_RF_ALL_ION_RCB)
#
# ------------------
#
# Plot MSE as a function of the number of trees used:
# plot(WT_PH5_RF_ALL_ION_RCB)
# plot(WT_PH7_RF_ALL_ION_RCB)
#
# ------------------
#
# Note, just typing "WT_PH5_RF_ALL_ION_RCB" into command terminal will give OOB error
# which is the misclassification error.
# http://stackoverflow.com/questions/28666066/get-the-accuracy-of-a-random-forest-in-r

##############################
# Random forest - Examples

# Example wild type random forest analysis
# @param keep.ls - specifies the interaction filter settings (supports multiple elements)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param pH.ls - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_rf_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                         pH.ls = c(5, 7)) {
  results = list()
  
  oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if ((!(pH_str %in% names(oc_by_pH))) ||
        (!(pH_str %in% names(elec_by_pH)))) {
      next
    }
    
    oc = oc_by_pH[[pH_str]]
    elec = elec_by_pH[[pH_str]]
    
    for (keep in keep.ls) {
      res_id = paste0("WT_PH", pH, "_RF_", keep, "_RCB")
      rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep
      )
      results[[res_id]] = get_rf(
        oc = oc,
        elec = elec,
        rdata_path = rdata_path,
        keep = keep
      )
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  # Return results!
  return(results)
}

# Example wild type probability correlation analysis
# @param keep - specifies the interaction filter settings (supports multiple elements)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param pH - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_rf_prob_cor_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                                  pH.ls = c(5, 7)) {
  rf.fits = do_rf_wt_rcb(keep.ls, pH.ls)
  
  results = list()
  
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if (!(pH_str %in% names(elec_by_pH))) {
      next
    }
    
    elec = elec_by_pH[[pH_str]]
    
    for (keep in keep.ls) {
      rfpcor_res_id = paste0("WT_PH", pH, "_RFPCOR_", keep, "_RCB")
      rf_res_id = paste0("WT_PH", pH, "_RF_", keep, "_RCB")
      
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep
      )
      
      results[[rfpcor_res_id]] = get_rf_prob_cor(
        rf.fit = rf.fits[[rf_res_id]],
        rdata_path = get_rf_prob_cor_path_from_rf(rf_rdata_path),
        elec = elec,
        keep = keep
      )
      
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Example wild type forest floor analysis
# @param keep - specifies the interaction filter settings (supports multiple elements)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param pH - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_ff_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                         pH.ls = c(5, 7)) {
  rf.fits = do_rf_wt_rcb(keep.ls, pH.ls)
  
  results = list()
  
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if (!(pH_str %in% names(elec_by_pH))) {
      next
    }
    
    elec = elec_by_pH[[pH_str]]
    
    for (keep in keep.ls) {
      ff_res_id = paste0("WT_PH", pH, "_FF_", keep, "_RCB")
      rf_res_id = paste0("WT_PH", pH, "_RF_", keep, "_RCB")
      
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep
      )
      
      results[[ff_res_id]] = get_ff(
        rf.fit = rf.fits[[rf_res_id]],
        rdata_path = get_ff_path_from_rf(rf_rdata_path),
        elec = elec,
        keep = keep
      )
      
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Example wild type feature contribution correlation analysis
# @param keep - specifies the interaction filter settings (supports multiple elements)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param pH - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_ff_cor_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                             pH.ls = c(5, 7)) {
  ff_ls = do_ff_wt_rcb(keep.ls, pH.ls)
  
  results = list()
  
  oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if (!(pH_str %in% names(elec_by_pH))) {
      next
    }
    
    oc = oc_by_pH[[pH_str]]
    elec = elec_by_pH[[pH_str]]
    
    for (keep in keep.ls) {
      ff_cor_res_id = paste0("WT_PH", pH, "_FFCOR_", keep, "_RCB")
      ff_res_id = paste0("WT_PH", pH, "_FF_", keep, "_RCB")
      
      # Determine output path
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep
      )
      ff_rdata_path = get_ff_path_from_rf(rf_rdata_path)
      ff_cor_rdata_path = get_ff_cor_path_from_ff(ff_rdata_path)
      
      results[[ff_cor_res_id]] = list(
        ff_cor = get_ff_cor(
          ff = ff_ls[[ff_res_id]],
          rdata_path = ff_cor_rdata_path,
          elec = elec,
          oc = oc,
          keep = keep
        ),
        ff_cor_rdata_path = ff_cor_rdata_path
      )
      
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Example wild type feature contribution scatter plots
# @param keep - specifies the interaction filter settings (supports multiple elements)
#   ALL.ION - all ion-ion interactions are kept
#   LOOP.ION - only ion-ion interactions between loop regions are kept
#   XORLOOP.ION - only ion-ion -interactions from loop to non-loop regions are kept
# @param pH - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_ff_scatter_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                                 pH.ls = c(5, 7)) {
  ff_cor_ls = do_ff_cor_wt_rcb(keep.ls, pH.ls)
  
  oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if ((!(pH_str %in% names(oc_by_pH))) ||
        (!(pH_str %in% names(elec_by_pH)))) {
      next
    }
    
    # Create image subdirectory
    img_dir = file.path(get_relax_c_beta_dir(),
                        "capt",
                        paste0("pH", pH),
                        "wt",
                        "plot_ff_scatter")
    dir.create(img_dir, showWarnings = FALSE)
    
    oc = oc_by_pH[[pH_str]]
    meta = elec_by_pH[[pH_str]]$meta
    
    # Obtain mapping to human readable name
    human.feat.names = get_human_int_feature_name(meta)
    
    for (keep in keep.ls) {
      ff_cor_res_id = paste0("WT_PH", pH, "_FFCOR_", keep, "_RCB")
      ff_cor_res = ff_cor_ls[[ff_cor_res_id]]
      knn_path = get_rf_knn_from_cor(ff_cor_res$ff_cor_rdata_path)
      knn_ls = load_rdata(knn_path)
      
      for (feat.name in names(knn_ls)) {
        knn_data = knn_ls[[feat.name]]
        human.feat.name = human.feat.names[feat.name]
        
        png_name = paste0(feat.name, ".", keep, ".fccor.scatter.png")
        png_path = file.path(img_dir, png_name)
        print(paste0("Generating scatter plot for: ", png_path))
        main = paste0(human.feat.name, " Random Forest Feature Contribution")
        main = paste0(main,
                      "\nGOF=",
                      knn_data$gof,
                      ", Spear^2=",
                      knn_data$spear,
                      ", pH=",
                      pH_str)
        xlab = paste0(human.feat.name, " Electrostatic Potential Energy")
        ylab = expression(paste(Delta, " Open Probability"))
        p = plot_ff_scatter(
          feat.ls = knn_data,
          oc = oc,
          main = main,
          xlab = xlab,
          ylab = ylab
        )
        
        ggsave(filename = png_path, plot = p)
      } # end iteration over features (interaction pairs)
    } # end iteration over keep
  } # end iteration over pH
}

# Generate feature contribution bar plots
do_ff_bar_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                             pH.ls = c(5, 7)) {
  ff_cor_ls = do_ff_cor_wt_rcb(keep.ls, pH.ls)
  
  for (pH in pH.ls) {
    # Create image subdirectory
    img_dir = file.path(get_relax_c_beta_dir(),
                        "capt",
                        paste0("pH", pH),
                        "wt",
                        "plot_ff_bar")
    dir.create(img_dir, showWarnings = FALSE)
    
    for (keep in keep.ls) {
      ff_cor_res_id = paste0("WT_PH", pH, "_FFCOR_", keep, "_RCB")
      # Skip missing pH levels
      if (!(ff_cor_res_id %in% names(ff_cor_ls))) {
        print(paste0("Identifier: ", ff_cor_res_id, " not found."))
        next
      }
      
      ff_cor = ff_cor_ls[[ff_cor_res_id]]
      
      main = paste0("Random Forest Most Influential Interactions, pH=", pH)
      numvar = nrow(ff_cor)
      topk = min(35, numvar)
      main = paste0(main, "\n(Top ", topk, " of ", numvar, ")")
      p = plot_ff_bar(ff_cor = ff_cor,
                      main = main,
                      topk = topk)
      png_name = paste0("ff_bar.pH", pH, ".", keep, ".png")
      png_path = file.path(img_dir, png_name)
      print(paste0("Generating: ", png_path))
      ggsave(filename = png_path, plot = p)
    } # end iteration over keep
  } # end iteration over pH
}

# Calls varImpPlot for each setting
do_rf_importance_wt_rcb <- function(keep.ls = RF_INT_KEEP,
                                    pH.ls = c(5, 7)) {
  rf.fits = do_rf_wt_rcb(keep.ls, pH.ls)
  
  results = list()
  
  elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if (!(pH_str %in% names(elec_by_pH))) {
      next
    }
    
    # Create image subdirectory
    img_dir = file.path(get_relax_c_beta_dir(),
                        "capt",
                        paste0("pH", pH),
                        "wt",
                        "plot_rf_imp")
    dir.create(img_dir, showWarnings = FALSE)
    
    elec = elec_by_pH[[pH_str]]
    meta = elec$meta
    
    for (keep in keep.ls) {
      rf_res_id = paste0("WT_PH", pH, "_RF_", keep, "_RCB")
      rf.fit = rf.fits[[rf_res_id]]
      
      main = paste0("Random Forest Variable Importance, pH=", pH)
      numvar = nrow(rf.fit$importance)
      topk = min(30, numvar)
      main = paste0(main, "\n(Top ", topk, " of ", numvar, ")")
      
      png_name = paste0("rf_imp.pH", pH, ".", keep, ".png")
      png_path = file.path(img_dir, png_name)
      
      print(paste0("Generating: ", png_path))
      png(filename = png_path,
          width = 1024,
          height = 1024)
      plot_rf_importance(
        rf.fit = rf.fit,
        meta = meta,
        main = main,
        topk = topk
      )
      dev.off()
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

###################################################################
# Boxplot interaction
###################################################################

# Generates a comparative box plot of the interaction strength of a specific
# electrostatic interaction pair. Comparisons are made between the open and
# and closed ensembles. The interaction pair is specified by reside
# sequence numbers: <res_num_a>.<res_num_b>
#
# E.g:
# boxplot_interaction(WT_PH5_ELEC_RCB$design, "97.227", WT_PH5_OC_RCB, WT_PH5_ELEC_RCB$meta, "wt", 5, 1.5, FALSE)
#
# @param design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples
# @param int_name - The interaction pair name in format "<res_num_a>.<res_num_b>"
# @param meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSe - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#     $MED.OPEN - Median among observed open samples at each interaction pair, NA if no observed samples
#     $MED.CLOSE - Median among observed closed samples at each interaction pair, NA if no observed samples
#     $MIN.OPEN - Minimum value among observed open samples at each interaction pair, NA if no observed samples
#     $Q1.OPEN - 25% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MEA.OPEN - Mean value among observed open samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.OPEN - 75% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MAX.OPEN - Maximum value among observed open samples at each interaction pair, NA if no observed samples
#     $MIN.CLOSE - Minimum value among observed closed samples at each interaction pair, NA if no observed samples
#     $Q1.CLOSE - 25% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MEA.CLOSE - Mean value among observed closed samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.CLOSE - 75% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MAX.CLOSE - Maximum value among observed closed samples at each interaction pair, NA if no observed samples
#     $FAVORS - At each interaction pair, values are: (note: NA -> not observed)
#       "OPEN" if (MED.OPEN < MED.CLOSE) or (CLOSE is NA and MED.OPEN < 0) or (OPEN is NA and MED.CLOSE > 0)
#       "CLOSE" if (MED.CLOSE < MED.OPEN) or (OPEN is NA and MED.CLOSE < 0) or (CLOSE is NA and MED.OPEN > 0)
#       "NONE" o/w
# @param sim - String with the simulation identifier (used for plot title)
# @param pH - The pH value for the captured interactions (used for plot title)
# @param range - Parameter to boxplot(...), determines how far the plot whiskers extend
#   out from the box
# @param outline - Parameter to boxplot(...)
#   If TRUE: outliers are rendered; If FALSE: no outliers are rendered
# @param filt.na - If TRUE: NA's are removed prior to plotting.
#   If FALSE: NAs are replaced with 0's.
boxplot_interaction <-
  function(design,
           int_name,
           oc,
           meta,
           sim = "wt",
           pH = 5,
           range = 1.5,
           outline = FALSE,
           filt.na = TRUE) {
    d.int = design[, int_name]
    
    # Filter NA's
    if (filt.na) {
      na_filt = !is.na(d.int)
      d.int = d.int[na_filt]
      oc = oc[na_filt]
    }
    else {
      # Replace with zeros and re-scale
      # This is what random forest sees
      d.int[is.na(d.int)] = 0
      d.int = 1000.0 * d.int
    }
    
    # Create title
    nfo = meta[int_name, ]
    main = paste0(
      nfo$RESA.NAME,
      "-",
      nfo$RESA.NO,
      "-",
      nfo$RESA.LOOP,
      " : ",
      nfo$RESB.NAME,
      "-",
      nfo$RESB.NO,
      "-",
      nfo$RESB.LOOP,
      " @ ",
      sim,
      ", pH=",
      pH
    )
    
    boxplot(
      x = list(d.int[oc], d.int[!oc]),
      names = c(
        paste0("open,\nn=", nfo$N.OPEN),
        paste0("close,\nn=", nfo$N.CLOSE)
      ),
      range = range,
      outline = outline,
      main = main
    )
    return(nfo)
  }

###################################################################
# Backbone analysis
###################################################################

# Prints percentage of unique backbones in parameter data sets
# Examples:
# common_backbone_summary(rownames(WT_PH5_ELEC_RCB$design), rownames(WT_PH7_ELEC_RCB$design), WT_PH5_OC_RCB, WT_PH7_OC_RCB)
# common_backbone_summary(rownames(WT_PH5_ELEC_RCB$design)[1:2500], rownames(WT_PH7_ELEC$design)[1:2500], WT_PH5_OC[1:2500], WT_PH7_OC[1:2500])
common_backbone_summary <-
  function(names1, names2, oc1, oc2, rtrim = 19) {
    base.names1 = substr(x = names1, 1, nchar(names1) - rtrim)
    base.names2 = substr(x = names2, 1, nchar(names2) - rtrim)
    common.names = intersect(base.names1, base.names2)
    
    perc_unique1 = 1.0 - (length(names1[base.names1 %in% common.names]) / length(names1))
    perc_unique2 = 1.0 - (length(names2[base.names2 %in% common.names]) / length(names2))
    print(paste0("Data set 1, % unique backbones: ", 100.0 * perc_unique1))
    print(paste0("Data set 2, % unique backbones: ", 100.0 * perc_unique2))
    
    rank1 = rep(1, length(names2))
    rank1 = cumsum(rank1)
    names(rank1) = base.names1
    rank2 = rep(1, length(names2))
    rank2 = cumsum(rank2)
    names(rank2) = base.names2
    
    com.rank1 = rank1[common.names]
    com.rank2 = rank2[common.names]
    
    signed_rank_delta = com.rank1 - com.rank2
    print("Signed energy rank delta:")
    print(summary(signed_rank_delta))
    
    abs_rank_delta = abs(com.rank1 - com.rank2)
    print("Absolute energy rank delta:")
    print(summary(abs_rank_delta))
    
    # OC analysis - are there any that flipped?
    base.oc1 = oc1
    names(base.oc1) = base.names1
    base.oc2 = oc2
    names(base.oc2) = base.names2
    
    com.oc1 = base.oc1[common.names]
    com.oc2 = base.oc2[common.names]
    flip = com.oc1 != com.oc2
    
    n.flip = sum(flip)
    n.flip.open = sum(com.oc2[flip])
    n.flip.close = n.flip - n.flip.open
    perc.flip = n.flip / length(base.oc1)
    prop.flip.open = n.flip.open / n.flip
    prop.flip.close = n.flip.close / n.flip
    print(paste0(
      "% samples that switch open/close state from set 1 -> 2: ",
      100.0 * perc.flip
    ))
    print(paste0("proportion that flip open: ", prop.flip.open))
    print(paste0("proportion that flip close: ", prop.flip.close))
  }

###################################################################
# Mann-Whitney
###################################################################

# @return name of Mann-Whitney rdata file
get_mann_whitney_rdata_name <- function(base.name,
                                        replace.na = FALSE,
                                        replace.val = 0.0,
                                        keep = c("ALL.ION", "LOOP.ION", "XORLOOP.ION"),
                                        max.energy.rank = MAX_ENERGY_RANK) {
  return(
    paste0(
      base.name,
      ".nar",
      sum(replace.na),
      ".nav",
      replace.val,
      ".mw.r",
      max.energy.rank,
      paste0(".", keep[1]),
      ".rdata"
    )
  )
}

# @return Path to Mann-Whitney rdata
get_mann_whitney_path <-
  function(base_output_dir = get_relax_c_beta_dir(),
           pH = 5,
           sim_id = "wt",
           replace.na = FALSE,
           replace.val = 0.0,
           keep = c("ALL.ION", "LOOP.ION", "XORLOOP.ION"),
           max_energy_rank = MAX_ENERGY_RANK) {
    pH_str = paste0("pH", pH)
    base.name = paste0(sim_id, ".", pH_str)
    file.path(
      base_output_dir,
      "capt",
      pH_str,
      sim_id,
      get_mann_whitney_rdata_name(base.name,
                                  replace.na,
                                  replace.val,
                                  keep,
                                  max_energy_rank)
    )
  }

##############################
# Mann-Whitney - ignore NAs

# WT RCB

WT_PH5_MW_IGNOR_NA_ALL_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id = "wt",
  replace.na =
    FALSE,
  keep = "ALL.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IGNOR_NA_ALL_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id = "wt",
  replace.na =
    FALSE,
  keep = "ALL.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH5_MW_IGNOR_NA_LOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id = "wt",
  replace.na =
    FALSE,
  keep = "LOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IGNOR_NA_LOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id = "wt",
  replace.na =
    FALSE,
  keep = "LOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH5_MW_IGNOR_NA_XORLOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id =
    "wt",
  replace.na =
    FALSE,
  keep = "XORLOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IGNOR_NA_XORLOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id =
    "wt",
  replace.na =
    FALSE,
  keep = "XORLOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)
##############################
# Mann-Whitney - impute NAs

# WT RCB

WT_PH5_MW_IMPUT_NA_ALL_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id = "wt",
  replace.na =
    TRUE,
  keep = "ALL.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IMPUT_NA_ALL_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id = "wt",
  replace.na =
    TRUE,
  keep = "ALL.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH5_MW_IMPUT_NA_LOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id = "wt",
  replace.na =
    TRUE,
  keep = "LOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IMPUT_NA_LOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id = "wt",
  replace.na =
    TRUE,
  keep = "LOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH5_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id =
    "wt",
  replace.na =
    TRUE,
  keep = "XORLOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

WT_PH7_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB = get_mann_whitney_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id =
    "wt",
  replace.na =
    TRUE,
  keep = "XORLOOP.ION",
  max_energy_rank =
    MAX_ENERGY_RANK
)

# @param design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples.
# @param oc - logical vector with TRUE = open, FALSE = close.
# @param meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSe - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#   See description in load_data_utils::elec_csv_2_rdata() for full column listings
# @param rdata_path - Output path for Mann-Whitney rdata
# @param replace.na - If TRUE, NAs are replaced with replace.val, else they are removed
# @param replace.val - Only used if replace.na is TRUE - all NAs are replaced with this value
# @param keep - specifies the interaction filter settings (only first element is used)
#     ALL.ION - all ion-ion interactions are kept
#     LOOP.ION - only ion-ion interactions between loop regions are kept
#     XORLOOP.ION - only ion-ion interactions from loop to non-loop regions are kept
# @param overwrite - If true, will always recompute analysis and overwrite anything already
#   stored to disk. Else, will check for previously computed rdata and return the cached result.
# @return data.frame with the following columns
#     $PVAL - The p-value for the interaction pair
#     $FDR - The false discovery rate adjusted p-value
#     $BON - The Bonferroni adjusted p-value
#     ... all columns from parameter meta data.frame are also included
get_mann_whitney <-
  function(design,
           oc,
           meta,
           rdata_path,
           replace.na = FALSE,
           replace.val = 0.0,
           keep = c("ALL.ION", "LOOP.ION", "XORLOOP.ION"),
           overwrite = SHOULD_OVERWRITE) {
    print(paste0("Mann-Whitney testing using output path: ", rdata_path))
    # Check if cached model exists
    if (!overwrite && file.exists(rdata_path)) {
      print(paste0(
        "Skipping Mann-Whitney testing as rdata already exists: ",
        rdata_path
      ))
      return(load_rdata(rdata_path))
    }
    
    # Filter covariates so only ion-ion & loop-loop interactions are kept
    # Default to "ALL.ION" - all ion-ion interactions are kept and at least one residue is part of a loop
    filt = meta$ION.ION &
      ((meta$RESA.LOOP != "NONE") | (meta$RESB.LOOP != "NONE"))
    keep = keep[1]
    if (keep == "LOOP.ION") {
      # Only ion-ion interactions between loop regions are kept
      filt = meta$LOOP.ION
    } else if (keep == "XORLOOP.ION") {
      # only ion-ion interactions from loop to non-loop regions are kept
      filt = meta$ION.ION &
        xor((meta$RESA.LOOP != "NONE"), (meta$RESB.LOOP != "NONE"))
    }
    # Also, only keep inter-loop interactions
    filt = filt & (meta$RESA.LOOP != meta$RESB.LOOP)
    
    loop_ion_names = rownames(meta)[filt]
    design = design[, loop_ion_names]
    
    # Total number of interaction pairs
    n.ints = ncol(design)
    
    # Allocate results vector
    PVAL = rep(NA, n.ints)
    names(PVAL) = colnames(design)
    
    # Check if we need to replace missing values
    if (replace.na) {
      print(paste0("...Replacing NAs with ", replace.val))
      design[is.na(design)] = replace.val
    }
    
    # Pre-separate open and closed samples. Hopefully enough memory for this.
    # The idea is to avoid repeatedly filtering within the loop.
    design.open = design[oc, ]
    design.close = design[!oc, ]
    
    # Avoid testing samples with less than this number of samples
    # in both the open and closed classes.
    MIN_NUM_SAMPLES_PER_CLASS = 30
    
    print("...Testing interactions")
    for (i.int in 1:n.ints) {
      x.open = design.open[, i.int]
      y.close = design.close[, i.int]
      
      # Filter missing values
      x.open = x.open[!is.na(x.open)]
      y.close = y.close[!is.na(y.close)]
      
      # Check if we have enough values to proceed
      if (length(x.open) < MIN_NUM_SAMPLES_PER_CLASS) {
        next
      }
      if (length(y.close) < MIN_NUM_SAMPLES_PER_CLASS) {
        next
      }
      # Non-parametric test the two populations
      result = wilcox.test(
        x = x.open,
        y = y.close,
        alternative = "two.sided",
        exact = FALSE,
        paired = FALSE
      )
      # Store p-value
      PVAL[i.int] = result$p.value
    }
    
    # Correct for multiple hypothesis testing
    FDR = p.adjust(PVAL, method = "fdr")
    BON = p.adjust(PVAL, method = "bonferroni")
    
    # Store into data.frame
    results = data.frame(PVAL, FDR, BON, meta[filt, ])
    # Sort by false discovery rate
    results = results[order(FDR), ]
    # Save results
    save_rdata(results, rdata_path)
    csv_path = get_csv_path_from_rdata(rdata_path)
    print(paste0("...Saving to: ", csv_path))
    write.csv(x = results, file = csv_path)
    return (results)
  }

##############################
# Mann-Whitney - Examples

# Example wildtype analysis
do_mann_whitney_wt_rcb <- function() {
  # Iterate over keep-filter
  keep = c("ALL.ION", "LOOP.ION", "XORLOOP.ION")
  
  # Parallel arrays for output path based on keep-filter and pH, etc.
  keep_pH5_ignor_paths_rcb = c(
    WT_PH5_MW_IGNOR_NA_ALL_ION_PATH_RCB,
    WT_PH5_MW_IGNOR_NA_LOOP_ION_PATH_RCB,
    WT_PH5_MW_IGNOR_NA_XORLOOP_ION_PATH_RCB
  )
  
  keep_pH7_ignor_paths_rcb = c(
    WT_PH7_MW_IGNOR_NA_ALL_ION_PATH_RCB,
    WT_PH7_MW_IGNOR_NA_LOOP_ION_PATH_RCB,
    WT_PH7_MW_IGNOR_NA_XORLOOP_ION_PATH_RCB
  )
  
  keep_pH5_imput_paths_rcb = c(
    WT_PH5_MW_IMPUT_NA_ALL_ION_PATH_RCB,
    WT_PH5_MW_IMPUT_NA_LOOP_ION_PATH_RCB,
    WT_PH5_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB
  )
  
  keep_pH7_imput_paths_rcb = c(
    WT_PH7_MW_IMPUT_NA_ALL_ION_PATH_RCB,
    WT_PH7_MW_IMPUT_NA_LOOP_ION_PATH_RCB,
    WT_PH7_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB
  )
  
  for (i in 1:length(keep)) {
    # Ignore NAs
    
    wt_pH5_mw_ignor_na_rcb = get_mann_whitney(
      WT_PH5_ELEC_RCB$design,
      WT_PH5_OC_RCB,
      WT_PH5_ELEC_RCB$meta,
      keep_pH5_ignor_paths_rcb[i],
      keep = keep[i],
      overwrite = SHOULD_OVERWRITE
    )
    
    wt_pH7_mw_ignor_na_rcb = get_mann_whitney(
      WT_PH7_ELEC_RCB$design,
      WT_PH7_OC_RCB,
      WT_PH7_ELEC_RCB$meta,
      keep_pH7_ignor_paths_rcb[i],
      keep = keep[i],
      overwrite = SHOULD_OVERWRITE
    )
    
    # Impute NAs
    
    wt_pH5_mw_imput_na_rcb = get_mann_whitney(
      WT_PH5_ELEC_RCB$design,
      WT_PH5_OC_RCB,
      WT_PH5_ELEC_RCB$meta,
      keep_pH5_imput_paths_rcb[i],
      replace.na = TRUE,
      keep = keep[i],
      overwrite = SHOULD_OVERWRITE
    )
    
    wt_pH7_mw_imput_na_rcb = get_mann_whitney(
      WT_PH7_ELEC_RCB$design,
      WT_PH7_OC_RCB,
      WT_PH7_ELEC_RCB$meta,
      keep_pH7_imput_paths_rcb[i],
      replace.na = TRUE,
      keep = keep[i],
      overwrite = SHOULD_OVERWRITE
    )
  }
}

###################################################################
# Consensus variable importance analysis
###################################################################

# @return Path to random forest rdata
get_cons_var_imp_path <-
  function(base_output_dir = get_relax_c_beta_dir(),
           pH = 5,
           sim_id = "wt",
           keep = "ALL.ION",
           topk = 30) {
    rdata_path = file.path(
      base_output_dir,
      "capt",
      paste0("pH", pH),
      sim_id,
      paste0("var.imp.",
             keep[1],
             ".topk",
             topk,
             ".rdata")
    )
  }

# Obtains the consensus among top-k most important variables for
# distinguishing open and closed OmpG states using several different
# methods:
#   1. Random forest feature contribution (median difference)
#   2. Random forest mean decrease in accuracy
#   3. Random forest mean decrease in gini
#   4. Mann-whitney FDR ranking
# The variables in the intersection of all methods are reported
# along with their average rank (hence, reported number of
# consensus interactions will always be less than topk)
# @param ff_cor - feature contribution data.frame
# @param rf.fit - trained random forest model
# @param mw - mann-whitney data.frame - assumed to be sorted by FDR
# @param topk - The top-k interactions for each variable importance
#   metric are intersected against each other and then their
#   average rank is reported.
# @param rdata_path - Output path for Mann-Whitney rdata
# @param feat_name_fnc - Human interaction feature name
# @param overwrite - If true, will always recompute analysis and overwrite
#   anything already stored to disk. Else, will check for previously
#   computed rdata and return the cached result.
# @return data.frame with intersection set and mean rank
get_cons_var_imp <- function(ff_cor,
                             rf.fit,
                             mw,
                             rdata_path,
                             feat_name_fnc = get_human_int_feature_name,
                             topk = 30,
                             overwrite = SHOULD_OVERWRITE) {
  print(paste0("Consensus variable importance using output path: ", rdata_path))
  # Check if cached model exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Order ff_cor by 'abs_med_diff'
  ff_order = order(ff_cor$abs_med_diff, decreasing = TRUE)
  ff_cor = ff_cor[ff_order,]
  ff_cor = ff_cor[1:topk,]
  
  # Random forest mean importance
  rf.imp = importance(x = rf.fit, scale = FALSE)
  
  # Random forest mean decrease in accuracy
  rf.mda = rf.imp[, "MeanDecreaseAccuracy"]
  rf.mda.order = order(rf.mda, decreasing = TRUE)
  rf.mda = rf.mda[rf.mda.order]
  rf.mda = rf.mda[1:topk]
  
  # Random forest mean decrease in gini
  rf.mdg = rf.imp[, "MeanDecreaseGini"]
  rf.mdg.order = order(rf.mdg, decreasing = TRUE)
  rf.mdg = rf.mdg[rf.mdg.order]
  rf.mdg = rf.mdg[1:topk]
  
  # Mann-Whitney should already be sorted by FDR
  mw = mw[1:topk,]
  
  # Determine intersection among all sets
  common_vars = intersect(rownames(ff_cor), rownames(mw))
  common_vars = intersect(common_vars, names(rf.mda))
  common_vars = intersect(common_vars, names(rf.mdg))
  
  mean_ranks = rep(0.0, length(common_vars))
  names(mean_ranks) = common_vars
  
  # Compute mean rank
  for (varname in common_vars) {
    ff_cor_rank = which(rownames(ff_cor) == varname)
    mw_rank = which(rownames(mw) == varname)
    rf.mda.rank = which(names(rf.mda) == varname)
    rf.mdg.rank = which(names(rf.mdg) == varname)
    mean_ranks[varname] = mean(c(ff_cor_rank, mw_rank, rf.mda.rank, rf.mdg.rank))
  }
  mean_ranks = sort(mean_ranks)
  
  # Munge data.frames
  # HACK - all original summary meta columns are in uppercase
  # but all forest floor specific columns are lowercase, so
  # select that subset, prefix it, and use it as the base
  # output data.frame
  rank_names = names(mean_ranks)
  output = ff_cor[rank_names, ]
  ff_col_sel = colnames(ff_cor) == tolower(colnames(ff_cor))
  colnames(output)[ff_col_sel] = paste0("ff.", colnames(ff_cor)[ff_col_sel])
  
  output = data.frame(
    var.name = feat_name_fnc(output),
    mean.rank = mean_ranks,
    rf.mda = rf.mda[rank_names],
    rf.mdg = rf.mdg[rank_names],
    mw.pval = mw[rank_names, "PVAL"],
    mw.fdr = mw[rank_names, "FDR"],
    mw.bon = mw[rank_names, "BON"],
    output
  )
  
  # Save data
  save_rdata(output, rdata_path)
  csv_path = get_csv_path_from_rdata(rdata_path)
  print(paste0("...Saving to: ", csv_path))
  write.csv(x = output, file = csv_path)
  return(output)
}

# Example wild type analysis for consensus variable importance
do_cons_var_imp_wt_rcb <-
  function(keep.ls = c("ALL.ION"),
           pH.ls = c(5, 7),
           topk.ls = c(30)) {
    # Forest floor feature contributions
    ff_cor_ls = do_ff_cor_wt_rcb(keep.ls, pH.ls)
    
    # Random forest
    rf.fits = do_rf_wt_rcb(keep.ls, pH.ls)
    
    # Mann-whitney
    mw_paths_wt_rcb = list(
      "5" = list(
        "ALL.ION" = WT_PH5_MW_IMPUT_NA_ALL_ION_PATH_RCB,
        "LOOP.ION" = WT_PH5_MW_IMPUT_NA_LOOP_ION_PATH_RCB,
        "XORLOOP.ION" = WT_PH5_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB
      ),
      "7" = list(
        "ALL.ION" = WT_PH7_MW_IMPUT_NA_ALL_ION_PATH_RCB,
        "LOOP.ION" = WT_PH7_MW_IMPUT_NA_LOOP_ION_PATH_RCB,
        "XORLOOP.ION" = WT_PH7_MW_IMPUT_NA_XORLOOP_ION_PATH_RCB
      )
    )
    
    # Electrostatic data
    elec_by_pH = list("5" = WT_PH5_ELEC_RCB, "7" = WT_PH7_ELEC_RCB)
    
    # Open-close data
    oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
    
    # Output data
    results = list()
    
    # Iterate over pH
    for (pH in pH.ls) {
      pH_str = as.character(pH)
      # Skip missing pH levels
      if (!(pH_str %in% names(elec_by_pH))) {
        next
      }
      
      elec = elec_by_pH[[pH_str]]
      oc = oc_by_pH[[pH_str]]
      
      # Iterate over keep setting
      for (keep in keep.ls) {
        rf_res_id = paste0("WT_PH", pH, "_RF_", keep, "_RCB")
        rf.fit = rf.fits[[rf_res_id]]
        
        ff_cor_res_id = paste0("WT_PH", pH, "_FFCOR_", keep, "_RCB")
        ff_cor = ff_cor_ls[[ff_cor_res_id]]
        ff_cor = ff_cor$ff_cor
        
        mw_path = mw_paths_wt_rcb[[pH_str]][[keep]]
        mw = get_mann_whitney(elec$design,
                              oc,
                              elec$meta,
                              mw_path,
                              replace.na = TRUE,
                              keep = keep)
        
        # Iterate over topk setting
        for (topk in topk.ls) {
          var_imp_rdata_path = get_cons_var_imp_path(
            base_output_dir = get_relax_c_beta_dir(),
            pH = pH,
            keep = keep,
            sim_id = "wt",
            topk = topk
          )
          var_imp_res_id = paste0("WT_PH", pH, "_VARIMP", topk, "_", keep, "_RCB")
          
          results[[var_imp_res_id]] = get_cons_var_imp(
            ff_cor = ff_cor,
            rf.fit = rf.fit,
            mw = mw,
            rdata_path = var_imp_rdata_path,
            topk = topk
          )
          
          
        } #end iteration over topk
      } # end iteration over keep
    } # end iteration over pH
    
    return(results)
  }

###################################################################
# Bootstrap stability analysis
###################################################################

library(boot)

# Number of bootstrap trials for each bootstrap estimate
NUM_BOOTSTRAP_TRIALS = 10000

# Number of parallel estimates
NUM_BOOTSTRAP_CPUS = 18

# Confidence interval (must be greater than 0 and less than 1)
BOOTSTRAP_CONF = 0.95

# The default bootstrap statistic to use
# Options - "MED_DIFF" | "MEA_DIFF" | "MEDQ_DOPEN" | "MEDQ_DCLOSE"
BOOTSTRAP_STAT_TYPE = "MED_DIFF"

# @return Name of bootstrap rdata file
get_boot_rdata_name <- function(base.name,
                                replace.na = FALSE,
                                replace.val = 0.0,
                                num.trials = NUM_BOOTSTRAP_TRIALS,
                                max_energy_rank = MAX_ENERGY_RANK,
                                stat_type = BOOTSTRAP_STAT_TYPE,
                                conf = BOOTSTRAP_CONF) {
  return(
    paste0(
      base.name,
      ".nar",
      sum(replace.na),
      ".nav",
      replace.val,
      ".nt",
      num.trials,
      ".bt.r",
      max_energy_rank,
      paste0(".", stat_type[1]),
      ".ci",
      conf,
      ".rdata"
    )
  )
}

# @return Path to bootstrap confidence interval rdata
#   Path is meant to store result of boot.ci() call
get_boot_path <- function(base_output_dir = get_relax_c_beta_dir(),
                          pH = 5,
                          sim_id = "wt",
                          replace.na = FALSE,
                          replace.val = 0.0,
                          num.trials = NUM_BOOTSTRAP_TRIALS,
                          max_energy_rank = MAX_ENERGY_RANK,
                          stat_type = BOOTSTRAP_STAT_TYPE,
                          conf = BOOTSTRAP_CONF) {
  pH_str = paste0("pH", pH)
  base.name = paste0(sim_id, ".", pH_str)
  file.path(
    base_output_dir,
    "capt",
    pH_str,
    sim_id,
    get_boot_rdata_name(
      base.name = base.name,
      replace.na = replace.na,
      replace.val = replace.val,
      num.trials = num.trials,
      max_energy_rank = max_energy_rank,
      stat_type = stat_type[1],
      conf = conf
    )
  )
}

# @return path to false-coverage rate data.frame from boot_path
get_boot_fcr_path <- function(boot_path, flt_type) {
  return(insert_attrib_in_rdata(rdata_path = boot_path,
                                attrib = paste0(flt_type, ".fcr")))
}

# @param boot_rdata_path - Path to confidence interval data frame.
#   Generated via call to get_boot_path()
# @param stat_type - Statistic type (eg. MED_DIFF)
# @return Directory for storing replicate data for each CI in data frame.
#   Note: replicate data is generated via boot() calls
get_boot_rep_dir <- function(boot_rdata_path, stat_type) {
  return(file.path(dirname(boot_rdata_path), paste0("boot.reps.", stat_type)))
}

# @param boot_rdata_path - Path to confidence interval data frame.
#   Generated via call to get_boot_path()
# @param element_id - String with name of element
# @return Name of replicate data file
#   Note: replicate data is generated via boot() calls
get_boot_rep_fname <- function(boot_rdata_path, element_id) {
  fid = get_fid_from_rdata(boot_rdata_path)
  return(paste0(fid, ".", element_id, ".rdata"))
}

# @param boot_rdata_path - Path to confidence interval data frame.
#   Generated via call to get_boot_path()
# @param stat_type - Statistic type (eg. MED_DIFF)
# @param element_id - String with name of element
# @return Path to replicate data derived from CI boot data path
#   Note: replicate data is generated via boot() calls
get_boot_rep_path <-
  function(boot_rdata_path, stat_type, element_id) {
    dir = get_boot_rep_dir(boot_rdata_path, stat_type)
    fname = get_boot_rep_fname(boot_rdata_path, element_id)
    return(file.path(dir, fname))
  }

##############################
# Bootstrap - ignore NAs

# WT RCB

WT_PH5_BT_IGNOR_NA_PATH_RCB = get_boot_path(base_output_dir = get_relax_c_beta_dir(),
                                            pH = 5,
                                            sim_id = "wt")

WT_PH7_BT_IGNOR_NA_PATH_RCB = get_boot_path(base_output_dir = get_relax_c_beta_dir(),
                                            pH = 7,
                                            sim_id = "wt")

##############################
# Bootstrap - impute NAs

# WT RCB

WT_PH5_BT_IMPUT_NA_PATH_RCB = get_boot_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 5,
  sim_id = "wt",
  replace.na = TRUE
)

WT_PH7_BT_IMPUT_NA_PATH_RCB = get_boot_path(
  base_output_dir = get_relax_c_beta_dir(),
  pH = 7,
  sim_id = "wt",
  replace.na = TRUE
)

##############################
# Bootstrap - favor utility

# Determines labels for favoring elements
# @param bt.ci - Data.frame of confidence intervals and meta data
# @param destab - Logical vector for destabilizing elements
# @param stab - Logical vector for stabilizing elements
# @param fav.open - Logical vector indicating elements favoring open state
# @param fav.close - Logical vector indicating elements favoring close state
# @return modified bt.ci with labels for which state is favored and if
#   that interaction is stabilizing or destabilizing
boot_fav_util <-
  function(bt.ci, destab, stab, fav.open, fav.close) {
    if (length(fav.open) > 0) {
      print(paste0(sum(fav.open), " elements favor open conformation."))
      bt.ci$fav[fav.open] = "OPEN"
      # Check if any are destabilizing elements for the closed conformation
      destab.close = destab & fav.open
      if (length(destab.close) > 0) {
        bt.ci$destab[destab.close] = "CLOSE"
      }
      # Check if any are stabilizing elements for the open conformation
      stab.open = stab & fav.open
      if (length(stab.open) > 0) {
        bt.ci$stab[stab.open] = "OPEN"
      }
    }
    
    # Check if any elements favor the close conformation
    if (length(fav.close) > 0) {
      print(paste0(sum(fav.close), " elements favor closed conformation."))
      bt.ci$fav[fav.close] = "CLOSE"
      # Check if any are destabilizing elements for the open conformation
      destab.open = destab & fav.close
      if (length(destab.open) > 0) {
        bt.ci$destab[destab.open] = "OPEN"
      }
      # Check if any are stabilizing elements for the closed conformation
      stab.close = stab & fav.close
      if (length(stab.close) > 0) {
        bt.ci$stab[stab.close] = "CLOSE"
      }
    }
    
    return(bt.ci)
  }

##############################
# Bootstrap - Quant utility

# Determines the normalized quantile of param val within param density (pdf)
# @param dens_data - Input data set to generate a pdf density from
# @param val - The value to determine P(X < val) where P is the parameter
#   density
# @param scale - All values are scaled by this value for numerical stability
# @param init_subd - The initial number of subdivisions used for numerical
#   integration of the pdf. This value must be a positive integer!
# @param subd_mult - If integration failed, subdivision is increased to
#   subdivision *= subd_mult. This value must be a positive integer!
# @param max_iters - The maximum number of iterations to attempt a succesful
#   numerical integration. If exceeded, the last value of the integrate() call
#   is returned and a warning is output to console.
# @return The normalized quantile in [0,1] corresponding to P(X < val) where P
#   is pdf defined by dens parameter
boot_quant_util <- function(dens_data,
                            val,
                            scale = 100.0,
                            init_subd = 100L,
                            subd_mult = 2L,
                            max_iters = 5) {
  # Special case: early out if data is degenerate (all 0.0)
  if (all(dens_data == 0.0)) {
    if (val < 0.0) {
      return(0.0)
    } else if (val > 0.0) {
      return(1.0)
    } else {
      # val == 0.0
      return(0.5)
    }
  }
  # Scale inputs
  val = scale * val
  dens_data = scale * dens_data
  # Create pdf density
  dens = density(x = dens_data)
  # Determine normalized quantile
  quant = 0.0
  # Only integrate if val is within mapped (lo, hi) density range
  if (dens$x[1] < val) {
    quant = 1.0
    if (val < dens$x[length(dens$x)]) {
      # Obtain probability density
      pdf = approxfun(
        x = dens$x,
        y = dens$y,
        method = "linear",
        yleft = 0.0,
        yright = 0.0
      )
      # Integrate to obtain P(X < val)
      subd = init_subd
      quant.res = integrate(
        f = pdf,
        lower = dens$x[1],
        # can also try -Inf
        upper = val,
        subdivisions = subd,
        stop.on.error = FALSE
      )
      # Integration failed, increase number of subdivisions
      num_iters = 1
      while (quant.res$message != "OK") {
        # Exit loop if exceed iteration cap
        if (num_iters >= max_iters) {
          break
        }
        num_iters = num_iters + 1
        subd = subd * subd_mult
        quant.res = integrate(
          f = pdf,
          lower = dens$x[1],
          upper = val,
          subdivisions = subd,
          stop.on.error = FALSE
        )
      }
      # Log if final integration status is not okay!
      # @TODO - figure out best way to handle different messages
      # - roundoff error - probably okay
      # - divergent - is there a way to trap this and just manually resolve?
      # - not enough subdivisions - probably okay
      # - malformed integral - ...
      #if (quant.res$message != "OK") {
      #  print(paste0("WARNING, integration message: ", quant.res$message))
      #}
      quant = quant.res$value
      
      # Check if quant is finite
      if (!is.finite(quant)) {
        # Quantile is not finite - estimate based on position within
        # input density rdata.
        quant = rank(c(dens_data, val))
        quant = quant / length(quant)
        quant = quant[length(quant)]
        
        # Note sure how we could end up here, but as last resort
        # just assign 0.5!
        if (!is.finite(quant)) {
          print("WARNING, non-finite quantile!")
          quant = 0.5
        } # end second check to see if quant is finite
      } # end check if quant is finite
      
      # Make sure quantile is in [0, 1]
      quant = max(min(quant, 1.0), 0.0)
      
    } # end check if val param within density range
  }
  
  return(quant)
}

##############################
# Bootstrap - MED_DIFF and MEA_DIFF

# Corresponds to stat_type = "MED_DIFF"
# @return median difference: med.open - med.close
boot_median_diff <- function(data, i, n.open) {
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  med.open = median(data[i.open])
  med.close = median(data[i.close])
  med.diff = med.open - med.close
  return(med.diff)
}

# Corresponds to stat_type = "MEA_DIFF"
# @return mean difference: mea.open - mea.close
boot_mean_diff <- function(data, i, n.open) {
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  mea.open = mean(data[i.open])
  mea.close = mean(data[i.close])
  mea.diff = mea.open - mea.close
  return(mea.diff)
}

# Determines if any interactions favor one conformation
# over an other given the computed confidence interval.
# @param bt.ci - Data.frame of confidence intervals and meta data
# @param destab - Logical vector for destabilizing interactions
# @param stab - Logical vector for stabilizing interactions
# @param lo_col - Name of column containing lower bound of interval
# @parma hi_col - Name of column containing upper boud of interval
# @return bt.ci modified with favoring state and (de)stabilizing labels
boot_median_or_mean_diff_post <-
  function(bt.ci,
           destab,
           stab,
           lo_col = "lo",
           hi_col = "hi") {
    # See if any of the confidence intervals are interesting:
    
    # Check if any interactions favor the open conformation
    # Favor -> more stabilizing or less destabilizing relative to opposite conformation
    fav.open = (bt.ci[, hi_col] < 0)
    
    # Check if any interactions favor the close conformation
    fav.close = (bt.ci[, lo_col] > 0)
    
    # Mark favored conformations
    bt.ci = boot_fav_util(
      bt.ci = bt.ci,
      destab = destab,
      stab = stab,
      fav.open = fav.open,
      fav.close = fav.close
    )
    
    # Insert original (non-processed) observed data
    bt.ci[, "obs.med_diff"] = bt.ci$MED.OPEN - bt.ci$MED.CLOSE
    bt.ci[, "obs.mea_diff"] = bt.ci$MEA.OPEN - bt.ci$MEA.CLOSE
    
    return(bt.ci)
  }

##############################
# Bootstrap - MEDQ_DOPEN

# Corresponds to stat_type = "MEDQ_DOPEN"
# @return normalized quantile of med.close in open state distribution
boot_median_quant_distrib_open <- function(data, i, n.open) {
  # http://stats.stackexchange.com/questions/6940/calculating-p-value-from-an-arbitrary-distribution
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  med.close = median(data[i.close])
  return(boot_quant_util(dens_data = data[i.open], val = med.close))
}

# Determines if any interactions favor one conformation
# over an other given the computed confidence interval.
# @param bt.ci - Data.frame of confidence intervals and meta data
# @param destab - Logical vector for destabilizing interactions
# @param stab - Logical vector for stabilizing interactions
# @param lo_col - Name of column containing lower bound of interval
# @parma hi_col - Name of column containing upper boud of interval
# @return bt.ci modified with favoring state and (de)stabilizing labels
boot_median_quant_distrib_open_post <-
  function(bt.ci,
           destab,
           stab,
           lo_col = "lo",
           hi_col = "hi") {
    # In this case, we want to know if the confidence interval crosses
    # the 0.5 quantile. If it doesn't, then we have some support
    # that this interaction is favored in one conformation relative
    # to the other.
    
    # Check for favoritism:
    # Close median is generally above 0.5 quantile in open distribution
    fav.open = bt.ci[, lo_col] > 0.5
    # # Close median is generally below 0.5 quantile in open distribution
    fav.close = bt.ci[, hi_col] < 0.5
    
    return(
      boot_fav_util(
        bt.ci = bt.ci,
        destab = destab,
        stab = stab,
        fav.open = fav.open,
        fav.close = fav.close
      )
    )
  }

##############################
# Bootstrap - MEDQ_DCLOSE

# Corresponds to stat_type = "MEDQ_DCLOSE"
# @return normalized quantile of med.open in closed state distribution
boot_median_quant_distrib_close <- function(data, i, n.open) {
  # http://stats.stackexchange.com/questions/6940/calculating-p-value-from-an-arbitrary-distribution
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  med.open = median(data[i.open])
  return(boot_quant_util(dens_data = data[i.close], val = med.open))
}

# Determines if any interactions favor one conformation
# over an other given the computed confidence interval.
# @param bt.ci - Data.frame of confidence intervals and meta data
# @param destab - Logical vector for destabilizing interactions
# @param stab - Logical vector for stabilizing interactions
# @param lo_col - Name of column containing lower bound of interval
# @parma hi_col - Name of column containing upper boud of interval
boot_median_quant_distrib_close_post <-
  function(bt.ci,
           destab,
           stab,
           lo_col = "lo",
           hi_col = "hi") {
    # In this case, we want to know if the confidence interval crosses
    # the 0.5 quantile. If it doesn't, then we have some support
    # that this interaction is favored in one conformation relative
    # to the other.
    
    # Check for favoritism:
    # Open median is generally below 0.5 quantile in closed distribution
    fav.open = bt.ci[, hi_col] < 0.5
    # # Open median is generally above 0.5 quantile in closed distribution
    fav.close = bt.ci[, lo_col] > 0.5
    
    return(
      boot_fav_util(
        bt.ci = bt.ci,
        destab = destab,
        stab = stab,
        fav.open = fav.open,
        fav.close = fav.close
      )
    )
  }

# Filter which performs no filtering (keeps everything)
boot_no_filt <- function(meta) {
  filt = rep(TRUE, nrow(meta))
  return(filt)
}

# Interaction filter which only keeps ion-ion interactions
#   where at least one of the ions is within a loop
boot_int_ion_ion_with_intra_filt <- function(meta) {
  # First keep only ion-ion interactions
  filt = meta$ION.ION
  # Now, make sure at least one of the ions is within a loop region
  filt = filt &
    ((meta$RESA.LOOP != "NONE") | (meta$RESB.LOOP != "NONE"))
  return(filt)
}

# Interaction filter which only keeps ion-ion interactions,
#   where at least one of the ions is within a loop and
#   both ions are not within the same loop
boot_int_ion_ion_no_intra_filt <- function(meta) {
  filt = boot_int_ion_ion_with_intra_filt(meta)
  # Finally, remove any ion-ion interactions where both are in same loop
  filt = filt & (meta$RESA.LOOP != meta$RESB.LOOP)
  return(filt)
}

# Default interaction filter
boot_int_filt <- boot_no_filt

# SRS filter which only considers loop atoms
boot_srs_loop_filt <- function(meta) {
  filt = (meta$RES.LOOP != "NONE")
  return(filt)
}

# Default SRS (single residue score) filter
boot_srs_filt <- boot_no_filt

# Default interaction stability
boot_int_stab <- function(bt.ci,
                          lo_col = "lo",
                          hi_col = "hi") {
  POS_IONS = c('ARG', 'HIS', 'LYS')
  NEG_IONS = c('ASP', 'GLU')
  
  # Flag destabilizing interactions
  destab = ((bt.ci$RESA.NAME %in% POS_IONS) &
              (bt.ci$RESB.NAME %in% POS_IONS)) |
    ((bt.ci$RESA.NAME %in% NEG_IONS) &
       (bt.ci$RESB.NAME %in% NEG_IONS))
  # Flag stabilizing interactions
  stab = !destab
  return(list(destab = destab, stab = stab))
}

# SRS stability based on median value
boot_srs_median_stab <- function(bt.ci,
                                 lo_col = "lo",
                                 hi_col = "hi") {
  # Use medians to determine if destabilizing
  fav.open = (bt.ci[, hi_col] < 0)
  fav.close = (bt.ci[, lo_col] > 0)
  
  stab = fav.open & bt.ci$MED.OPEN < 0
  stab[fav.close] = bt.ci$MED.CLOSE[fav.close] < 0
  
  destab = fav.open & bt.ci$MED.OPEN > 0
  destab[fav.close] = bt.ci$MED.CLOSE[fav.close] > 0
  
  return(list(destab = destab, stab = stab))
}

##############################
# Bootstrap - Compute CI

# @param design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples.
# @param oc - logical vector with TRUE = open, FALSE = close.
# @param meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSe - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#   See description in load_data_utils::elec_csv_2_rdata() for full column listings
# @param rdata_path - Output path for bootstrap confidence interval rdata
# @param replace.na - If TRUE, NAs are replaced with replace.val, else they are removed
# @param replace.val - Only used if replace.na is TRUE - all NAs are replaced with this value
# @param overwrite - If true, will always recompute analysis and overwrite anything already
#   stored to disk. Else, will check for previously computed rdata and return the cached result.
# @param num.trials - The number of bootstrap trials to run for each CI estimate
# @param num.cpus - The number of cpus to use for parallelization
# @param stat_type - Current options for bootstrap statistic are
#   "MED_DIFF" - statistic computes Med.open - Med.close (median difference)
#   "MEA_DIFF" - statistic computes Mea.open - Mea.close (mean difference)
#   "MEDQ_DOPEN" - statistic computes P(X.open < Med.close) where P is open probability density
#   "MEDQ_DCLOSE" - statistic computes P(X.close < Med.open) where P is close probability density
#   "MED_DIFF_SRS" - median difference but with single residue instead of interaction data
# @param conf - The confidence interval range (scalar) in [0,1]
# @return data.frame with the following columns
#   $lo - the lower bound for the 95% CI
#   $hi - the upper bound for the 95% CI
#   $t0 - result of applying statistic to observed (possibly imputed) data
#   $obs - observed value of statistic according to input meta data (filtered NAs)
#   $fav - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more stabilizing or less destabilizing for that conformation (according to 95% CI)
#   $stab - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more stabilizing for that conformation (according to 95% CI)
#   $destab - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more destabilizing for that conformation (according to 95% CI)
#   ... also all columns from meta parameter are included in output
get_boot_ci <- function(design,
                        oc,
                        meta,
                        rdata_path,
                        replace.na = FALSE,
                        replace.val = 0.0,
                        overwrite = SHOULD_OVERWRITE,
                        num.trials = NUM_BOOTSTRAP_TRIALS,
                        num.cpus = NUM_BOOTSTRAP_CPUS,
                        stat_type = BOOTSTRAP_STAT_TYPE,
                        conf = BOOTSTRAP_CONF) {
  # Make sure confidence interval is valid
  stopifnot(conf >= 0)
  stopifnot(conf < 1)
  
  # Statistic options
  stat_opts = list(
    "MED_DIFF" = list(
      fn = boot_median_diff,
      flt = boot_int_filt,
      stb = boot_int_stab,
      post = boot_median_or_mean_diff_post
    ),
    "MEA_DIFF" = list(
      fn = boot_mean_diff,
      flt = boot_int_filt,
      stb = boot_int_stab,
      post = boot_median_or_mean_diff_post
    ),
    "MEDQ_DOPEN" = list(
      fn = boot_median_quant_distrib_open,
      boot_int_filt = boot_int_filt,
      stb = boot_int_stab,
      post = boot_median_quant_distrib_open_post
    ),
    "MEDQ_DCLOSE" = list(
      fn = boot_median_quant_distrib_close,
      flt = boot_int_filt,
      stb = boot_int_stab,
      post = boot_median_quant_distrib_close_post
    ),
    "MED_DIFF_SRS" = list(
      fn = boot_median_diff,
      flt = boot_srs_filt,
      stb = boot_srs_median_stab,
      post = boot_median_or_mean_diff_post
    )
  )
  
  # Check statistic
  stopifnot(stat_type %in% names(stat_opts))
  
  # Obtain statistic function
  stat_fn = stat_opts[[stat_type]]$fn
  # Obtain sample filter
  stat_flt = stat_opts[[stat_type]]$flt
  # Obtain stability callback
  stat_stb = stat_opts[[stat_type]]$stb
  # Obtain post function called after all replicates have been computed
  stat_post = stat_opts[[stat_type]]$post
  
  print(paste0(
    "Bootstrap confidence intervals using output path: ",
    rdata_path
  ))
  # Check if cached model exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping bootstrap as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Create directory for storing bootstrap replicates
  rep_dir = get_boot_rep_dir(rdata_path, stat_type)
  print(paste0("...Creating bootstrap replicate directory at: ", rep_dir))
  dir.create(rep_dir, showWarnings = FALSE)
  
  # Filter covariates if only certain elements are desired:
  filt = stat_flt(meta)
  
  # Apply filter
  meta.filt = meta[filt, ]
  elem_names = rownames(meta.filt)
  design = design[, elem_names]
  
  # Total number of elements
  n.elem = ncol(design)
  print(paste0("...After filtering, ", n.elem, " elements remain."))
  
  # Check if we need to replace missing values
  if (replace.na) {
    print(paste0("...Replacing NAs with ", replace.val))
    design[is.na(design)] = replace.val
  }
  
  # Pre-separate open and closed samples. Hopefully enough memory for this.
  # The idea is to avoid repeatedly filtering within the loop.
  design.open = design[oc, ]
  design.close = design[!oc, ]
  
  # The output bootstrap confidence intervals data frame
  bt.ci = data.frame(
    lo = rep(NA, n.elem),
    hi = rep(NA, n.elem),
    t0 = rep(NA, n.elem),
    fav = rep("NONE", n.elem),
    stab = rep("NONE", n.elem),
    destab = rep("NONE", n.elem),
    meta.filt,
    stringsAsFactors = FALSE
  )
  
  # Avoid testing samples with less than this number of samples
  # in both the open and closed classes.
  MIN_NUM_SAMPLES_PER_CLASS = 30
  
  # @HACK: strip CI from path name as replicate data should not
  # be dependent on interval
  boot_rep_base_rdata_path = gsub("\\.ci.*\\.rdata$", ".rdata", rdata_path)
  
  print("...Performing bootstrap")
  for (i.int in 1:n.elem) {
    # Determine path to replicate data
    element_id = colnames(design)[i.int]
    rep_fname = get_boot_rep_fname(boot_rep_base_rdata_path, element_id)
    rep_rdata_path = file.path(rep_dir, rep_fname)
    
    # Check if replicates already exist for this element
    rep_data = c()
    if (!overwrite && file.exists(rep_rdata_path)) {
      print(
        paste0(
          "...Skipping bootstrap replicate as rdata already exists: ",
          rep_rdata_path
        )
      )
      rep_data = load_rdata(rep_rdata_path)
    }
    else {
      # No cached replicate data exists, must compute
      x.open = design.open[, i.int]
      y.close = design.close[, i.int]
      
      # Filter missing values
      x.open = x.open[!is.na(x.open)]
      y.close = y.close[!is.na(y.close)]
      
      n.open = length(x.open)
      n.close = length(y.close)
      
      # Check if we have enough values to proceed
      if (n.open < MIN_NUM_SAMPLES_PER_CLASS) {
        next
      }
      if (n.close < MIN_NUM_SAMPLES_PER_CLASS) {
        next
      }
      # Create bootstrap replicate data
      rep_data = boot(
        data = c(x.open, y.close),
        statistic = stat_fn,
        strata = rep(x = c(1, 2),
                     times = c(n.open,
                               n.close)),
        R = num.trials,
        parallel = "multicore",
        ncpus = num.cpus,
        n.open = n.open
      )
      
      # Save replicate data
      save_rdata(rep_data, rep_rdata_path)
    }
    # Compute confidence intervals
    # @TODO - look into BCA (bootstrap control accelerated)
    
    ci = boot.ci(rep_data, conf = conf, type = "basic")
    if (!is.null(ci)) {
      ci = ci$basic
      bt.ci$lo[i.int] = ci[length(ci) - 1]
      bt.ci$hi[i.int] = ci[length(ci)]
    } else {
      # http://stackoverflow.com/questions/16652852/avoid-error-in-r-boot-ci-function-when-all-values-in-sampled-set-are-equal
      # Replace with observed value if CI cannot be computed due
      # to lack of uniqueness
      bt.ci$lo[i.int] = rep_data$t0
      bt.ci$hi[i.int] = rep_data$t0
    }
    # Store observed value of statistic applied to input data
    # See ?boot
    bt.ci$t0[i.int] = rep_data$t0
  }
  
  # Remove any rows which were skipped due to not enough samples
  bt.ci = bt.ci[!is.na(bt.ci$lo) & !is.na(bt.ci$hi), ]
  
  # Determine stabilizing and destabilizing elements
  stb_ls = stat_stb(bt.ci)
  destab = stb_ls$destab
  stab = stb_ls$stab
  
  # Post replicate call - meant to check if any of the confidence intervals
  # are interesting (favor one conformation over the other)
  bt.ci = stat_post(bt.ci = bt.ci,
                    destab = destab,
                    stab = stab)
  
  # Save results
  save_rdata(bt.ci, rdata_path)
  write.csv(x = bt.ci,
            quote = FALSE,
            file = get_csv_path_from_rdata(rdata_path))
  return (bt.ci)
}

##############################
# Bootstrap - Compute FCR

# Determine false coverage rate (FCR)
# Assumes 0 is the critical point of interest
# https://en.wikipedia.org/wiki/False_coverage_rate
# https://en.wikipedia.org/wiki/False_discovery_rate
# @param bt.ci - data.frame with the following columns
#   $lo - the lower bound for the 95% CI
#   $hi - the upper bound for the 95% CI
#   $t0 - result of applying statistic to observed (possibly imputed) data
#   $obs - observed value of statistic according to input meta data (filtered NAs)
#   $fav - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more stabilizing or less destabilizing for that conformation (according to 95% CI)
#   $stab - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more stabilizing for that conformation (according to 95% CI)
#   $destab - can be 'NONE', 'OPEN', or 'CLOSE'. If 'OPEN' or 'CLOSE', then the interaction
#     is more destabilizing for that conformation (according to 95% CI)
#   ... also all columns from meta parameter are included in output
# @param bt.ci_rdata_path - The path used for loading bt.ci from rdata. Needed for
#   obtaining the bootstrap replicate data
# @param alpha - Significance level for FCR correction
# @param flt - Filter samples prior to p-value and corresponding FDR correction
# @param flt_desc - Character vector with description of filter type, used for determining
#   default output path
# @param rdata_path - Output path for FCR data
# @param stat_type - Current options for bootstrap statistic are
#   "MED_DIFF" - statistic computes Med.open - Med.close (median difference)
#   "MEA_DIFF" - statistic computes Mea.open - Mea.close (mean difference)
#   "MED_DIFF_SRS" - median difference but with single residue instead of interaction data
# @param stb - Stability filter - use boot_int_stab | boot_srs_median_stab
#   determines if element is stabilizing or destabilizing
# @param post - Post-process method - used for assigning open or close favorability
# @return list with element 'bt' as data.frame with appended columns
#   $pval - p-value defined as min(% >= 0, % <= 0)
#   $fdr - false discovery rate adjusted p-value
#   $bon - Bonferroni correction
#   $lo.fcr - lower bound of false coverage rate interval
#   $hi.fcr - upper bound of false coverage rate interval
#   $fav.fcr - favorability according to fcr
# list will also append .ci suffix to non-fcr confidence interval columns lo, hi, fav
# list will also contain 'path' element for the parameter rdata_path used
adjust_boot_ci <- function(bt.ci,
                           bt.ci_rdata_path,
                           alpha = 0.05,
                           stat_type = BOOTSTRAP_STAT_TYPE,
                           stb = boot_int_stab,
                           post = boot_median_or_mean_diff_post,
                           flt = boot_int_ion_ion_with_intra_filt,
                           flt_desc = "ion.only",
                           rdata_path = get_boot_fcr_path(bt.ci_rdata_path, flt_desc),
                           overwrite = SHOULD_OVERWRITE) {
  # Only statistics with 0-critical point are allowed
  stopifnot(stat_type %in% c("MED_DIFF", "MEA_DIFF", "MED_DIFF_SRS"))
  
  print(paste0("Adjusting bootstrap FCR using output path: ", rdata_path))
  # Check if cached data exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping FCR adjustment as rdata already exists: ",
      rdata_path
    ))
    return(list(bt = load_rdata(rdata_path),
                path = rdata_path))
  }
  
  # Filter covariates if only certain elements are desired:
  filt = flt(bt.ci)
  bt.ci = bt.ci[filt, ]
  
  # Total number of elements
  n.elem = nrow(bt.ci)
  print(paste0("...After filtering, ", n.elem, " elements remain."))
  
  # Rename (1 - alpha) confidence interval columns
  lo_col_ix = which(colnames(bt.ci) == 'lo')
  hi_col_ix = which(colnames(bt.ci) == 'hi')
  fav_col_ix = which(colnames(bt.ci) == 'fav')
  colnames(bt.ci)[lo_col_ix] = 'lo.ci'
  colnames(bt.ci)[hi_col_ix] = 'hi.ci'
  colnames(bt.ci)[fav_col_ix] = 'fav.ci'
  # Drop redundant stab and destab columns
  bt.ci = bt.ci[,!(colnames(bt.ci) %in% c("stab", "destab"))]
  
  # @HACK: strip CI from path name as replicate data should not
  # be dependent on interval
  boot_rep_base_rdata_path = gsub("\\.ci.*\\.rdata$", ".rdata", bt.ci_rdata_path)
  # Determine directory where replicates are stored
  rep_dir = get_boot_rep_dir(bt.ci_rdata_path, stat_type)
  
  # The output FCR data frame
  na_vec = rep(NA, n.elem)
  none_vec = rep("NONE", n.elem)
  bt.fcr = data.frame(
    lo.fcr = na_vec,
    hi.fcr = na_vec,
    fdr = na_vec,
    bon = na_vec,
    pval = na_vec,
    fav = none_vec,
    stab = none_vec,
    destab = none_vec,
    conf.fcr = na_vec,
    conf.ci = rep(1.0 - alpha, n.elem),
    bt.ci,
    stringsAsFactors = FALSE
  )
  
  # Utility to load the i-th set of replicates
  load_rep <- function(i.elem) {
    # Determine path to replicate data
    element_id = rownames(bt.ci)[i.elem]
    rep_fname = get_boot_rep_fname(boot_rep_base_rdata_path, element_id)
    rep_rdata_path = file.path(rep_dir, rep_fname)
    # Check if replicates exist for this element
    rep_data = NULL
    if (!file.exists(rep_rdata_path)) {
      print(paste0("...Skipping as no replicates exist for: ",
                   rep_rdata_path))
    } else {
      # Load replicate data
      rep_data = load_rdata(rep_rdata_path)
    }
    return(rep_data)
  }
  
  print("...Computing p-value")
  for (i.elem in 1:n.elem) {
    rep_data = load_rep(i.elem)
    if (is.null(rep_data)) {
      next
    }
    # Munge observed stat and boot stats into single set
    reps = c(rep_data$t0, rep_data$t[, 1])
    
    # Calculate p-value as min(% <= 0, % >= 0)
    n = length(reps)
    bt.fcr$pval[i.elem] = min(sum(reps <= 0.0) / n, sum(reps >= 0.0) / n)
  }
  
  print("...Adjusting p-value")
  bt.fcr$fdr = p.adjust(p = bt.fcr$pval, method = "fdr")
  bt.fcr$bon = p.adjust(p = bt.fcr$pval, method = "bonferroni")
  
  # Determine number of "significant" elements according to alpha
  k = sum(bt.fcr$fdr < alpha)
  # Adjust confidence region
  conf = 1.0 - (k * alpha / n.elem)
  bt.fcr$conf.fcr = rep(conf, n.elem)
  print(paste0("...Computing ", conf, " FCR"))
  for (i.elem in 1:n.elem) {
    rep_data = load_rep(i.elem)
    if (is.null(rep_data)) {
      next
    }
    ci = boot.ci(rep_data, conf = conf, type = "basic")
    if (!is.null(ci)) {
      ci = ci$basic
      bt.fcr$lo.fcr[i.elem] = ci[length(ci) - 1]
      bt.fcr$hi.fcr[i.elem] = ci[length(ci)]
    } else {
      # http://stackoverflow.com/questions/16652852/avoid-error-in-r-boot-ci-function-when-all-values-in-sampled-set-are-equal
      # Replace with observed value if CI cannot be computed due
      # to lack of uniqueness
      bt.fcr$lo.fcr[i.elem] = rep_data$t0
      bt.fcr$hi.fcr[i.elem] = rep_data$t0
    }
  }
  
  # Determine stabilizing and destabilizing elements
  stb_ls = stb(bt.fcr, lo_col = "lo.fcr", hi_col = "hi.fcr")
  destab = stb_ls$destab
  stab = stb_ls$stab
  
  # Post replicate call - meant to check if any of the confidence intervals
  # are interesting (favor one conformation over the other)
  bt.fcr = post(
    bt.ci = bt.fcr,
    destab = destab,
    stab = stab,
    lo_col = "lo.fcr",
    hi_col = "hi.fcr"
  )
  
  # Rename favor .fcr suffix
  fav_col_ix = which(colnames(bt.fcr) == 'fav')
  colnames(bt.fcr)[fav_col_ix] = 'fav.fcr'
  
  # Save results
  save_rdata(bt.fcr, rdata_path)
  write.csv(x = bt.fcr,
            quote = FALSE,
            file = get_csv_path_from_rdata(rdata_path))
  return(list(bt = bt.fcr, path = rdata_path))
}

##############################
# Bootstrap - Plots

# http://stackoverflow.com/questions/4559229/drawing-pyramid-plot-using-r-and-ggplot2
# https://learnr.wordpress.com/2009/09/24/ggplot2-back-to-back-bar-charts/
# http://docs.ggplot2.org/0.9.3.1/geom_errorbar.html
# Generates bar plots of e.g. the median net stability
# @param bt.ci - Data.frame with boot strap confidence intervals
# @param main - Title of plot
# @param xlab - X-axis label
# @param ylab - Y-axis label
# @param topk - Only top-k features are plotted
# @param fill_lo - The low color value for the bar plot fill
# @param fill_hi - The high color vlaue for the bar plot fill
# See https://coolors.co/
# http://www.colorhexa.com/
plot_boot_ci_bar <- function(bt.ci,
                             feat_name_fnc = get_human_int_feature_name,
                             main = "Bootstrap Top Median Net Stability",
                             xlab = "OmpG Residue Interactions",
                             ylab = "Median Net Stability",
                             topk = min(35, nrow(bt.ci)),
                             # Yankees Blue
                             fill_lo = "#132B43",
                             # Blue Jeans
                             fill_hi = "#56B1F7",
                             # Mandarin
                             col_err = "#E57A44") {
  stopifnot(topk >= 1)
  library(ggplot2)
  # Order by decreasing median feature contribution
  bt.ci$abs_t0 = abs(bt.ci$t0)
  
  # Filter topK
  ix_order = order(bt.ci$abs_t0, decreasing = TRUE)
  bt.ci = bt.ci[ix_order, ]
  bt.ci = bt.ci[1:topk, ]
  feat.names = feat_name_fnc(bt.ci)
  bt.ci$feature = feat.names # as.factor(feat.names)
  
  # Determine error limits
  ix_neg_t0 = bt.ci$t0 < 0.0
  bt.ci$abs_lo = abs(bt.ci$lo)
  bt.ci$abs_hi = abs(bt.ci$hi)
  bt.ci$abs_lo[ix_neg_t0] = bt.ci$abs_hi[ix_neg_t0]
  bt.ci$abs_hi[ix_neg_t0] = abs(bt.ci$lo[ix_neg_t0])
  
  limits <-
    aes(ymax = abs_hi, ymin = abs_lo)
  
  dodge <- position_dodge(width = 0.9)
  
  p = ggplot(bt.ci,
             aes(
               fill = abs_lo,
               y = abs_lo,
               x = reorder(feature, abs_lo, function(x) {
                 -x
               })
             )) +
    geom_bar(stat = "identity",
             position = dodge) +
    scale_fill_gradient(low = fill_lo, high = fill_hi) +
    labs(x = xlab,
         y = ylab,
         title = main) +
    theme(axis.text.x = element_text(angle = -90),
          legend.position = "none") +
    geom_errorbar(limits,
                  position = dodge,
                  width = 0.25,
                  color = col_err)
  return(p)
}

##############################
# Bootstrap - Examples

# Example wildtype analysis
do_boot_wt_rcb_default <- function() {
  WT_PH5_BT_IGNOR_RCB = get_boot_ci(
    design = WT_PH5_ELEC_RCB$design,
    oc = WT_PH5_OC_RCB,
    meta = WT_PH5_ELEC_RCB$meta,
    rdata_path = WT_PH5_BT_IGNOR_NA_PATH_RCB,
    overwrite = SHOULD_OVERWRITE
  )
  
  WT_PH7_BT_IGNOR_RCB = get_boot_ci(
    design = WT_PH7_ELEC_RCB$design,
    oc = WT_PH7_OC_RCB,
    meta = WT_PH7_ELEC_RCB$meta,
    rdata_path = WT_PH7_BT_IGNOR_NA_PATH_RCB,
    overwrite = SHOULD_OVERWRITE
  )
  
  WT_PH5_BT_IMPUT_RCB = get_boot_ci(
    design = WT_PH5_ELEC_RCB$design,
    oc = WT_PH5_OC_RCB,
    meta = WT_PH5_ELEC_RCB$meta,
    rdata_path = WT_PH5_BT_IMPUT_NA_PATH_RCB,
    replace.na = TRUE,
    overwrite = SHOULD_OVERWRITE
  )
  
  WT_PH7_BT_IMPUT_RCB = get_boot_ci(
    design = WT_PH7_ELEC_RCB$design,
    oc = WT_PH7_OC_RCB,
    meta = WT_PH7_ELEC_RCB$meta,
    rdata_path = WT_PH7_BT_IMPUT_NA_PATH_RCB,
    replace.na = TRUE,
    overwrite = SHOULD_OVERWRITE
  )
}

# Call to perform all wiltdtype boot statistics with and w/o NA replacement
do_boot_wt_rcb_all <-
  function(stat_types = c("MEDQ_DCLOSE", "MEDQ_DOPEN", "MED_DIFF", "MEA_DIFF"),
           replace.nas = c(TRUE, FALSE),
           confs = c(0.95, 0.90, 0.85),
           should_collect = TRUE) {
    base_output_dir = get_relax_c_beta_dir()
    sim_id = "wt"
    output = list()
    
    for (stat_type in stat_types) {
      for (replace.na in replace.nas) {
        for (conf in confs) {
          bt_pH5_path = get_boot_path(
            base_output_dir = base_output_dir,
            pH = 5,
            sim_id = sim_id,
            replace.na = replace.na,
            stat_type = stat_type,
            conf = conf
          )
          
          bt_pH7_path = get_boot_path(
            base_output_dir = base_output_dir,
            pH = 7,
            sim_id = sim_id,
            replace.na = replace.na,
            stat_type = stat_type,
            conf = conf
          )
          
          bt_pH5 = get_boot_ci(
            design = WT_PH5_ELEC_RCB$design,
            oc = WT_PH5_OC_RCB,
            meta = WT_PH5_ELEC_RCB$meta,
            rdata_path = bt_pH5_path,
            replace.na = replace.na,
            stat_type = stat_type,
            conf = conf
          )
          
          bt_pH7 = get_boot_ci(
            design = WT_PH7_ELEC_RCB$design,
            oc = WT_PH7_OC_RCB,
            meta = WT_PH7_ELEC_RCB$meta,
            rdata_path = bt_pH7_path,
            replace.na = replace.na,
            stat_type = stat_type,
            conf = conf
          )
          
          if (should_collect) {
            pH5_id = get_fid_from_rdata(bt_pH5_path)
            output[[pH5_id]] = list(bt = bt_pH5, path = bt_pH5_path)
            pH7_id = get_fid_from_rdata(bt_pH7_path)
            output[[pH7_id]] = list(bt = bt_pH7, path = bt_pH7_path)
          }
        }
      }
    }
    return(output)
  }

# FCR adjustment
do_boot_adjust_wt_rcb_all <- function(alpha = 0.05,
                                      stat_types = c("MED_DIFF"),
                                      should_collect = TRUE) {
  bt.ci.ls = do_boot_wt_rcb_all(
    stat_types = stat_types,
    replace.nas = TRUE,
    confs = 1.0 - alpha,
    should_collect = TRUE
  )
  
  output = list()
  
  for (fid in names(bt.ci.ls)) {
    bt.ci = bt.ci.ls[[fid]]$bt
    bt.ci_rdata_path = bt.ci.ls[[fid]]$path
    
    bt.fcr.ls = adjust_boot_ci(bt.ci = bt.ci,
                               bt.ci_rdata_path = bt.ci_rdata_path,
                               alpha = alpha)
    
    if (should_collect) {
      fid = get_fid_from_rdata(bt.fcr.ls$path)
      output[[fid]] = bt.fcr.ls$bt
    }
  }
  
  return(output)
}

# Creates bootstrap confidence interval bar plots for wild type
do_boot_ci_bar_wt_rcb_all <- function() {
  stat_type = "MED_DIFF"
  ci = BOOTSTRAP_CONF
  ci_str = as.character(as.integer(100.0 * ci))
  
  bt_ls = list(
    "5" = get_boot_ci(
      design = WT_PH5_ELEC_RCB$design,
      oc = WT_PH5_OC_RCB,
      meta = WT_PH5_ELEC_RCB$meta,
      rdata_path = WT_PH5_BT_IMPUT_NA_PATH_RCB,
      replace.na = TRUE,
      conf = ci,
      stat_type = stat_type,
      overwrite = SHOULD_OVERWRITE
    ),
    "7" = get_boot_ci(
      design = WT_PH7_ELEC_RCB$design,
      oc = WT_PH7_OC_RCB,
      meta = WT_PH7_ELEC_RCB$meta,
      rdata_path = WT_PH7_BT_IMPUT_NA_PATH_RCB,
      replace.na = TRUE,
      conf = ci,
      stat_type = stat_type,
      overwrite = SHOULD_OVERWRITE
    )
  )
  
  target_cols = c("stab", "destab")
  target_vals = c("OPEN", "CLOSE")
  target_descs = list(stab = "Stabilizing", destab = "Destabilizing")
  target_palletes = list(
    stab = list(
      # Yankees Blue
      fill_lo = "#132B43",
      # Blue Jeans
      fill_hi = "#56B1F7",
      # Mandarin
      col_err = "#E57A44"
    ),
    destab = list(
      fill_lo = "darkred",
      fill_hi = "darksalmon",
      # Pastel Purple
      col_err = "#AA9FB1"
    )
  )
  res_flts = c("ALL", "ION", "POL", "ION.POL")
  res_flt_descs = list(
    ALL = "",
    ION = "Ionic",
    POL = "Polar",
    ION.POL = "Ionic or Polar"
  )
  
  misc_flts = c("ALL", "LOOP.LOOP", "XORLOOP")
  misc_flt_descs = list(ALL = "",
                        LOOP.LOOP = "Loop-to-Loop",
                        XORLOOP = "Loop-to-Non-Loop")
  
  # Utility for generating plot according to data sets used
  do_plot <- function(bt.ci,
                      pH_str,
                      res_flt,
                      res_flt_desc,
                      misc_flt,
                      misc_flt_desc,
                      target_col,
                      target_val,
                      target_desc,
                      target_pallete,
                      img_dir) {
    ix_target = bt.ci[, target_col] == target_val
    bt.ci = bt.ci[ix_target, ]
    if (res_flt == "ION") {
      print("... keeping ionic residues")
      bt.ci = bt.ci[bt.ci$ION.ION, ]
    } else if (res_flt == "POL") {
      print("... keeping polar residues")
      bt.ci = bt.ci[bt.ci$POL.POL, ]
    } else if (res_flt == "ION.POL") {
      print("... keeping ionic and polar residues")
      logi_ion_pol = ((bt.ci$RESA.ION |
                         bt.ci$RESA.POL) &
                        (bt.ci$RESB.ION | bt.ci$RESB.POL))
      bt.ci = bt.ci[logi_ion_pol, ]
    }
    
    if (misc_flt == "LOOP.LOOP") {
      bt.ci = bt.ci[bt.ci$LOOP.LOOP, ]
    }
    else if (misc_flt == "XORLOOP") {
      logi_loop_barrel = xor((bt.ci$RESA.LOOP != "NONE"),
                             (bt.ci$RESB.LOOP != "NONE"))
      bt.ci = bt.ci[logi_loop_barrel, ]
    }
    
    if (nrow(bt.ci) <= 0) {
      print("No records found. Skipping.")
      return(NULL)
    }
    
    main = paste0(
      target_val,
      ": Median Net ",
      target_desc,
      ", pH=",
      pH_str,
      "\nBootstrap ",
      ci_str,
      "% CI"
    )
    topk = min(30, nrow(bt.ci))
    
    flt_sep = ""
    if (nzchar(res_flt_desc) && nzchar(misc_flt_desc)) {
      flt_sep = ", "
    }
    flt_pad = ""
    flt_main_sep = ""
    if (nzchar(res_flt_desc) || nzchar(misc_flt_desc)) {
      flt_main_sep = ":"
      flt_pad = " "
    }
    
    flt_base_desc = paste0(res_flt_desc,
                           flt_sep,
                           misc_flt_desc)
    flt_main_desc = paste0(flt_base_desc, flt_main_sep, flt_pad)
    flt_lab_desc = paste0(flt_pad,
                          flt_base_desc)
    
    main = paste0(main,
                  "\n(",
                  flt_main_desc,
                  "Top ",
                  topk,
                  " of ",
                  nrow(bt.ci),
                  ")")
    xlab = paste0("OmpG", flt_lab_desc, " Residue Interactions")
    ylab = paste0("Median Net ", target_desc)
    
    png_name = paste0(
      "bt.ci.",
      ci_str,
      ".",
      res_flt,
      ".",
      misc_flt,
      ".",
      target_col,
      ".",
      target_val,
      ".",
      stat_type,
      ".png"
    )
    png_path = file.path(img_dir, png_name)
    print(paste0("Generating: ", png_path))
    p = plot_boot_ci_bar(
      bt.ci = bt.ci,
      main = main,
      xlab = xlab,
      ylab = ylab,
      topk = topk,
      fill_lo = target_pallete$fill_lo,
      fill_hi = target_pallete$fill_hi,
      col_err = target_pallete$col_err
    )
    ggsave(filename = png_path, plot = p)
  }
  
  # Iterate over data set types (pH, etc.)
  for (pH_str in names(bt_ls)) {
    bt.ci = bt_ls[[pH_str]]
    
    # Create image base subdirectory
    img_base_dir = file.path(get_relax_c_beta_dir(),
                             "capt",
                             paste0("pH", pH_str),
                             "wt",
                             "plot_boot_ci_bar")
    dir.create(img_base_dir, showWarnings = FALSE)
    
    # All interactions
    # Only ion-ion
    # Only ion-polar
    # Select 'stab' or 'destab'
    for (target_col in target_cols) {
      target_desc = target_descs[[target_col]]
      target_pallete = target_palletes[[target_col]]
      # Select 'OPEN' or 'CLOSE'
      for (target_val in target_vals) {
        # Select filter 'ALL', 'ION', 'POL', 'ION.POL'
        for (res_flt in res_flts) {
          res_flt_desc = res_flt_descs[[res_flt]]
          # Select interactions "ALL", "LOOP.LOOP", "XORLOOP"
          for (misc_flt in misc_flts) {
            misc_flt_desc = misc_flt_descs[[misc_flt]]
            img_dir_name = paste0("r", res_flt, ".m", misc_flt)
            img_dir = file.path(img_base_dir, img_dir_name)
            dir.create(img_dir, showWarnings = FALSE)
            do_plot(
              bt.ci = bt.ci,
              pH_str = pH_str,
              res_flt = res_flt,
              res_flt_desc = res_flt_desc,
              misc_flt = misc_flt,
              misc_flt_desc = misc_flt_desc,
              target_col = target_col,
              target_val = target_val,
              target_desc = target_desc,
              target_pallete = target_pallete,
              img_dir = img_dir
            )
          } # end iteration over interaction filter
        } # end iteration over residue filter
      } # end iteration over target value (OPEN or CLOSE)
    } # end iteration over target column (stab or destab)
  } # end iteration over pH
}

###################################################################
# Ionic interaction correlation analysis
###################################################################

# @return Path to correlation rdata
get_cor_path <- function(base_output_dir = get_relax_c_beta_dir(),
                         pH = 5,
                         sim_id = "wt") {
  pH_str = paste0("pH", pH)
  base.name = paste0(sim_id, ".", pH_str)
  file.path(base_output_dir,
            "capt",
            pH_str,
            sim_id,
            paste0("corr.", pH_str, ".rdata"))
}

# WT RCB

WT_PH5_COR_PATH_RCB = get_cor_path(base_output_dir = get_relax_c_beta_dir(),
                                   pH = 5,
                                   sim_id = "wt")

WT_PH7_COR_PATH_RCB = get_cor_path(base_output_dir = get_relax_c_beta_dir(),
                                   pH = 7,
                                   sim_id = "wt")

# @param design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples.
# @param oc - logical vector with TRUE = open, FALSE = close.
# @param meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSe - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#   See description in load_data_utils::elec_csv_2_rdata() for full column listings
# @param rdata_path - Output rdata path
# @param overwrite - If true, will always recompute analysis and overwrite anything already
#   stored to disk. Else, will check for previously computed rdata and return the cached result.
# @return list with members
#   $pears_all, pears_open, pears_close, spear_all, spear_open, spear_close
get_cor <- function(design,
                    oc,
                    meta,
                    rdata_path,
                    overwrite = SHOULD_OVERWRITE) {
  print(paste0("Correlation using base output path: ", rdata_path))
  
  # Check if cached results exist
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping correlation as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Correlation - probably want all values present, so "impute"
  design[is.na(design)] = 0.0
  
  # Determine paths for pearson v spearman as well as open v close v all
  cor_dir = dirname(rdata_path)
  fid = get_fid_from_rdata(rdata_path)
  # Pearson
  rdata_path_pears_all = file.path(cor_dir, paste0(fid, ".pears.all.rdata"))
  rdata_path_pears_open = file.path(cor_dir, paste0(fid, ".pears.open.rdata"))
  rdata_path_pears_close = file.path(cor_dir, paste0(fid, ".pears.close.rdata"))
  # Spearman
  rdata_path_spear_all = file.path(cor_dir, paste0(fid, ".spear.all.rdata"))
  rdata_path_spear_open = file.path(cor_dir, paste0(fid, ".spear.open.rdata"))
  rdata_path_spear_close = file.path(cor_dir, paste0(fid, ".spear.close.rdata"))
  
  # Allocate output list
  out = list(
    pears_all = c(),
    pears_open = c(),
    pears_close = c(),
    spear_all = c(),
    spear_open = c(),
    spear_close = c()
  )
  
  # Filter covariates to keep only ion-ion interactions where at least one of the
  # ions is within a loop region and none are within the same loop region:
  
  # First keep only ion-ion interactions
  filt = meta$ION.ION
  # Now, make sure at least one of the ions is within a loop region
  filt = filt &
    ((meta$RESA.LOOP != "NONE") | (meta$RESB.LOOP != "NONE"))
  # Finally, remove any ion-ion interactions where both are in same loop
  filt = filt & (meta$RESA.LOOP != meta$RESB.LOOP)
  
  # Apply filter
  meta.filt = meta[filt, ]
  interaction_names = rownames(meta.filt)
  design = design[, interaction_names]
  
  # Total number of interaction pairs
  n.ints = ncol(design)
  print(paste0("...After filtering, ", n.ints, " interaction pairs remain."))
  
  # Create open design matrix
  design.open = design[oc, ]
  # Create closed design matrix
  design.close = design[!oc, ]
  
  # Helper method to perform correlation and save result
  do_cor <- function(d, meth, rdpath) {
    o = cor(x = d, method = meth)
    save_rdata(data = o, rdata_path = rdpath)
    write.csv(x = o, file = get_csv_path_from_rdata(rdpath))
    return(o)
  }
  
  out$pears_all = do_cor(design, "pearson", rdata_path_pears_all)
  out$pears_open = do_cor(design.open, "pearson", rdata_path_pears_open)
  out$pears_close = do_cor(design.close, "pearson", rdata_path_pears_close)
  
  out$spear_all = do_cor(design, "spearman", rdata_path_spear_all)
  out$spear_open = do_cor(design.open, "spearman", rdata_path_spear_open)
  out$spear_close = do_cor(design.close, "spearman", rdata_path_spear_close)
  
  save_rdata(data = out, rdata_path = rdata_path)
  return(out)
}

##############################
# Correlation - Examples

# Example wildtype analysis
do_cor_wt_rcb <- function() {
  WT_PH5_COR_RCB = get_cor(
    design = WT_PH5_ELEC_RCB$design,
    oc = WT_PH5_OC_RCB,
    meta = WT_PH5_ELEC_RCB$meta,
    rdata_path = WT_PH5_COR_PATH_RCB,
    overwrite = SHOULD_OVERWRITE
  )
  
  WT_PH7_COR_RCB = get_cor(
    design = WT_PH7_ELEC_RCB$design,
    oc = WT_PH7_OC_RCB,
    meta = WT_PH7_ELEC_RCB$meta,
    rdata_path = WT_PH7_COR_PATH_RCB,
    overwrite = SHOULD_OVERWRITE
  )
}

###################################################################
# SRS - Single Residue Scoring
###################################################################

# @return Corresponding SRS rdata path from elec path
get_srs_path_from_elec <-
  function(elec_path, erank = MAX_ENERGY_RANK) {
    d = dirname(elec_path)
    return(file.path(d,
                     paste0("srs.", erank, ".rdata")))
  }

# @return list with elements
#   $design - a matrix where each column is the residue score for each sample
#   $meta - meta information associated with each atom such as
#     loop identifier, summary stats, and open vs close favorability
get_srs <- function(design,
                    oc,
                    meta,
                    rdata_path,
                    overwrite = SHOULD_OVERWRITE) {
  print(paste0("Single residue scoring for: ", rdata_path))
  
  # Check if cached results exist
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping single residue scoring as rdata already exists: ",
      rdata_path
    ))
    return(load_rdata(rdata_path))
  }
  
  # Replace any missing values with 0
  design[is.na(design)] = 0.0
  
  # Pre-mult by 0.5 to address double counting
  design = 0.5 * design
  
  # Create mapping from residue name to column index
  res_ids = sort(as.numeric(unique(unlist(
    strsplit(colnames(design), ".", fixed = TRUE)
  ))))
  names(res_ids) = as.character(res_ids)
  n_res = length(res_ids)
  res_ids_to_ix = rank(res_ids)
  names(res_ids_to_ix) = names(res_ids)
  
  # Create mapping from residue pair to individual residues
  pair_id_to_res_ids = do.call(rbind, strsplit(colnames(design), ".", fixed = TRUE))
  rownames(pair_id_to_res_ids) = colnames(design)
  colnames(pair_id_to_res_ids) = c("RESA.NO", "RESB.NO")
  
  # Allocate output design matrix
  out.design = matrix(data = 0.0,
                      nrow = nrow(design),
                      ncol = n_res)
  colnames(out.design) = names(res_ids)
  rownames(out.design) = rownames(design)
  
  # Allocate output meta data.frame
  empty.na.vec = rep(NA, n_res)
  empty.z.vec = rep(0, n_res)
  out.meta = data.frame(
    RES.NO = as.integer(res_ids),
    RES.NAME = empty.na.vec,
    RES.LOOP = empty.na.vec,
    N = empty.z.vec,
    N.OPEN = empty.z.vec,
    N.CLOSE = empty.z.vec,
    RES.ION = empty.na.vec,
    RES.POL = empty.na.vec,
    MED.OPEN = empty.na.vec,
    MED.CLOSE = empty.na.vec,
    MIN.OPEN = empty.na.vec,
    Q1.OPEN = empty.na.vec,
    MEA.OPEN = empty.na.vec,
    Q3.OPEN = empty.na.vec,
    MAX.OPEN = empty.na.vec,
    MIN.CLOSE = empty.na.vec,
    Q1.CLOSE = empty.na.vec,
    MEA.CLOSE = empty.na.vec,
    Q3.CLOSE = empty.na.vec,
    MAX.CLOSE = empty.na.vec,
    FAVORS = rep("NONE", n_res),
    stringsAsFactors = FALSE
  )
  rownames(out.meta) = names(res_ids)
  
  # Iterate over interaction pairs
  stopifnot(nrow(pair_id_to_res_ids) == ncol(design))
  for (i_pair in 1:nrow(pair_id_to_res_ids)) {
    # Design contribution
    observed_vals = design[, i_pair]
    target_a_nm = pair_id_to_res_ids[i_pair, 1]
    target_b_nm = pair_id_to_res_ids[i_pair, 2]
    out.design[, target_a_nm] = out.design[, target_a_nm] + observed_vals
    out.design[, target_b_nm] = out.design[, target_b_nm] + observed_vals
    
    # Meta contribution
    target_a_ix = res_ids_to_ix[target_a_nm]
    target_b_ix = res_ids_to_ix[target_b_nm]
    pair_nm = colnames(design)[i_pair]
    
    # Meta: N, N.OPEN, N.CLOSE (observed counts)
    meta.N = meta[pair_nm, "N"]
    meta.N.OPEN = meta[pair_nm, "N.OPEN"]
    meta.N.CLOSE = meta[pair_nm, "N.CLOSE"]
    out.meta$N[target_a_ix] = max(out.meta$N[target_a_ix], meta.N)
    out.meta$N.OPEN[target_a_ix] = max(out.meta$N.OPEN[target_a_ix], meta.N.OPEN)
    out.meta$N.CLOSE[target_a_ix] = max(out.meta$N.OPEN[target_a_ix], meta.N.CLOSE)
    
    out.meta$N[target_b_ix] = max(out.meta$N[target_b_ix], meta.N)
    out.meta$N.OPEN[target_b_ix] = max(out.meta$N.OPEN[target_b_ix], meta.N.OPEN)
    out.meta$N.CLOSE[target_b_ix] = max(out.meta$N.OPEN[target_b_ix], meta.N.CLOSE)
    
    # Meta: RES.NAME, RES.LOOP, RES.ION, RES.POL
    # @TODO - This really onlys need to be set once
    out.meta$RES.NAME[target_a_ix] = meta[pair_nm, "RESA.NAME"]
    out.meta$RES.LOOP[target_a_ix] = meta[pair_nm, "RESA.LOOP"]
    out.meta$RES.ION[target_a_ix] = meta[pair_nm, "RESA.ION"]
    out.meta$RES.POL[target_a_ix] = meta[pair_nm, "RESA.POL"]
    
    out.meta$RES.NAME[target_b_ix] = meta[pair_nm, "RESB.NAME"]
    out.meta$RES.LOOP[target_b_ix] = meta[pair_nm, "RESB.LOOP"]
    out.meta$RES.ION[target_b_ix] = meta[pair_nm, "RESB.ION"]
    out.meta$RES.POL[target_b_ix] = meta[pair_nm, "RESB.POL"]
  }
  
  # Now, compute summary information
  
  # Helper method: formats a summary element as a vector
  # @param s - a summary object
  # @ix - the element to extract
  get_summary_feat <- function(s, ix) {
    out = strsplit(x = s[ix,],
                   split = ":",
                   fixed = TRUE)
    out = as.data.frame(out, stringsAsFactors = FALSE)
    # Suppress "Warning: NAs introduced by coercion"
    out = suppressWarnings(as.numeric(out[2,]))
    return(out)
  }
  
  # Helper method: formats summary as a data.frame
  get_summary <- function(d, suffix) {
    IX_MIN = 1
    IX_Q1 = 2
    IX_MED = 3
    IX_MEA = 4
    IX_Q3 = 5
    IX_MAX = 6
    s = summary(d)
    MIN = get_summary_feat(s, IX_MIN)
    Q1 = get_summary_feat(s, IX_Q1)
    MED = get_summary_feat(s, IX_MED)
    MEA = get_summary_feat(s, IX_MEA)
    Q3 = get_summary_feat(s, IX_Q3)
    MAX = get_summary_feat(s, IX_MAX)
    s_df = data.frame(
      MIN = MIN,
      Q1 = Q1,
      MED = MED,
      MEA = MEA,
      Q3 = Q3,
      MAX = MAX,
      stringsAsFactors = FALSE
    )
    colnames(s_df) = paste0(colnames(s_df), suffix)
    return(s_df)
  }
  
  # Compute and format open and close summaries
  s_df.open = get_summary(out.design[oc,], ".OPEN")
  s_df.close = get_summary(out.design[!oc,], ".CLOSE")
  
  # Copy summary information over
  out.meta[, colnames(s_df.open)] = s_df.open[, ]
  out.meta[, colnames(s_df.close)] = s_df.close[, ]
  
  # Compute favorability based on median score
  favors.open = out.meta$MED.OPEN < out.meta$MED.CLOSE
  favors.close = out.meta$MED.CLOSE < out.meta$MED.OPEN
  out.meta$FAVORS[favors.open] = "OPEN"
  out.meta$FAVORS[favors.close] = "CLOSE"
  
  # Combine into list structure
  data = list()
  data[["design"]] = out.design
  data[["meta"]] = out.meta
  
  # Save results
  save_rdata(data, rdata_path)
  design_csv_path = get_csv_path_from_rdata(insert_attrib_in_rdata(rdata_path, "design"))
  meta_csv_path = get_csv_path_from_rdata(insert_attrib_in_rdata(rdata_path, "meta"))
  write.csv(x = out.design, file = design_csv_path)
  write.csv(x = out.meta, file = meta_csv_path)
  return(data)
}

##############################
# SRS - Examples

###### Generate raw data ######

# Example wildtype analysis
do_srs_wt_rcb <- function() {
  WT_PH5_SRS_RCB = get_srs(
    design = WT_PH5_ELEC_RCB$design,
    oc = WT_PH5_OC_RCB,
    meta = WT_PH5_ELEC_RCB$meta,
    rdata_path = get_srs_path_from_elec(WT_PH5_ELEC_CSV_PATH_RCB),
    overwrite = SHOULD_OVERWRITE
  )
  
  WT_PH7_SRS_RCB = get_srs(
    design = WT_PH7_ELEC_RCB$design,
    oc = WT_PH7_OC_RCB,
    meta = WT_PH7_ELEC_RCB$meta,
    rdata_path = get_srs_path_from_elec(WT_PH7_ELEC_CSV_PATH_RCB),
    overwrite = SHOULD_OVERWRITE
  )
  
  return(list(WT_PH5_SRS_RCB = WT_PH5_SRS_RCB,
              WT_PH7_SRS_RCB = WT_PH7_SRS_RCB))
}

###### Random forest ######

# Default single residue score random forest version
RF_SRS_VERSION = RF_VERSION_IN_BAG

# Default single residue score keep setting for random forest
# Possible values are:
#     SRS.ALL - all residues kept
#     SRS.ION - only ionic residues kept
#     SRS.POLAR - only polar residues kept
#     SRS.ION.POLAR - only ionic or polar residues kept
#     SRS.LOOP.ION - only ionic residues within loops are kept
#     SRS.LOOP.ION.POLAR - only ionic or polar residues within loops are kept
# Can be vector like so: c("SRS.ALL", "SRS.ION", "SRS.POLAR") and methods
# will correctly either use all args or only take the first arg.
RF_SRS_KEEP = c("SRS.ION.POLAR", "SRS.LOOP.ION.POLAR", "SRS.LOOP.ION")

# Example single residue score wild type random forest analysis
# @param keep.ls - specifies the filter settings (supports multiple elements)
# @param pH.ls - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_srs_rf_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                             pH.ls = c(5, 7),
                             rf_version = RF_SRS_VERSION) {
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  # List for storing results
  results = list()
  
  # Map from pH to open-close logicals
  oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
  
  # Base output directory
  base_output_dir = get_relax_c_beta_dir()
  
  for (pH in pH.ls) {
    pH_str = as.character(pH)
    # Skip missing pH levels
    if ((!(pH_str %in% names(oc_by_pH)))) {
      next
    }
    
    oc = oc_by_pH[[pH_str]]
    srs_id = paste0("WT_PH", pH_str, "_SRS_RCB")
    elec = srs_ls[[srs_id]]
    
    for (keep in keep.ls) {
      rf_id = paste0("WT_PH", pH, "_SRS_RF_", keep, "_RCB")
      # Determine random forest rdata path
      rdata_path = get_rf_path(
        base_output_dir = base_output_dir,
        pH = pH,
        sim_id = "wt",
        keep = keep,
        rf_version = rf_version
      )
      # Generate random forest!
      results[[rf_id]] = get_rf(
        oc = oc,
        elec = elec,
        rdata_path = rdata_path,
        keep = keep,
        rf_version = rf_version
      )
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  # Return results!
  return(results)
}

# Example single residue score wild type probability correlation analysis
# @param keep.ls - specifies the filter settings (supports multiple elements)
# @param pH.ls - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_srs_rf_prob_cor_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                                      pH.ls = c(5, 7),
                                      rf_version = RF_SRS_VERSION) {
  # Obtain random forest fits
  rf.fits = do_srs_rf_wt_rcb(keep.ls, pH.ls, rf_version)
  
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  # Allocate output list
  results = list()
  
  for (pH in pH.ls) {
    # Identifier for SRS data
    srs_id = paste0("WT_PH", pH, "_SRS_RCB")
    if (!(srs_id %in% names(srs_ls))) {
      print(paste0("Warning: missing SRS identifier: ", srs_id))
      print("...Skipping")
      next
    }
    
    elec = srs_ls[[srs_id]]
    
    for (keep in keep.ls) {
      rfpcor_res_id = paste0("WT_PH", pH, "_SRS_RFPCOR_", keep, "_RCB")
      rf_res_id = paste0("WT_PH", pH, "_SRS_RF_", keep, "_RCB")
      
      stopifnot(rf_res_id %in% names(rf.fits))
      
      # Determine random forest rdata path
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep,
        rf_version = rf_version
      )
      
      results[[rfpcor_res_id]] = get_rf_prob_cor(
        rf.fit = rf.fits[[rf_res_id]],
        rdata_path = get_rf_prob_cor_path_from_rf(rf_rdata_path),
        elec = elec,
        keep = keep,
        rf_version = rf_version
      )
      
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Example single residue score wild type forest floor analysis
# @param keep.ls - specifies the filter settings (supports multiple elements)
# @param pH.ls - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_srs_ff_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                             pH.ls = c(5, 7),
                             rf_version = RF_SRS_VERSION) {
  # Obtain random forest fits
  rf.fits = do_srs_rf_wt_rcb(keep.ls, pH.ls, rf_version)
  
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  # Allocate output list
  results = list()
  
  for (pH in pH.ls) {
    # Identifier for SRS data
    srs_id = paste0("WT_PH", pH, "_SRS_RCB")
    if (!(srs_id %in% names(srs_ls))) {
      print(paste0("Warning: missing SRS identifier: ", srs_id))
      print("...Skipping")
      next
    }
    
    elec = srs_ls[[srs_id]]
    
    for (keep in keep.ls) {
      ff_res_id = paste0("WT_PH", pH, "_SRS_FF_", keep, "_RCB")
      rf_res_id = paste0("WT_PH", pH, "_SRS_RF_", keep, "_RCB")
      
      stopifnot(rf_res_id %in% names(rf.fits))
      
      # Determine random forest rdata path
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep,
        rf_version = rf_version
      )
      
      results[[ff_res_id]] = get_ff(
        rf.fit = rf.fits[[rf_res_id]],
        rdata_path = get_ff_path_from_rf(rf_rdata_path),
        elec = elec,
        keep = keep,
        rf_version = rf_version
      )
      
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Example single residue score wild type forest floor analysis
# @param keep.ls - specifies the filter settings (supports multiple elements)
# @param pH.ls - specifies which pH levels to keep (only 5 or 7 supported)
# @return list with results for each keep and pH combo
do_srs_ff_cor_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                                 pH.ls = c(5, 7),
                                 rf_version = RF_SRS_VERSION) {
  # Obtain random forest fits
  ff_ls = do_srs_ff_wt_rcb(keep.ls, pH.ls, rf_version)
  
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  # Map from pH to open-close logicals
  oc_by_pH = list("5" = WT_PH5_OC_RCB, "7" = WT_PH7_OC_RCB)
  
  # Allocate output list
  results = list()
  
  for (pH in pH.ls) {
    # Identifier for SRS data
    srs_id = paste0("WT_PH", pH, "_SRS_RCB")
    if (!(srs_id %in% names(srs_ls))) {
      print(paste0("Warning: missing SRS identifier: ", srs_id))
      print("...Skipping")
      next
    }
    
    oc = oc_by_pH[[as.character(pH)]]
    elec = srs_ls[[srs_id]]
    
    for (keep in keep.ls) {
      ff_cor_res_id = paste0("WT_PH", pH, "_SRS_FFCOR_", keep, "_RCB")
      ff_res_id = paste0("WT_PH", pH, "_SRS_FF_", keep, "_RCB")
      
      stopifnot(ff_res_id %in% names(ff_ls))
      
      # Determine output path
      rf_rdata_path = get_rf_path(
        base_output_dir = get_relax_c_beta_dir(),
        pH = pH,
        sim_id = "wt",
        keep = keep,
        rf_version = rf_version
      )
      ff_rdata_path = get_ff_path_from_rf(rf_rdata_path)
      ff_cor_rdata_path = get_ff_cor_path_from_ff(ff_rdata_path)
      
      results[[ff_cor_res_id]] = get_ff_cor(
        ff = ff_ls[[ff_res_id]],
        rdata_path = ff_cor_rdata_path,
        elec = elec,
        oc = oc,
        keep = keep,
        rf_version = rf_version
      )
    } # end iteration over keep.ls
  } # end iteration over pH.ls
  
  return(results)
}

# Utility to obtain a human readable srs name
# Converts <resno> to <res1LetterCode><resno>
# Eg. "68" becomes "R68"
get_human_srs_feature_name <- function(df) {
  res.code = AA3to1[df$RES.NAME]
  res.no = df$RES.NO
  feat.names = paste0(res.code, res.no)
  names(feat.names) = rownames(df)
  return(feat.names)
}

# Generate feature contribution bar plots
do_srs_ff_bar_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                                 pH.ls = c(5, 7),
                                 rf_version = RF_SRS_VERSION) {
  # Obtain forest floor correlation data
  ff_cor_ls = do_srs_ff_cor_wt_rcb(keep.ls = keep.ls,
                                   pH.ls = pH.ls,
                                   rf_version = rf_version)
  
  for (pH in pH.ls) {
    # Create image subdirectory
    img_dir = file.path(get_relax_c_beta_dir(),
                        "capt",
                        paste0("pH", pH),
                        "wt",
                        "plot_srs_ff_bar")
    dir.create(img_dir, showWarnings = FALSE)
    
    for (keep in keep.ls) {
      ff_cor_res_id = paste0("WT_PH", pH, "_SRS_FFCOR_", keep, "_RCB")
      # Skip missing pH levels
      if (!(ff_cor_res_id %in% names(ff_cor_ls))) {
        print(paste0("Identifier: ", ff_cor_res_id, " not found."))
        next
      }
      
      ff_cor = ff_cor_ls[[ff_cor_res_id]]
      
      main = paste0("Random Forest Most Influential Residues, pH=", pH)
      numvar = nrow(ff_cor)
      topk = min(35, numvar)
      main = paste0(main, "\n(Top ", topk, " of ", numvar, ")")
      p = plot_ff_bar(
        ff_cor = ff_cor,
        feat_name_fnc = get_human_srs_feature_name,
        main = main,
        xlab = "OmpG Residue",
        ylab = "Influence Score",
        topk = topk
      )
      
      png_name = paste0("ff_srs_bar.pH", pH, ".", keep, ".", rf_version, ".png")
      png_path = file.path(img_dir, png_name)
      print(paste0("Generating: ", png_path))
      ggsave(filename = png_path, plot = p)
    } # end iteration over keep
  } # end iteration over pH
}

# Calls varImpPlot for each setting
do_srs_rf_importance_wt_rcb <- function(keep.ls = RF_SRS_KEEP,
                                        pH.ls = c(5, 7),
                                        rf_version = RF_SRS_VERSION) {
  # Obtain random forest trained models
  rf.fits = do_srs_rf_wt_rcb(keep.ls, pH.ls, rf_version)
  
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  for (pH in pH.ls) {
    # Identifier for SRS data
    srs_id = paste0("WT_PH", pH, "_SRS_RCB")
    if (!(srs_id %in% names(srs_ls))) {
      print(paste0("Warning: missing SRS identifier: ", srs_id))
      print("...Skipping")
      next
    }
    
    elec = srs_ls[[srs_id]]
    
    # Create image subdirectory
    img_dir = file.path(get_relax_c_beta_dir(),
                        "capt",
                        paste0("pH", pH),
                        "wt",
                        "plot_srs_rf_imp")
    dir.create(img_dir, showWarnings = FALSE)
    
    meta = elec$meta
    
    for (keep in keep.ls) {
      rf_res_id = paste0("WT_PH", pH, "_SRS_RF_", keep, "_RCB")
      rf.fit = rf.fits[[rf_res_id]]
      
      main = paste0("Random Forest Variable Importance, pH=", pH)
      numvar = nrow(rf.fit$importance)
      topk = min(30, numvar)
      main = paste0(main, "\n(Top ", topk, " of ", numvar, ")")
      
      png_name = paste0("rf_srs_imp.pH", pH, ".", keep, ".", rf_version, ".png")
      png_path = file.path(img_dir, png_name)
      
      print(paste0("Generating: ", png_path))
      png(filename = png_path,
          width = 1024,
          height = 1024)
      plot_rf_importance(
        rf.fit = rf.fit,
        feat_name_fnc = get_human_srs_feature_name,
        meta = meta,
        main = main,
        topk = topk
      )
      dev.off()
    } # end iteration over keep.ls
  } # end iteration over pH.ls
}

###### Bootstrap ######

# Default bootstrap stat for SRS
BOOTSTRAP_SRS_STAT_TYPE = "MED_DIFF_SRS"

# Default bootstrap base output dir for SRS
BOOTSTRAP_SRS_BASE_OUTPUT_DIR = get_relax_c_beta_dir()

# Wildtype bootstrap output path for SRS, pH 5
WT_PH5_SRS_BT_RCB_RDATA_PATH = get_boot_path(
  base_output_dir = BOOTSTRAP_SRS_BASE_OUTPUT_DIR,
  pH = 5,
  sim_id = "wt",
  stat_type = BOOTSTRAP_SRS_STAT_TYPE
)

# Wildtype bootstrap output path for SRS, pH 7
WT_PH7_SRS_BT_RCB_RDATA_PATH = get_boot_path(
  base_output_dir = BOOTSTRAP_SRS_BASE_OUTPUT_DIR,
  pH = 7,
  sim_id = "wt",
  stat_type = BOOTSTRAP_SRS_STAT_TYPE
)

# Example wildtype bootstrap CI analysis
do_srs_boot_wt_rcb <- function() {
  # Obtain single residue scores
  srs_ls = do_srs_wt_rcb()
  
  WT_PH5_SRS_BT_RCB = get_boot_ci(
    design = srs_ls[["WT_PH5_SRS_RCB"]]$design,
    oc = WT_PH5_OC_RCB,
    meta = srs_ls[["WT_PH5_SRS_RCB"]]$meta,
    rdata_path = WT_PH5_SRS_BT_RCB_RDATA_PATH,
    stat_type = BOOTSTRAP_SRS_STAT_TYPE
  )
  
  WT_PH7_SRS_BT_RCB = get_boot_ci(
    design = srs_ls[["WT_PH7_SRS_RCB"]]$design,
    oc = WT_PH7_OC_RCB,
    meta = srs_ls[["WT_PH7_SRS_RCB"]]$meta,
    rdata_path = WT_PH7_SRS_BT_RCB_RDATA_PATH,
    stat_type = BOOTSTRAP_SRS_STAT_TYPE
  )
  
  return(list("5" = WT_PH5_SRS_BT_RCB, "7" = WT_PH7_SRS_BT_RCB))
}

# FCR adjustment
do_srs_boot_adjust_wt_rcb_all <- function(should_collect = TRUE) {
  alpha = 0.05
  stat_type = BOOTSTRAP_SRS_STAT_TYPE
  bt.ci.ls = do_srs_boot_wt_rcb()
  bt.ci_rdata_paths = list("5" = WT_PH5_SRS_BT_RCB_RDATA_PATH,
                           "7" = WT_PH7_SRS_BT_RCB_RDATA_PATH)
  
  output = list()
  
  for (pH_str in names(bt.ci.ls)) {
    bt.ci = bt.ci.ls[[pH_str]]
    bt.ci_rdata_path = bt.ci_rdata_paths[[pH_str]]
    
    bt.fcr.ls = adjust_boot_ci(
      bt.ci = bt.ci,
      bt.ci_rdata_path = bt.ci_rdata_path,
      alpha = alpha,
      stat_type = stat_type,
      stb = boot_srs_median_stab,
      flt = boot_srs_filt,
      flt_desc = "no.flt"
    )
    
    if (should_collect) {
      fid = get_fid_from_rdata(bt.fcr.ls$path)
      output[[fid]] = bt.fcr.ls$bt
    }
  }
  
  return(output)
}

# Creates bootstrap confidence interval bar plots for wild type
do_srs_boot_ci_bar_wt_rcb_all <- function() {
  stat_type = "MED_DIFF"
  ci = BOOTSTRAP_CONF
  ci_str = as.character(as.integer(100.0 * ci))
  
  bt_ls = do_srs_boot_wt_rcb()
  
  target_cols = c("stab", "destab")
  target_vals = c("OPEN", "CLOSE")
  target_descs = list(stab = "Stabilizing", destab = "Destabilizing")
  target_palletes = list(
    stab = list(
      # Yankees Blue
      fill_lo = "#132B43",
      # Blue Jeans
      fill_hi = "#56B1F7",
      # Mandarin
      col_err = "#E57A44"
    ),
    destab = list(
      fill_lo = "darkred",
      fill_hi = "darksalmon",
      # Pastel Purple
      col_err = "#AA9FB1"
    )
  )
  res_flts = c("ALL", "ION", "POL", "ION.POL")
  res_flt_descs = list(
    ALL = "",
    ION = "Ionic",
    POL = "Polar",
    ION.POL = "Ionic or Polar"
  )
  misc_flts = c("ALL", "LOOP")
  misc_flt_descs = list(ALL = "",
                        LOOP = "Loop")
  
  # Utility for generating plot according to data sets used
  do_plot <- function(bt.ci,
                      pH_str,
                      res_flt,
                      res_flt_desc,
                      misc_flt,
                      misc_flt_desc,
                      target_col,
                      target_val,
                      target_desc,
                      target_pallete,
                      img_dir) {
    ix_target = bt.ci[, target_col] == target_val
    bt.ci = bt.ci[ix_target, ]
    if (res_flt == "ION") {
      print("... keeping ionic residues")
      bt.ci = bt.ci[bt.ci$RES.ION, ]
    } else if (res_flt == "POL") {
      print("... keeping polar residues")
      bt.ci = bt.ci[bt.ci$RES.POL, ]
    } else if (res_flt == "ION.POL") {
      print("... keeping ionic and polar residues")
      logi_ion_pol = bt.ci$RES.ION | bt.ci$RES.POL
      bt.ci = bt.ci[logi_ion_pol, ]
    }
    
    if (misc_flt == "LOOP") {
      logi_loop = (bt.ci$RES.LOOP != "NONE")
      bt.ci = bt.ci[logi_loop, ]
    }
    
    if (nrow(bt.ci) <= 0) {
      print("No records found. Skipping.")
      return(NULL)
    }
    
    main = paste0(
      target_val,
      ": Median Net ",
      target_desc,
      ", pH=",
      pH_str,
      "\nBootstrap ",
      ci_str,
      "% CI"
    )
    topk = min(30, nrow(bt.ci))
    flt_sep = ""
    if (nzchar(res_flt_desc) && nzchar(misc_flt_desc)) {
      flt_sep = ", "
    }
    flt_pad = ""
    flt_main_sep = ""
    if (nzchar(res_flt_desc) || nzchar(misc_flt_desc)) {
      flt_main_sep = ":"
      flt_pad = " "
    }
    
    flt_base_desc = paste0(res_flt_desc,
                           flt_sep,
                           misc_flt_desc)
    flt_main_desc = paste0(flt_base_desc, flt_main_sep, flt_pad)
    flt_lab_desc = paste0(flt_pad,
                          flt_base_desc)
    
    main = paste0(main,
                  "\n(",
                  flt_main_desc,
                  "Top ",
                  topk,
                  " of ",
                  nrow(bt.ci),
                  ")")
    xlab = paste0("OmpG", flt_lab_desc, " Residues")
    ylab = paste0("Median Net ", target_desc)
    
    png_name = paste0(
      "bt.ci.",
      ci_str,
      ".",
      res_flt,
      ".",
      misc_flt,
      ".",
      target_col,
      ".",
      target_val,
      ".",
      stat_type,
      ".png"
    )
    png_path = file.path(img_dir, png_name)
    print(paste0("Generating: ", png_path))
    p = plot_boot_ci_bar(
      bt.ci = bt.ci,
      feat_name_fnc = get_human_srs_feature_name,
      main = main,
      xlab = xlab,
      ylab = ylab,
      topk = topk,
      fill_lo = target_pallete$fill_lo,
      fill_hi = target_pallete$fill_hi,
      col_err = target_pallete$col_err
    )
    ggsave(filename = png_path, plot = p)
  }
  
  # Iterate over data set types (pH, etc.)
  for (pH_str in names(bt_ls)) {
    bt.ci = bt_ls[[pH_str]]
    
    # Create image subdirectory
    img_base_dir = file.path(
      get_relax_c_beta_dir(),
      "capt",
      paste0("pH", pH_str),
      "wt",
      "plot_srs_boot_ci_bar"
    )
    dir.create(img_base_dir, showWarnings = FALSE)
    
    # All interactions
    # Only ion-ion
    # Only ion-polar
    # Select 'stab' or 'destab'
    for (target_col in target_cols) {
      target_desc = target_descs[[target_col]]
      target_pallete = target_palletes[[target_col]]
      # Select 'OPEN' or 'CLOSE'
      for (target_val in target_vals) {
        # Select filter 'ALL', 'ION', 'POL', 'ION.POL'
        for (res_flt in res_flts) {
          res_flt_desc = res_flt_descs[[res_flt]]
          # Select interactions "ALL", "LOOP.BARREL"
          for (misc_flt in misc_flts) {
            misc_flt_desc = misc_flt_descs[[misc_flt]]
            img_dir_name = paste0("r", res_flt, ".m", misc_flt)
            img_dir = file.path(img_base_dir, img_dir_name)
            dir.create(img_dir, showWarnings = FALSE)
            do_plot(
              bt.ci = bt.ci,
              pH_str = pH_str,
              res_flt = res_flt,
              res_flt_desc = res_flt_desc,
              misc_flt = misc_flt,
              misc_flt_desc = misc_flt_desc,
              target_col = target_col,
              target_val = target_val,
              target_desc = target_desc,
              target_pallete = target_pallete,
              img_dir = img_dir
            )
          } # end iteration over misc filter
        } # end iteration over residue filter
      } # end iteration over target value (OPEN or CLOSE)
    } # end iteration over target column (stab or destab)
  } # end iteration over pH
}
