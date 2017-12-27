###################################################################
# sasa_analysis.R
#
# Analysis relating to solvent-accesible surface area (SASA)
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

# @return root directory of project
get_root_dir <- function() {
  file.path(get_script_dir(), '..', '..')
}

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "load_data_utils.R"))

###################################################################
# Specific path utilities
###################################################################

# Files in this directory are merged 2iwv and 2iww templates,
# but c_beta carbon is *not* fixed during minimization
# @return path to relax_c_beta directory
get_relax_c_beta_dir <- function() {
  file.path(get_root_dir(), "ompg", "output", "relax_c_beta")
}

# @return path to energy capture file
get_capt_path <- function(base_output_dir = get_relax_c_beta_dir(),
                          pH = 5,
                          sim_id = "wt") {
  file.path(base_output_dir, "capt", paste0("pH", pH), sim_id, "ener.csv")
}

# @return path to CASTp open/close file
get_castp_path <- function(base_output_dir = get_relax_c_beta_dir(),
                           pH = 5,
                           sim_id = "wt",
                           probe = 2.75) {
  file.path(base_output_dir,
            "castp",
            paste0("pH", pH),
            sim_id,
            paste0("oc", probe, ".csv"))
}

# @return path to SASA base directory
get_sasa_base_dir <-
  function(base_output_dir = get_relax_c_beta_dir(),
           pH = 5,
           sim_id = "wt") {
    file.path(base_output_dir,
              "sasa",
              paste0("pH", pH),
              sim_id)
  }

# @return path to SASA contrib directory
get_sasa_contrib_dir <-
  function(sasa_base_dir = get_sasa_base_dir()) {
    file.path(sasa_base_dir, "contrib")
  }

# @return path to per-atom SASA data.frame
get_sasa_df_path <- function(sasa_base_dir = get_sasa_base_dir()) {
  file.path(sasa_base_dir, "sasa.rdata")
}

# @return path to per-ionizable-residue SASA scores data.frame
get_ion_res_sasa_scores_path <-
  function(sasa_base_dir = get_sasa_base_dir()) {
    file.path(sasa_base_dir, "sasa.ion.res.scores.rdata")
  }

###################################################################
# Global data paths
###################################################################

# Probe radius in Angstroms for open/close analysis
# 2.75 Angstroms = K
# 1.75 Angstroms = Cl
PROBE = 2.75

# WT - Relax C Beta (RCB)
WT_PH5_ENER_CSV_PATH_RCB = get_capt_path(get_relax_c_beta_dir(), 5, "wt")
WT_PH7_ENER_CSV_PATH_RCB = get_capt_path(get_relax_c_beta_dir(), 7, "wt")
WT_PH5_OC_CSV_PATH_RCB = get_castp_path(get_relax_c_beta_dir(), 5, "wt", PROBE)
WT_PH7_OC_CSV_PATH_RCB = get_castp_path(get_relax_c_beta_dir(), 7, "wt", PROBE)

###################################################################
# Load binary data
###################################################################

# Switch to toggle overwrite behavior
SHOULD_OVERWRITE = FALSE

##############################
# Energy data sets

# Capped energy rank ultimately to limit number of interactions considered
MAX_ENERGY_RANK = 5000

# WT RCB
WT_PH5_ENER_RCB = ener_csv_2_rdata(WT_PH5_ENER_CSV_PATH_RCB, MAX_ENERGY_RANK, SHOULD_OVERWRITE)
WT_PH7_ENER_RCB = ener_csv_2_rdata(WT_PH7_ENER_CSV_PATH_RCB, MAX_ENERGY_RANK, SHOULD_OVERWRITE)

##############################
# Open-close data sets

# WT RCB
WT_PH5_OC_RCB = oc_csv_2_rdata(WT_PH5_OC_CSV_PATH_RCB,
                               WT_PH5_ENER_RCB$NAME,
                               SHOULD_OVERWRITE)

WT_PH7_OC_RCB = oc_csv_2_rdata(WT_PH7_OC_CSV_PATH_RCB,
                               WT_PH7_ENER_RCB$NAME,
                               SHOULD_OVERWRITE)

###################################################################
# Munge SASA
###################################################################

# Combines all contrib files into a single data.frame
# @param sasa_contrib_dir - Directory storing individual solvent-accessible
#   surface area files
# @param rdata_path - Output path to store cached data. If cached
#   data exists, it is returned unless overwrite = TRUE
# @param overwrite - If TRUE, any cached data stored on disk is
#   overwritten
# @return data.frame with columns
#   $AtomNo - Atom number
#   $AtomName - Name of atom
#   $ResName - Name of residue
#   $ResNo - Residuen number
#   $Vert - Vert column (?)
#   $<fid>* - A column for each contrib file contain the SA column for that
#     contrib file. SA should stand for solvent accessible surface area.
get_sasa_df <- function(sasa_contrib_dir = get_sasa_contrib_dir(),
                        rdata_path = get_sasa_df_path(),
                        overwrite = SHOULD_OVERWRITE) {
  # Check for cached data
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Obtain paths to each contrib file
  contrib_fnames = list.files(path = sasa_contrib_dir,
                              pattern = "*.contrib")
  
  contrib_paths = file.path(sasa_contrib_dir, contrib_fnames)
  
  # Obtain column names to match those of the open/close file
  library(tools)
  # Strip the last two extensions
  contrib_ids = file_path_sans_ext(file_path_sans_ext(x = contrib_fnames))
  
  # Contrib format
  contrib_cols = c("RecName",
                   "AtomNo",
                   "AtomName",
                   "ResName",
                   "ResNo",
                   "Vert",
                   "SA",
                   "MS")
  num_tok = length(contrib_cols)
  
  # Output data.frame
  sasa_df = data.frame()
  
  # Process each contrib file
  for (i in 1:length(contrib_paths)) {
    contrib_path = contrib_paths[i]
    print(paste0(i, ": Processing ", contrib_path))
    # Read all lines into character vector
    contrib_lines = readLines(con = contrib_path)
    # Strip comment lines
    first_char = substr(x = contrib_lines,
                        start = 1,
                        stop = 1)
    is_comment = (first_char == "#")
    contrib_lines = contrib_lines[!is_comment]
    # Split lines by whitespace
    contrib_tokens = strsplit(x = contrib_lines, split = "\\s+")
    # Make sure all records have same length
    stopifnot(all(sapply(X = contrib_tokens, FUN = length) == num_tok))
    # Convert to character matrix
    contrib_mat = matrix(unlist(contrib_tokens),
                         ncol = num_tok,
                         byrow = TRUE)
    colnames(contrib_mat) = contrib_cols
    # Determine if we need to generate meta data columns
    if (0 == length(sasa_df)) {
      sasa_df = as.data.frame(contrib_mat[, 2:(num_tok - 2)], stringsAsFactors =
                                FALSE)
      sasa_df$AtomNo = as.integer(sasa_df$AtomNo)
      sasa_df$ResNo = as.integer(sasa_df$ResNo)
      sasa_df$Vert = as.integer(sasa_df$Vert)
    }
    # Append solvent accessible surface area column
    stopifnot(nrow(sasa_df) == nrow(contrib_mat))
    sasa_col = contrib_ids[i]
    sasa_df[, sasa_col] = as.numeric(contrib_mat[, "SA"])
  }
  
  # Save results
  save_rdata(sasa_df, rdata_path)
  write.csv(x = sasa_df,
            quote = FALSE,
            file = get_csv_path_from_rdata(rdata_path))
  return (sasa_df)
}

###################################################################
# Residue SASA scores
###################################################################

# Computes aggregate SASA scores at each residue. Scores computed
# include:
#   1.) mean "functional group" SASA - only examines O and N atoms
#   2.) mean side chain SASA (excluding functional group)
#   3.) mean side chain + backbone SASA (excluding functional group)
#   4.) log ratio of 1:2
#   5.) log ratio of 1:3
# @param sasa_df - data.frame with columns
#   $AtomNo - Atom number
#   $AtomName - Name of atom
#   $ResName - Name of residue
#   $ResNo - Residue number
#   $Vert - Vert column (?)
#   $<fid>* - A column for each contrib file contain the SA column for that
#     contrib file. SA should stand for solvent accessible surface area.
# @param eps - padding to avoid log(0) and divide by 0
# @param rdata_path - path to write output residue scores
# @param overwrite - If TRUE, any cached data stored on disk is
#   overwritten
# @return list with elements
#   $Meta - data.frame with meta-info columns ResName and ResNo
#   $FgSasa - matrix with functional group sas scores
#   $ScSasa - matrix with side chain sas scores (excluding functional group)
#   $ScBbSasa - matrix with side chain + backbone sas scores (excluding functional group)
#   $LogFgToScSasa - matrix with log ratio of FgSasa to ScSas
#   $LogFgToScBbSasa - matrix with log ratio of FgSasa to ScBbSas
get_ion_res_sasa_scores <- function(sasa_df = get_sasa_df(),
                                    eps = 0.0001,
                                    rdata_path = get_ion_res_sasa_scores_path(),
                                    overwrite = SHOULD_OVERWRITE) {
  # Check for cached data
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Names of ionizable residues
  ion_res_names = c("ASP",
                    "GLU",
                    "HIS",
                    "ARG",
                    "LYS")
  
  # Names of backbone atoms
  bb_atom_names = c("C", "CA", "O", "N")
  
  # Map for names of "functional group" atoms
  # (in quotes because only looking at oxygens or nitrogens)
  fg_atom_names = list(
    "ASP" = c("OD1", "OD2"),
    "GLU" = c("OE1", "OE2"),
    "HIS" = c("ND1", "NE2"),
    "ARG" = c("NH1", "NH2"),
    "LYS" = c("NZ")
  )
  
  # Filter atoms that are not part of ionizable residues
  sasa_df = sasa_df[sasa_df$ResName %in% ion_res_names,]
  # Name for meta data columns in SASA frame
  sasa_meta_col_names = c("AtomNo", "AtomName", "ResName", "ResNo", "Vert")
  sasa_samp_col_names = colnames(sasa_df)[!(colnames(sasa_df) %in% sasa_meta_col_names)]
  sasa_samp_indices = which(colnames(sasa_df) %in% sasa_samp_col_names)
  n_samp = length(sasa_samp_col_names)
  
  # Determine residue identifiers
  res_nos = unique(sasa_df$ResNo)
  
  # Allocate meta frame
  n_res = length(res_nos)
  meta = data.frame(ResName = rep(NA, n_res),
                    ResNo = res_nos)
  
  # Utility allocates a matrix for storing results
  alloc_results_mat <- function() {
    m = matrix(data = 0.0,
               nrow = n_res,
               ncol = n_samp)
    colnames(m) = sasa_samp_col_names
    rownames(m) = res_nos
    return(m)
  }
  
  # Allocate results matrices
  # @TODO - save individual matrices and restore from disk
  results.FgSasa = alloc_results_mat()
  results.ScSasa = alloc_results_mat()
  results.ScBbSasa = alloc_results_mat()
  results.LogFgToScSasa = alloc_results_mat()
  results.LogFgToScBbSasa = alloc_results_mat()
  
  # Process each residue
  for (i in 1:length(res_nos)) {
    res_no = res_nos[i]
    res_df = sasa_df[sasa_df$ResNo == res_no,]
    res_name = res_df$ResName[1]
    meta$ResName[i] = res_name
    print(paste0(i, ": Processing residue ", res_name, "|", res_no))
    
    # Compute mean fg SASA
    is_fg_atom = res_df$AtomName %in% fg_atom_names[[res_name]]
    fg_df = res_df[is_fg_atom, sasa_samp_indices]
    fg_sasa = colMeans(fg_df)
    results.FgSasa[i, ] = fg_sasa
    
    # Compute mean sc - fg SASA
    is_bb_atom = res_df$AtomName %in% bb_atom_names
    sc_df = res_df[(!is_fg_atom) & (!is_bb_atom), sasa_samp_indices]
    sc_sasa = colMeans(sc_df)
    results.ScSasa[i, ] = sc_sasa
    
    # Compute mean sc + bb - fg SASA
    sc_bb_df = res_df[!is_fg_atom, sasa_samp_indices]
    sc_bb_sasa = colMeans(sc_bb_df)
    results.ScBbSasa[i, ] = sc_bb_sasa
    
    # Compute log [ fg : sc - fg ] SASA
    fg_eps_sasa = fg_sasa + eps
    sc_eps_sasa = sc_sasa + eps
    log_fg_to_sc_sasa = log(fg_eps_sasa / sc_eps_sasa)
    stopifnot(all(is.finite(log_fg_to_sc_sasa)))
    results.LogFgToScSasa[i, ] = log_fg_to_sc_sasa
    
    # Compute log [ fg : sc + bb - fg ] SASA
    sc_bb_eps_sasa = sc_bb_sasa + eps
    log_fg_to_sc_bb_sasa = log(fg_eps_sasa / sc_bb_eps_sasa)
    stopifnot(all(is.finite(log_fg_to_sc_bb_sasa)))
    results.LogFgToScBbSasa[i, ] = log_fg_to_sc_bb_sasa
  } # end iteration over residues
  
  # Write as list
  out.ls = list(
    Meta = meta,
    FgSasa = results.FgSasa,
    ScSasa = results.ScSasa,
    ScBbSasa = results.ScBbSasa,
    LogFgToScSasa = results.LogFgToScSasa,
    LogFgToScBbSasa = results.LogFgToScBbSasa
  )
  
  save_rdata(out.ls, rdata_path)
  
  return(out.ls)
}

###################################################################
# Bootstrap
###################################################################

##############################
# Libraries

library(boot)
library(matrixStats)
library(parallel)

##############################
# Globals

# Default number of bootstrap trials
NUM_BOOTSTRAP_TRIALS = 10000

# Default number of parallel estimates
NUM_BOOTSTRAP_CPUS = detectCores()

# Default significance level
BOOTSTRAP_ALPHA = 0.05

# Default bootstrap data type - for residue level
BOOTSTRAP_DATA_TYPE = "LogFgToScBbSasa"

##############################
# Paths

# @return path to SASA bootstrap summary data
get_sasa_boot_sum_path <-
  function(sasa_base_dir = get_sasa_base_dir(),
           data_type = BOOTSTRAP_DATA_TYPE,
           stat_type = "med_diff",
           alpha = BOOTSTRAP_ALPHA) {
    file.path(
      sasa_base_dir,
      paste0("sasa.boot.",
             data_type, ".",
             stat_type, ".",
             alpha, ".rdata")
    )
  }

# @return path to SASA bootstrap replicate data path
get_sasa_boot_rep_path <-
  function(sasa_boot_sum_path = get_sasa_boot_sum_path()) {
    return(insert_attrib_in_rdata(sasa_boot_sum_path, "reps"))
  }

##############################
# Misc utils

# Appends loop assignment to residue data.frame
# @param meta - Contains meta data about residue types
# @return data.frame with Loop column with loop labels
add_loop_labels <- function(meta) {
  # Conditionally add loop column
  if (!("Loop" %in% colnames(meta))) {
    meta = data.frame(meta, Loop = rep(NA, nrow(meta)))
  }
  
  # Make sure to assign default value
  meta$Loop = "NONE"
  
  # Utility to assign loop label to parameter region
  assign_loop <- function(meta, start_res_no, end_res_no, loop_id) {
    logi = ((meta$ResNo >= start_res_no) & (meta$ResNo <= end_res_no))
    meta$Loop[logi] = loop_id
    return(meta)
  }
  
  # Assign specific loop labels
  meta = assign_loop(meta, 18, 29, "LOOP_1")
  meta = assign_loop(meta, 54, 65, "LOOP_2")
  meta = assign_loop(meta, 97, 106, "LOOP_3")
  meta = assign_loop(meta, 177, 188, "LOOP_5")
  meta = assign_loop(meta, 217, 234, "LOOP_6")
  meta = assign_loop(meta, 259, 267, "LOOP_7")
  return(meta)
}

##############################
# Stats

# Computes median difference between open and closed
# @param obs - Matrix of observations where each row is a
#   multivariate observation
# @param i - Observation indices - corresponds to the rows
#   of data to use for computing statistic set
# @param n.open - The first n.open rows of data correspond to open
#   samples
# @return vector of median differences
sasa_med_diff_stat <- function(obs, i, n.open) {
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  obs.open = obs[i.open,]
  obs.close = obs[i.close,]
  med.open = colMedians(x = obs.open)
  med.close = colMedians(x = obs.close)
  med.diff = med.open - med.close
  return(med.diff)
}

# Computes mean difference between open and closed
# @param obs - Matrix of observations where each row is a
#   multivariate observation
# @param i - Observation indices - corresponds to the rows
#   of data to use for computing statistic set
# @param n.open - The first n.open rows of data correspond to open
#   samples
# @return vector of mean differences
sasa_mea_diff_stat <- function(obs, i, n.open) {
  stopifnot(n.open > 0)
  i.open = i[i <= n.open]
  i.close = i[i > n.open]
  stopifnot(i.open > 0)
  stopifnot(i.close > 0)
  obs.open = obs[i.open,]
  obs.close = obs[i.close,]
  mea.open = colMeans(x = obs.open)
  mea.close = colMeans(x = obs.close)
  mea.diff = mea.open - mea.close
  return(mea.diff)
}

# Determines if accessiblity is consistently greater in open vs closed state
# @param meta - data.frame contain $lo_id, $hi_id, and $acc_id columns
# @param lo_id - name of column containing lower confidence interval
# @param hi_id - name of column containing higher confidence interval
# @param acc_id - name of column to write accessiblity label
# @return meta data.frame with target accessiblity column
#   filled with "OPEN", "CLOSE", or "NONE" indicating which
#   state has greater accessibility
sasa_diff_acc <- function(meta, lo_id, hi_id, acc_id) {
  # Check if open state is more accessible
  acc_open = (meta[, lo_id] > 0)
  meta[acc_open, acc_id] = "OPEN"
  
  # Check if close state is more accessible
  acc_close = (meta[, hi_id] < 0)
  meta[acc_close, acc_id] = "CLOSE"
  return(meta)
}

##############################
# Core

# Computes bootstrap statistics on the solvent-accessible surface
# area of parameter SASA data.frame
# @param sasa_scores_mat - Aggregate residue SASA scores matrix where
#   each column is a sample and each row is a residue
# @param meta - Data frame with residue name and number for each row
#   of sasa_scores_mat
# @param oc - Open/close status of each sample column in sasa_df
# @param stat - The statistic to compute
# @param acc - Method for determining if accessibility is favored in
#   open vs close state
# @param alpha - The base significance level: (1-alpha) CI
# @param num_trials - Number of bootstrap replicates to generate
# @param num_cpus - Number of parallel CPUs to use
# @param sum_rdata_path - Path to output bootstrap summary data
# @param rep_rdata_path - Path to output bootstrap replicate data
# @param overwrite - If TRUE, any cached data stored on disk is
#   overwritten
# @return data.frame with bootstrap summary data
get_sasa_boot <-
  function(sasa_scores_mat = get_ion_res_sasa_scores()[[BOOTSTRAP_DATA_TYPE]],
           meta = get_ion_res_sasa_scores()$Meta,
           oc = WT_PH5_OC_RCB,
           stat = sasa_med_diff_stat,
           acc = sasa_diff_acc,
           alpha = BOOTSTRAP_ALPHA,
           num_trials = NUM_BOOTSTRAP_TRIALS,
           num_cpus = NUM_BOOTSTRAP_CPUS,
           sum_rdata_path = get_sasa_boot_sum_path(alpha = alpha),
           rep_rdata_path = get_sasa_boot_rep_path(sum_rdata_path),
           overwrite = SHOULD_OVERWRITE) {
    stopifnot(alpha > 0.0 && alpha < 1.0)
    
    #########################################
    # Check if cached model exists
    
    if (!overwrite && file.exists(sum_rdata_path)) {
      print(paste0(
        "Skipping bootstrap as rdata already exists: ",
        sum_rdata_path
      ))
      return(load_rdata(sum_rdata_path))
    }
    
    # Append loop labels
    meta = add_loop_labels(meta)
    
    #########################################
    # Separate open and close samples
    
    n.open = sum(oc)
    n.close = length(oc) - n.open
    names.open = names(oc)[oc]
    names.close = names(oc)[!oc]
    sasa_open = sasa_scores_mat[, names.open]
    sasa_close = sasa_scores_mat[, names.close]
    
    # Generate input bootstrap observation data
    # Need to transpose data.frame as boot() expects observations in rows
    open_t = t(sasa_open)
    close_t = t(sasa_close)
    # Place open before close to allow stratified statistic computation
    obs = rbind(open_t, close_t)
    rownames(obs) = c(names.open, names.close)
    colnames(obs) = meta$ResNo
    
    #########################################
    # Compute bootstrap replicate data
    
    # Check if cached data exists
    rep_data = c()
    if (!overwrite && file.exists(rep_rdata_path)) {
      print(paste0(
        "Skipping bootstrap replication as rdata already exists: ",
        rep_rdata_path
      ))
      rep_data = load_rdata(rep_rdata_path)
    } else {
      # Do bootstrap
      rep_data = boot(
        data = obs,
        statistic = stat,
        strata = rep(x = c(1, 2),
                     times = c(n.open,
                               n.close)),
        R = num_trials,
        parallel = "multicore",
        ncpus = num_cpus,
        n.open = n.open
      )
    }
    
    # Check everything is expected size
    n_elem = nrow(meta)
    stopifnot(n_elem == length(rep_data$t0))
    stopifnot(n_elem == ncol(rep_data$t))
    
    #########################################
    # Compute 1-alpha confidence interval
    
    meta$T0 = rep_data$t0
    conf = 1.0 - alpha
    
    # Utility for appending confidence intervals
    do_ci <- function(meta, conf, conf_id, lo_id, hi_id, acc_id) {
      meta[, conf_id] = conf
      meta[, lo_id] = 0.0
      meta[, hi_id] = 0.0
      meta[, acc_id] = "NONE"
      for (i in 1:n_elem) {
        ci = boot.ci(
          boot.out = rep_data,
          conf = conf,
          type = "basic",
          index = i
        )
        
        if (!is.null(ci)) {
          ci = ci$basic
          meta[i, lo_id] = ci[length(ci) - 1]
          meta[i, hi_id] = ci[length(ci)]
        } else {
          # http://stackoverflow.com/questions/16652852/avoid-error-in-r-boot-ci-function-when-all-values-in-sampled-set-are-equal
          # Replace with observed value if CI cannot be computed due
          # to lack of uniqueness
          meta[i, lo_id] = rep_data$t0[i]
          meta[i, hi_id] = rep_data$t0[i]
        }
      }
      
      # Determine which state is more accessible
      meta = acc(meta, lo_id, hi_id, acc_id)
      
      return(meta)
    }
    
    meta = do_ci(meta, conf, "Conf.Ci", "Lo.Ci", "Hi.Ci", "Acc.Ci")
    
    #########################################
    # Compute p-value
    
    meta$Pval = 0.0
    for (i in 1:n_elem) {
      # Munge observed stat and boot stats into single set
      reps = c(rep_data$t0[i], rep_data$t[, i])
      # Calculate p-value as min(% <= 0, % >= 0)
      n = length(reps)
      meta$Pval[i] = min(sum(reps <= 0.0) / n, sum(reps >= 0.0) / n)
    }
    
    #########################################
    # Compute FDR
    
    meta$Fdr = p.adjust(p = meta$Pval, method = "fdr")
    meta$Bon = p.adjust(p = meta$Pval, method = "bonferroni")
    
    #########################################
    # Compute FCR
    
    # Determine number of "significant" elements according to alpha
    k = sum(meta$Fdr < alpha)
    # Adjust confidence region
    conf = 1.0 - (k * alpha / n_elem)
    
    meta = do_ci(meta, conf, "Conf.Fcr", "Lo.Fcr", "Hi.Fcr", "Acc.Fcr")
    
    #########################################
    # Save data
    
    # Summary
    save_rdata(meta, sum_rdata_path)
    write.csv(x = meta,
              quote = FALSE,
              file = get_csv_path_from_rdata(sum_rdata_path))
    
    # Reps
    save_rdata(rep_data, rep_rdata_path)
    
    return(meta)
  }

###################################################################
# Shell
###################################################################

##############################
# Compute SASA df

do_sasa_df <- function(base_output_dir = get_relax_c_beta_dir(),
                       pH.ls = c(5, 7),
                       sim_id.ls = c("wt")) {
  results = list()
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      sasa_base_dir = get_sasa_base_dir(base_output_dir = base_output_dir,
                                        pH = pH,
                                        sim_id = sim_id)
      sasa_contrib_dir = get_sasa_contrib_dir(sasa_base_dir = sasa_base_dir)
      rdata_path = get_sasa_df_path(sasa_base_dir = sasa_base_dir)
      
      print(paste0("Processing ", rdata_path))
      sasa_df = get_sasa_df(sasa_contrib_dir = sasa_contrib_dir,
                            rdata_path = rdata_path)
      
      result_id = paste0(sim_id, ".", "pH", pH)
      results[[result_id]] = sasa_df
    } # end iteration over simulation identifier
  } # end iteration over pH
  
  return(results)
}

##############################
# Compute SASA residue scores

do_ion_res_sasa_scores <-
  function(base_output_dir = get_relax_c_beta_dir(),
           pH.ls = c(5, 7),
           sim_id.ls = c("wt")) {
    results = list()
    for (pH in pH.ls) {
      for (sim_id in sim_id.ls) {
        sasa_df = do_sasa_df(
          base_output_dir = base_output_dir,
          pH.ls = pH,
          sim_id.ls = sim_id
        )[[1]]
        
        sasa_base_dir = get_sasa_base_dir(base_output_dir = base_output_dir,
                                          pH = pH,
                                          sim_id = sim_id)
        
        rdata_path = get_ion_res_sasa_scores_path(sasa_base_dir = sasa_base_dir)
        
        print(paste0("Processing ", rdata_path))
        
        ion_res_sasa_scores_ls = get_ion_res_sasa_scores(sasa_df = sasa_df,
                                                         rdata_path = rdata_path)
        
        result_id = paste0(sim_id, ".", "pH", pH)
        results[[result_id]] = ion_res_sasa_scores_ls
      } # end iteration over simulation identifier
    } # end iteration over pH
    
    return(results)
  }

##############################
# Compute SASA bootstrap

do_sasa_boot <- function(base_output_dir = get_relax_c_beta_dir(),
                         pH.ls = c(5, 7),
                         sim_id.ls = c("wt"),
                         stat_type.ls = c("med_diff", "mea_diff"),
                         data_type.ls = c("FgSasa",
                                          "ScSasa",
                                          "ScBbSasa",
                                          "LogFgToScSasa",
                                          "LogFgToScBbSasa")) {
  oc_ls = list(wt.pH5 = WT_PH5_OC_RCB,
               wt.pH7 = WT_PH7_OC_RCB)
  
  stat_ls = list(med_diff = sasa_med_diff_stat,
                 mea_diff = sasa_mea_diff_stat)
  
  results = list()
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      sasa_scores_mat = do_ion_res_sasa_scores(
        base_output_dir = base_output_dir,
        pH.ls = pH,
        sim_id.ls = sim_id
      )[[1]]
      
      meta = sasa_scores_mat$Meta
      
      ident = paste0(sim_id, ".", "pH", pH)
      oc = oc_ls[[ident]]
      
      sasa_base_dir = get_sasa_base_dir(base_output_dir = base_output_dir,
                                        pH = pH,
                                        sim_id = sim_id)
      for (stat_type in stat_type.ls) {
        stat = stat_ls[[stat_type]]
        
        for (data_type in data_type.ls) {
          sum_rdata_path = get_sasa_boot_sum_path(
            sasa_base_dir = sasa_base_dir,
            data_type = data_type,
            stat_type = stat_type
          )
          
          print(paste0("Processing ", sum_rdata_path))
          
          result = get_sasa_boot(
            sasa_scores_mat = sasa_scores_mat[[data_type]],
            meta = meta,
            oc = oc,
            stat = stat,
            sum_rdata_path = sum_rdata_path
          )
          
          final_ident = paste0(ident, ".", stat_type, ".", data_type)
          results[[final_ident]] = result
          
        } # end iteration over data_type
      } # end iteration over stat_type
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

do_sasa_atom_boot <-
  function(base_output_dir = get_relax_c_beta_dir(),
           pH.ls = c(5, 7),
           sim_id.ls = c("wt"),
           stat_type.ls = c("med_diff", "mea_diff")) {
    oc_ls = list(wt.pH5 = WT_PH5_OC_RCB,
                 wt.pH7 = WT_PH7_OC_RCB)
    
    stat_ls = list(med_diff = sasa_med_diff_stat,
                   mea_diff = sasa_mea_diff_stat)
    
    results = list()
    for (pH in pH.ls) {
      for (sim_id in sim_id.ls) {
        ident = paste0(sim_id, ".", "pH", pH)
        oc = oc_ls[[ident]]
        
        sasa_df = do_sasa_df(
          base_output_dir = base_output_dir,
          pH.ls = pH,
          sim_id.ls = sim_id
        )[[1]]
        
        # Column names of samples
        sample_names = names(oc)
        # Column names of meta data (drop Vert column)
        meta_names = colnames(sasa_df)[!(colnames(sasa_df) %in% c(sample_names, "Vert"))]
        
        # Keep only rows for ionizable residues
        is.ion = sasa_df$ResName %in% c("ASP", "GLU", "HIS", "ARG", "LYS")
        
        # Extract meta frame and scores matrix
        meta = sasa_df[is.ion, meta_names]
        sasa_scores_mat = as.matrix(sasa_df[is.ion, sample_names])
        rownames(sasa_scores_mat) = meta$AtomNo
        
        sasa_base_dir = get_sasa_base_dir(base_output_dir = base_output_dir,
                                          pH = pH,
                                          sim_id = sim_id)
        
        for (stat_type in stat_type.ls) {
          stat = stat_ls[[stat_type]]
          
          sum_rdata_path = get_sasa_boot_sum_path(
            sasa_base_dir = sasa_base_dir,
            data_type = "Atom",
            stat_type = stat_type
          )
          
          print(paste0("Processing ", sum_rdata_path))
          
          result = get_sasa_boot(
            sasa_scores_mat = sasa_scores_mat,
            meta = meta,
            oc = oc,
            stat = stat,
            sum_rdata_path = sum_rdata_path
          )
          
          final_ident = paste0(ident, ".", stat_type)
          results[[final_ident]] = result
          
        } # end iteration over stat_type
      } # end iteration over sim_id
    } # end iteration over pH
    
    return(results)
  }
