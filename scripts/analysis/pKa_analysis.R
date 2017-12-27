###################################################################
# pKa_analysis.R
#
# Utlities for wrangling and analysing pKa data
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
source(file.path(get_script_dir(), "electrostat_merge_utilities.R"))

# Global wild type simulation
SIM_WT = load_sim(
  base_output_dir = DEF_BASE_OUTPUT_DIR,
  sim_prefix = "wt",
  max_energy_rank = DEF_MAX_ENERGY_RANK,
  should_load_elec = FALSE
)

###################################################################
# Wrangle
###################################################################

# Path to folder containing propka outputs
# @param sim_id - simulation identifier
# @return path to relax directory
get_propka_dir <- function(sim_id = DEF_SIM_ID) {
  file.path(get_root_dir(),
            "ompg",
            "output",
            "multi_loop",
            sim_id,
            "propka")
}

# @return path for storing pKa wrangled datasets
get_pKa_dir <- function(base_output_dir = DEF_BASE_OUTPUT_DIR) {
  file.path(base_output_dir, "pKa")
}

# @return key for wrangled pKa data - e.g: can serve as lookup for list
#   of pKa datasets or as part of the output file name
get_pKa_key <-
  function(mut = "wt",
           pH = 7,
           template = "merge",
           k = DEF_MAX_ENERGY_RANK) {
    pKa_key = c()
    if (k > 0) {
      pKa_key = paste0(mut, ".pH", pH, ".", template, ".k", k, ".pka")
    } else {
      pKa_key = paste0(mut, ".pH", pH, ".", template, ".pka")
    }
    return(pKa_key)
  }

# @return path to wrangled pKa data
get_pKa_path <-
  function(base_output_dir = DEF_BASE_OUTPUT_DIR,
           mut = "wt",
           pH = 7,
           template = "merge") {
    fid = paste0(get_pKa_key(
      mut = mut,
      pH = pH,
      template = template,
      k = -1
    ),
    ".rdata")
    return(file.path(get_pKa_dir(base_output_dir = base_output_dir), fid))
  }

# Extracts simulation identifier from a file name
# e.g. "<dir>/wt.2iww.rot1.ofa0.85.ofna0.95..." will be mapped to "wt_2iww"
# @param fid - file name to extract simulation identifier from
# @return simulation identifier
to_sim_id <- function(fid) {
  # Strip directory information
  b = basename(fid)
  # Tokenize names according to "."
  l = strsplit(x = b, split = ".", fixed = TRUE)
  # Extract first token from each list element: mutant
  # http://stackoverflow.com/questions/22430365/extract-second-element-from-every-item-in-a-list
  mutant = sapply(X = l, "[", 1)
  # Extract second token from each list element: barrel
  barrel = sapply(X = l, "[", 2)
  # Concatenate to generate simulation identifier <mutant>_<barrel>
  sim_id = paste0(mutant, "_", barrel)
  return(sim_id)
}

# @return determines the propka file name used for sample energy identifier
ener_id_to_propka <- function(ener_ids) {
  library(tools)
  fid = file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(ener_ids))))
  return(paste0(fid, ".pka"))
}

# @param rdata_path - output rdata path
# @param k - energy rank cuttoff - only k lowest energy samples retained
# @param ener_ids - energy identifiers to load pKa data for
# @param ener_id_to_propka_func - function callback for converting an energy
#   identifer to the corresponding propka filename (sans directory)
# @param overwrite - If TRUE, will overwrite any cached rdata
# @return pKa matrix with column names <3LETTERAA>.<RESNO> and
#   row names equal to sample name which is same as names(oc) and
#   ener$NAMES. The value of each row is the propka pKa at each
#   of the ionizable residues. Filtering is performed such that only
#   NAMD modifiable residues (ASP, GLU, HIS, LYS) are retained.
load_pKa <- function(rdata_path = get_pKa_path(),
                     k = DEF_MAX_ENERGY_RANK,
                     ener_ids = SIM_WT$pH7$merge$ener$NAME,
                     ener_id_to_propka_func = ener_id_to_propka,
                     overwrite = DEF_SHOULD_OVERWRITE) {
  stopifnot(length(ener_ids) > 0)
  
  if (k > length(ener_ids)) {
    print("WARNING: k > length(ener_ids), setting k = length(ener_ids)")
    k = length(ener_ids)
  } else if (k <= 0) {
    k = length(ener_ids)
  }
  
  # @HACK - supersetting works easier with a csv path
  csv_path = get_csv_path_from_rdata(rdata_path)
  
  # Determine path to binary cache (energy croppped data set)
  rdata_path = get_cropped_rdata_path_from_csv(csv_path, k)
  
  # Skip file if rdata exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping pKa format as rdata already exists: ",
                 rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Check if superset exists
  if (!overwrite) {
    rdata_superset_path = get_superset_path_from_csv(csv_path, k)
    if (file.exists(rdata_superset_path)) {
      superset = load_rdata(rdata_superset_path)
      data = superset[1:k,]
      # Write binary, compressed data
      save_rdata(data, rdata_path)
      return(data)
    }
  }
  
  # Crop energy identifiers
  stopifnot(k > 0)
  stopifnot(k <= length(ener_ids))
  ener_ids = ener_ids[1:k]
  
  # Determine propka paths for each energy sample
  propka_dir = get_propka_dir(to_sim_id(ener_ids))
  propka_name = ener_id_to_propka_func(ener_ids)
  propka_path = file.path(propka_dir, propka_name)
  
  # Utility for processing a propka file
  # @return pKa vector with names in format <3LETTERAA>.<RESNO>
  parse_propka <- function(path, erank) {
    print(paste0(erank, ": Processing ", path))
    # Open file connection
    con <- file(path)
    open(con)
    
    lines = readLines(con = con)
    close(con)
    
    # Determine pKa output region
    # Add +2 to skip title line and column header line
    pKa_start = which(grepl(
      pattern = "SUMMARY OF THIS PREDICTION",
      x = lines,
      fixed = TRUE
    )) + 2
    # Subtract -3 to skip back to last residue reported
    pKa_stop = which(
      grepl(
        pattern = "Free energy of   folding (kcal/mol) as a function of pH (using neutral reference)",
        x = lines,
        fixed = TRUE
      )
    ) - 3
    pKa_lines = lines[pKa_start:pKa_stop]
    # Tokenize by whitespace
    pKa_tokens = strsplit(x = pKa_lines, "\\s+")
    
    # Extract relevant columns
    res_name = sapply(X = pKa_tokens, "[", 2)
    res_no = sapply(X = pKa_tokens, "[", 3)
    res_pKa = as.numeric(sapply(X = pKa_tokens, "[", 5))
    
    # Keep only ionizable residues which can be modified in NAMD
    keep = res_name %in% c("ASP", "GLU", "HIS", "LYS")
    res_name = res_name[keep]
    res_no = res_no[keep]
    res_pKa = res_pKa[keep]
    
    # Add names
    names(res_pKa) = paste0(res_name, ".", res_no)
    return(res_pKa)
  }
  
  # Process first element to determine output dimensions
  res_pKa = parse_propka(propka_path[1], 1)
  pKa_mat = matrix(
    data = NA,
    nrow = length(propka_path),
    ncol = length(res_pKa)
  )
  rownames(pKa_mat) = ener_ids
  colnames(pKa_mat) = names(res_pKa)
  pKa_mat[1,] = res_pKa
  na_vec = rep(NA, nrow(pKa_mat))
  
  # Parse remaining propka output
  if (length(propka_path) > 1) {
    for (erank in 2:length(propka_path)) {
      path = propka_path[erank]
      res_pKa = parse_propka(path = path, erank = erank)
      
      #
      # NOTE: It is possible for propka to not output a pKa value
      # (not quite sure of reason). Therefore, it is possible that
      # the output matrix, pKa_mat, is missing columns and so loop
      # must check each propka file to see if a new column should
      # be added.
      #
      
      missing = !(names(res_pKa) %in% colnames(pKa_mat))
      if (any(missing)) {
        fin_col_names = colnames(pKa_mat)
        fin_row_names = rownames(pKa_mat)
        missing_names = names(res_pKa)[missing]
        # YUCK!
        for (missing_name in missing_names) {
          pKa_mat = cbind(pKa_mat, na_vec)
          fin_col_names = c(fin_col_names, missing_name)
        }
        colnames(pKa_mat) = fin_col_names
        # Just to be safe
        rownames(pKa_mat) = fin_row_names
      }
      stopifnot(all(names(res_pKa) %in% colnames(pKa_mat)))
      
      # Store pKa values for this sample (missing values will be NA)
      pKa_mat[erank, names(res_pKa)] = res_pKa
    }
  }
  
  # Save rdata
  save_rdata(data = pKa_mat,
             rdata_path = rdata_path,
             b_create_subdirs = TRUE)
  return(pKa_mat)
}

###################################################################
# Boot
###################################################################

##############################
# Libraries

library(boot)
library(matrixStats)
library(parallel)

##############################
# Globals

# Default number of bootstrap trials
NUM_BOOTSTRAP_TRIALS = 20000

# Default number of parallel estimates
NUM_BOOTSTRAP_CPUS = detectCores()

# Default significance level
BOOTSTRAP_ALPHA = 0.05

##############################
# Stats

# Computes difference in proportion protonated between two states
# @param obs - Matrix of observations where each row is a
#   multivariate observation: 1 if residue is protonated, 0 o/w
# @param i - Observation indices - corresponds to the rows
#   of data to use for computing statistic set
# @param n.first - The first n.first rows of data correspond to
#   number of samples in first state. Remaining rows are assumed
#   to be from second state.
# @return vector of mean differences
pKa_prot_diff_stat <- function(obs, i, n.first) {
  stopifnot(n.first > 0)
  i.first = i[i <= n.first]
  i.second = i[i > n.first]
  stopifnot(length(i.first) == n.first)
  stopifnot(length(i.second) == (length(i) - n.first))
  obs.first = obs[i.first, ]
  obs.second = obs[i.second, ]
  prot.first = colSums(obs.first) / n.first
  stopifnot((prot.first >= 0.0) && (prot.first <= 1.0))
  prot.second = colSums(obs.second) / (length(i) - n.first)
  stopifnot((prot.second >= 0.0) && (prot.first <= 1.0))
  prot.diff = prot.first - prot.second
  return(prot.diff)
}

# Determines if difference is consistently greater in first vs second state
# @param meta - Data.frame contain $lo_id, $hi_id, and $fav_id columns
# @param lo_id - Name of column containing lower confidence interval
# @param hi_id - Name of column containing higher confidence interval
# @param fav_id - Name of column to write favorability label
# @param fav_label_a - Name used in summary data if first state is favored (no whitespace)
# @param fav_label_b - Name used in summary data if second state is favored (no whitespace)
# @param fav_label_none - Name used in summary data if neither state is favored (no whitespace)
# @return meta data.frame with target favorability column filled with favorability label
#   indicating which state is favored (if any)
pKa_prot_diff_fav <-
  function(meta,
           lo_id,
           hi_id,
           fav_id,
           fav_label_a,
           fav_label_b,
           fav_label_none) {
    # Assign default: neither state is favored
    meta[, fav_id] = fav_label_none
    
    # Check if open state is more favored
    fav_first = (meta[, lo_id] > 0)
    meta[fav_first, fav_id] = fav_label_a
    
    # Check if close state is more favored
    fav_second = (meta[, hi_id] < 0)
    meta[fav_second, fav_id] = fav_label_b
    return(meta)
  }

##############################
# Data paths

# @return path to pKa bootstrap summary data
get_pKa_boot_sum_path <-
  function(base_output_dir = DEF_BASE_OUTPUT_DIR,
           sum_label_a = "wt.pH7.merge.open",
           sum_label_b = "wt.pH7.merge.close",
           stat_label = "PROT_DIFF",
           alpha = BOOTSTRAP_ALPHA,
           num_trials = NUM_BOOTSTRAP_TRIALS) {
    file.path(
      get_pKa_dir(base_output_dir),
      paste0(
        "boot.",
        sum_label_a,
        ".V.",
        sum_label_b,
        ".",
        stat_label,
        ".a",
        alpha,
        ".nt",
        num_trials,
        ".rdata"
      )
    )
  }

# @return path to bootstrap replicate data path
get_pKa_boot_rep_path <-
  function(pKa_boot_sum_path = get_pKa_boot_sum_path()) {
    return(insert_attrib_in_rdata(pKa_boot_sum_path, "reps"))
  }

##############################
# Core

# @param pKa_mat_a - first matrix containing pKa data
#   (columns are amino acids, rows are samples)
# @param pKa_mat_b - second matrix containing pKa data in same format
#   as pKa_mat_a
# @param pH_a - the environmental pH of first state, an amino acid
#   is protonated if pKa < pH
# @param pH_b - the environmetnal pH of second state
# @param sum_label_a - Name of first data set (pKa_mat_a under pH).
#   Used for determining output path to bootstrap summary data. Not
#   needed if sum_rdata_path is provided
# @param sum_label_b - Name of second data set (pKa_mat_b under pH).
#   Used for determining output path to bootstrap summary data. Not
#   needed if sum_rdata_path is provided
# @param stat - Bootstrap statistic to compute confidence interval for
# @param stat_label - Name of statistic
# @param fav - Method for determining favorability towards first state v second state
# @param fav_label_a - Name used in summary data if first state is favored (no whitespace)
# @param fav_label_b - Name used in summary data if second state is favored (no whitespace)
#   (must be different from fav_label_a)
# @param fav_label_none - Name used in summary data if neither state is favored (no whitespace)
#   (must be different for fav_label_a and fav_label_b)
# @param rep_rdata_path - Path to bootstrap replicate data
# @param overwrite - TRUE if bootstrap data should be recomputed
# @return summary data.frame where each row is an amino acid with rowname: <3LETTERAA.RESNO>
#   and with columns:
#     $RES.NAME - residue name
#     $RES.NO - residue number
#     $STAT - statistic computed
#     $OBS - observed value of statistic
#     $CONF.CI - Span of confidence interval = 1 - alpha
#     $CONF.LO - Lower bound of confidence interval
#     $CONF.HI - Upper bound of confidence interval
#     $CONF.FAV - Which state is favored (more protonated)
#     $P.VALUE - heuristic p-value of confidence interval
#     $FDR - false discovery rate correction of each p-value
#     $BON - Bonferroni correction
#     $FCR.CI - Span of false confidence region (FCR)
#     $FCR.LO - Lower bound of FCR
#     $FCR.HI - Upper bound of FCR
#     $FCR.FAV - Which state is favored according to FCR (more protonated)
get_pKa_boot <-
  function(pKa_mat_a = load_pKa()[SIM_WT$pH7$merge$oc,],
           pKa_mat_b = load_pKa()[!SIM_WT$pH7$merge$oc,],
           pH_a = 7,
           pH_b = 7,
           sum_label_a = "wt.pH7.merge.open",
           sum_label_b = "wt.pH7.merge.close",
           stat = pKa_prot_diff_stat,
           stat_label = "PROT_DIFF",
           fav = pKa_prot_diff_fav,
           fav_label_a = "OPEN",
           fav_label_b = "CLOSE",
           fav_label_none = "NONE",
           alpha = BOOTSTRAP_ALPHA,
           num_trials = NUM_BOOTSTRAP_TRIALS,
           num_cpus = NUM_BOOTSTRAP_CPUS,
           sum_rdata_path = get_pKa_boot_sum_path(
             sum_label_a = sum_label_a,
             sum_label_b = sum_label_b,
             stat_label = stat_label,
             alpha = alpha,
             num_trials = num_trials
           ),
           rep_rdata_path = get_pKa_boot_rep_path(sum_rdata_path),
           overwrite = DEF_SHOULD_OVERWRITE)
  {
    stopifnot(alpha > 0.0 && alpha < 1.0)
    stopifnot(fav_label_a != fav_label_b)
    stopifnot(fav_label_none != fav_label_a)
    stopifnot(fav_label_none != fav_label_b)
    
    #########################################
    # Check if cached model exists
    
    if (!overwrite && file.exists(sum_rdata_path)) {
      print(paste0(
        "Skipping bootstrap as rdata already exists: ",
        sum_rdata_path
      ))
      return(load_rdata(sum_rdata_path))
    }
    
    #########################################
    # Only work with common columns
    
    common_cols = intersect(colnames(pKa_mat_a), colnames(pKa_mat_b))
    if ((length(common_cols) != ncol(pKa_mat_a)) ||
        (length(common_cols) != ncol(pKa_mat_b))) {
      print("Warning, removing columns not common to both data sets: ")
      all_cols = union(colnames(pKa_mat_a), colnames(pKa_mat_b))
      diff_cols = setdiff(all_cols, common_cols)
      print(diff_cols)
    }
    stopifnot(length(common_cols) > 0)
    
    # Make data set columns be in same order
    pKa_mat_a = pKa_mat_a[, common_cols]
    pKa_mat_b = pKa_mat_b[, common_cols]
    
    #########################################
    # Impute missing data: replace with median
    
    do_impute_colMedian <- function(mat) {
      logi.na = is.na(mat)
      if (any(logi.na)) {
        # Compute column medians
        colMeds = colMedians(x = mat, na.rm = TRUE)
        # Obtain index matrix of missing values
        ix = which(logi.na, arr.ind = TRUE)
        # Replace missing values with corresponding column median
        mat[ix] = colMeds[ix[, 2]]
      }
      return(mat)
    }
    
    # Am imputing for perf reasons:
    #   Ideally, we would just remove samples with missing data
    #   but this complicates the bootstrap statistic as we
    #   then can't do all columns in parallel. Essentially, propka
    #   failed to compute a pKa for a very small fraction of
    #   residues in a very small fraction of samples, so let's
    #   just impute these values by replacing with ensemble
    #   median.
    pKa_mat_a_impute = do_impute_colMedian(pKa_mat_a)
    pKa_mat_b_impute = do_impute_colMedian(pKa_mat_b)
    stopifnot(!any(is.na(pKa_mat_a_impute)))
    stopifnot(!any(is.na(pKa_mat_b_impute)))
    
    #########################################
    # Determine protonation states for each set
    
    do_prot <- function(pKa_mat, pH, fav_label) {
      mat_prot = matrix(data = 0.0, nrow = nrow(pKa_mat), ncol(pKa_mat))
      # Avoid possible name clashes when munging data by adding fav label
      rownames(mat_prot) = paste0(rownames(pKa_mat), ".", fav_label)
      colnames(mat_prot) = colnames(pKa_mat)
      is_prot = pKa_mat > pH
      mat_prot[is_prot] = 1.0
      return(mat_prot)
    }
    
    prot_mat_a = do_prot(pKa_mat_a_impute, pH_a, fav_label_a)
    prot_mat_b = do_prot(pKa_mat_b_impute, pH_b, fav_label_b)
    
    #########################################
    # Munge data sets for bootstrap
    
    n.first = nrow(prot_mat_a)
    n.second = nrow(prot_mat_b)
    obs = rbind(prot_mat_a, prot_mat_b)
    
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
                     times = c(n.first,
                               n.second)),
        R = num_trials,
        parallel = "multicore",
        ncpus = num_cpus,
        n.first = n.first
      )
    }
    
    # Check everything is expected size
    n_elem = ncol(obs)
    stopifnot(n_elem == length(rep_data$t0))
    stopifnot(n_elem == ncol(rep_data$t))
    
    #########################################
    # Initialize summary (meta) data.frame
    
    #     $RES.NAME - residue name
    #     $RES.NO - residue number
    #     $STAT - statistic computed
    #     $OBS - observed value of statistic
    
    # Tokenize residue identifiers
    res_tok = strsplit(x = colnames(obs),
                       split = ".",
                       fixed = TRUE)
    
    # Extract relevant columns
    res_name = sapply(X = res_tok, "[", 1)
    res_no = sapply(X = res_tok, "[", 2)
    
    meta = data.frame(
      RES.NAME = res_name,
      RES.NO = res_no,
      STAT = rep(stat_label, n_elem),
      OBS = rep_data$t0
    )
    
    
    #########################################
    # Compute 1-alpha confidence interval
    
    conf = 1.0 - alpha
    
    # Utility for appending confidence intervals
    do_ci <- function(meta, conf, conf_id, lo_id, hi_id, fav_id) {
      meta[, conf_id] = conf
      meta[, lo_id] = 0.0
      meta[, hi_id] = 0.0
      meta[, fav_id] = "NONE"
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
      
      # Determine which state is favored
      meta = fav(meta,
                 lo_id,
                 hi_id,
                 fav_id,
                 fav_label_a,
                 fav_label_b,
                 fav_label_none)
      
      return(meta)
    }
    
    meta = do_ci(meta, conf, "CONF.CI", "LO.CI", "HI.CI", "FAV.CI")
    
    #########################################
    # Compute p-value
    
    meta$P.VALUE = 0.0
    for (i in 1:n_elem) {
      # Munge observed stat and boot stats into single set
      reps = c(rep_data$t0[i], rep_data$t[, i])
      # Calculate p-value as min(% <= 0, % >= 0)
      n = length(reps)
      meta$P.VALUE[i] = min(sum(reps <= 0.0) / n, sum(reps >= 0.0) / n)
    }
    
    #########################################
    # Compute FDR
    
    meta$FDR = p.adjust(p = meta$P.VALUE, method = "fdr")
    meta$BON = p.adjust(p = meta$P.VALUE, method = "bonferroni")
    
    #########################################
    # Compute FCR
    
    # Determine number of "significant" elements according to alpha
    k = sum(meta$FDR < alpha)
    # Adjust confidence region
    conf = 1.0 - (k * alpha / n_elem)
    
    meta = do_ci(meta, conf, "CONF.FCR", "LO.FCR", "HI.FCR", "FAV.FCR")
    
    #########################################
    # Save data
    
    # Summary
    attr(x = meta, "boot.id") = get_fid_from_rdata(sum_rdata_path)
    save_rdata(meta, sum_rdata_path)
    write.csv(x = meta,
              quote = FALSE,
              file = get_csv_path_from_rdata(sum_rdata_path))
    
    # Reps
    save_rdata(rep_data, rep_rdata_path)
    
    
    return(meta)
  }

###################################################################
# Plot
###################################################################

library(ggplot2)

source(file.path(get_script_dir(), "theme_Publication.R"))

# @return Path to FCR bar plot png
get_pKa_fcr_bar_plot_path <-
  function(base_output_dir = DEF_BASE_OUTPUT_DIR,
           fid = attr(get_pKa_boot(), "boot.id"),
           ext = ".png") {
    file.path(get_pKa_dir(),
              paste0(fid, ext))
  }


# http://stackoverflow.com/questions/4559229/drawing-pyramid-plot-using-r-and-ggplot2
# https://learnr.wordpress.com/2009/09/24/ggplot2-back-to-back-bar-charts/
# http://docs.ggplot2.org/0.9.3.1/geom_errorbar.html
# Generates bar plots of e.g. the median net stability
# @param bt.ci - Data.frame with boot strap confidence intervals
# @param main - Title of plot
# @param submain = Sub-title of plot
# @param xlab - X-axis label
# @param ylab - Y-axis label
# @param topk - Only top-k features are plotted
# @param fill_lo - The low color value for the bar plot fill
# @param fill_hi - The high color vlaue for the bar plot fill
# @param col_err - The color for error bar line
# @param top_k - keeps only top-k ranked |TO| (observed statistic) samples, 0 is no filtering
# @return ggplot object if successful, NULL o/w
# See https://coolors.co/
# http://www.colorhexa.com/
plot_pKa_fcr_bar <- function(bt.ci = get_pKa_boot(),
                             main = paste0("Bootstrap Top Net Protonation, ",
                                           round(100.0 * bt.ci$CONF.FCR[1], 1),
                                           " % FCR"),
                             submain = attr(bt.ci, "boot.id"),
                             xlab = "OmpG Ionizable Residue",
                             ylab = "Net Protonation",
                             fill_lo = "cadetblue",
                             fill_hi = "greenyellow",
                             col_err = "darksalmon",
                             top_k = -1) {
  # @TODO - make optional
  # Keep only significant entries
  bt.ci = bt.ci[(bt.ci$FAV.FCR != "NONE"), ]
  if (nrow(bt.ci) < 1) {
    return(NULL)
  }
  
  # Keep only top-ranked "interesting" samples
  if ((top_k > 0) && (nrow(bt.ci) > top_k)) {
    abs.t0 = abs(bt.ci$OBS)
    ord = order(abs.t0, decreasing = TRUE)
    bt.ci = bt.ci[ord[1:top_k], ]
  }
  
  # Order by decreasing observed statistic
  ix_order = order(bt.ci$OBS, decreasing = TRUE)
  bt.ci = bt.ci[ix_order,]
  
  bt.ci$feature = rownames(bt.ci)
  
  limits <-
    aes(ymax = HI.FCR, ymin = LO.FCR)
  
  dodge <- position_dodge(width = 0.9)
  
  p = if (nzchar(submain)) {
    ggplot(bt.ci,
           aes(
             fill = OBS,
             y = OBS,
             x = reorder(feature, OBS, function(x) {
               -x
             })
           )) +
      geom_bar(stat = "identity",
               position = dodge) +
      scale_fill_gradient(low = fill_lo, high = fill_hi) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(bquote(atop(.(main), atop(italic(
        .(submain)
      ), "")))) +
      theme_Publication(base_size = 12) +
      theme(axis.text.x = element_text(angle = -90),
            legend.position = "none") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      geom_errorbar(limits,
                    position = dodge,
                    width = 0.25,
                    color = col_err)
  } else {
    ggplot(bt.ci,
           aes(
             fill = OBS,
             y = OBS,
             x = reorder(feature, OBS, function(x) {
               -x
             })
           )) +
      geom_bar(stat = "identity",
               position = dodge) +
      scale_fill_gradient(low = fill_lo, high = fill_hi) +
      labs(x = xlab,
           y = ylab) +
      ggtitle(main) +
      theme_Publication(base_size = 12) +
      theme(axis.text.x = element_text(angle = -90),
            legend.position = "none") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      geom_errorbar(limits,
                    position = dodge,
                    width = 0.25,
                    color = col_err)
  }
  return(p)
}

###################################################################
# Shell
###################################################################

#####################
# Wrangle

do_pKa_wrangle_for_sim <- function(sim = SIM_WT,
                                   sim_name = "wt") {
  pH.ls = c(5, 7)
  template.ls = c("merge", "b2iww", "b2iwv")
  # Sort from high to low for best performance:
  k.ls = c(5000)
  
  out.ls = list()
  
  for (pH in pH.ls) {
    pH.id = paste0("pH", pH)
    for (template in template.ls) {
      ener_ids = sim[[pH.id]][[template]][["ener"]][, "NAME"]
      for (k in k.ls) {
        pKa_path = get_pKa_path(mut = sim_name,
                                pH = pH,
                                template = template)
        pKa_key = get_pKa_key(
          mut = sim_name,
          pH = pH,
          template = template,
          k = k
        )
        pKa_data =  load_pKa(rdata_path = pKa_path,
                             k = k,
                             ener_ids = ener_ids)
        out.ls[[pKa_key]] = pKa_data
      } # end iteration over energy cutoffs
    } # end iteration over sub-simulation type
  } # end iteration over pH level
  
  return(out.ls)
}

#####################
# Boot

do_pKa_boot_for_sim <- function(sim = SIM_WT,
                                sim_name = "wt") {
  pKa_raw = do_pKa_wrangle_for_sim(sim = sim, sim_name = sim_name)
  
  out.ls = list()
  
  # Experiment: pH7 merge vs pH 5 merge
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.merge.k5000.pka")]],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.merge.k5000.pka")]],
    pH_a = 7,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH7.merge.k5000"),
    sum_label_b = paste0(sim_name, ".pH5.merge.k5000"),
    fav_label_a = "PH7",
    fav_label_b = "PH5"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 7 b2iwv vs pH 5 b2iww
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.b2iwv.k5000.pka")]],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.b2iww.k5000.pka")]],
    pH_a = 7,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH7.b2iwv.k5000"),
    sum_label_b = paste0(sim_name, ".pH5.b2iww.k5000"),
    fav_label_a = "PH7.2iwv",
    fav_label_b = "PH5.2iww"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 7 merge open vs pH 7 merge close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.merge.k5000.pka")]][sim$pH7$merge$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH7.merge.k5000.pka")]][!sim$pH7$merge$oc, ],
    pH_a = 7,
    pH_b = 7,
    sum_label_a = paste0(sim_name, ".pH7.merge.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH7.merge.k5000.close"),
    fav_label_a = "OPEN",
    fav_label_b = "CLOSE"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 5 merge open vs pH 5 merge close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH5.merge.k5000.pka")]][sim$pH5$merge$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.merge.k5000.pka")]][!sim$pH5$merge$oc, ],
    pH_a = 5,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH5.merge.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH5.merge.k5000.close"),
    fav_label_a = "OPEN",
    fav_label_b = "CLOSE"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 7 b2iwv open vs pH 7 b2iwv close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.b2iwv.k5000.pka")]][sim$pH7$b2iwv$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH7.b2iwv.k5000.pka")]][!sim$pH7$b2iwv$oc, ],
    pH_a = 7,
    pH_b = 7,
    sum_label_a = paste0(sim_name, ".pH7.b2iwv.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH7.b2iwv.k5000.close"),
    fav_label_a = "OPEN",
    fav_label_b = "CLOSE"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 5 b2iww open vs pH 5 b2iww close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH5.b2iww.k5000.pka")]][sim$pH5$b2iww$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.b2iww.k5000.pka")]][!sim$pH5$b2iww$oc, ],
    pH_a = 5,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH5.b2iww.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH5.b2iww.k5000.close"),
    fav_label_a = "OPEN",
    fav_label_b = "CLOSE"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 7 merge open vs pH 5 merge close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.merge.k5000.pka")]][sim$pH7$merge$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.merge.k5000.pka")]][!sim$pH5$merge$oc, ],
    pH_a = 7,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH7.merge.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH5.merge.k5000.close"),
    fav_label_a = "OPEN.PH7",
    fav_label_b = "CLOSE.PH5"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  # Experiment : pH 7 b2iwv open vs pH 5 b2iww close
  
  boot.res = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.b2iwv.k5000.pka")]][sim$pH7$b2iwv$oc, ],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.b2iww.k5000.pka")]][!sim$pH5$b2iww$oc, ],
    pH_a = 7,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH7.b2iwv.k5000.open"),
    sum_label_b = paste0(sim_name, ".pH5.b2iww.k5000.close"),
    fav_label_a = "OPEN.PH7.2iwv",
    fav_label_b = "CLOSE.PH5.2iww"
  )
  out.ls[[attr(boot.res, "boot.id")]] = boot.res
  
  return(out.ls)
}

#####################
# Plots

do_pKa_plots <- function(sim = SIM_WT,
                         sim_name = "wt") {
  boot_ls = do_pKa_boot_for_sim(sim = sim, sim_name = sim_name)
  
  # boot identifier is too large, so crop it
  hack_get_submain_from_boot_id <- function(boot.id, stat_type="PROT_DIFF") {
    # Strip "boot." prefix
    submain = strsplit(x=boot.id, split="boot.", fixed=TRUE)[[1]][2]
    # Keep everything before state_type
    submain2 = strsplit(x=submain, split=paste0(".", stat_type))[[1]][1]
    return(submain2)
  }
  
  for (bt.ci in boot_ls) {
    fid = attr(bt.ci, "boot.id")
    png_path = get_pKa_fcr_bar_plot_path(fid = fid, ext = ".png")
    svg_path = get_pKa_fcr_bar_plot_path(fid = fid, ext = ".svg")
    emf_path = get_pKa_fcr_bar_plot_path(fid = fid, ext = ".emf")
    
    p = plot_pKa_fcr_bar(bt.ci = bt.ci,
                         submain = hack_get_submain_from_boot_id(boot.id=fid))
    
    # Save plot if valid
    if (!is.null(p)) {
      ggplot2::ggsave(filename = png_path, plot = p)
      ggplot2::ggsave(filename = svg_path, plot = p)
      if (.Platform$OS.type == "windows") {
        # Extended meta format: Insert into powerpoint then "ungroup" image
        # to allow editing of vector graphic elements!
        ggplot2::ggsave(filename = emf_path, plot = p)
      } # end check for windows emf support
    } # end check if plot is valid
  } # end iteration over bootstrap results
}


do_pKa_plot_JACS <- function() {
  sim = SIM_WT
  sim_name = "wt"
  pKa_raw = do_pKa_wrangle_for_sim(sim = sim, sim_name = sim_name)
  
  bt.ci = get_pKa_boot(
    pKa_mat_a = pKa_raw[[paste0(sim_name, ".pH7.merge.k5000.pka")]],
    pKa_mat_b = pKa_raw[[paste0(sim_name, ".pH5.merge.k5000.pka")]],
    pH_a = 7,
    pH_b = 5,
    sum_label_a = paste0(sim_name, ".pH7.merge.k5000"),
    sum_label_b = paste0(sim_name, ".pH5.merge.k5000"),
    fav_label_a = "PH7",
    fav_label_b = "PH5"
  )
  
  # Note, set top_k = -1 to plot all significant residues
  top_k = 20
  p = plot_pKa_fcr_bar(
    bt.ci = bt.ci,
    main = "Net Protonation, pH 7 vs pH 5",
    #submain = paste0(round(100.0 * bt.ci$CONF.FCR[1], 1), " % FCR"),
    submain = "",
    xlab = "OmpG Residue",
    ylab = "Net Protonation",
    fill_lo = "grey",
    fill_hi = "grey",
    col_err = "black",
    top_k = top_k
  )
  
  if (!is.null(p)) {
    DOCS_DIR = file.path(get_script_dir(), "..", "..", "docs")
    JACS_DIR = file.path(DOCS_DIR, "raw", "OmpG.JACS")
    OUT_PATH = file.path(JACS_DIR, paste0("prot.pH7.v.pH5.top",top_k,".pdf"))
    dir.create(JACS_DIR, showWarnings = FALSE)
    print(paste0("Saving ", OUT_PATH))
    cowplot::ggsave(filename = OUT_PATH, plot = p)
  }
}

#####################
# Protonation fraction

# Simply reports protonation fraction for wildtype merged population at pH level
get_prot_prop_merge <- function(pH = 5,
                                res_id = "GLU.31") {
  rdata_path = get_pKa_path(pH = pH,
                            mut = "wt",
                            template = "merge")
  
  ener_ids = NULL
  if (pH == 7) {
    ener_ids = SIM_WT$pH7$merge$ener$NAME
  } else if (pH == 5) {
    ener_ids = SIM_WT$pH5$merge$ener$NAME
  } else {
    stopifnot(FALSE, "Unrecognized pH level")
  }
  
  pKa_mat = load_pKa(rdata_path = rdata_path,
                     ener_ids = ener_ids)
  pKa_mat = pKa_mat[, res_id, drop = FALSE]
  
  mat_prot = matrix(data = 0.0, nrow = nrow(pKa_mat), ncol(pKa_mat))
  rownames(mat_prot) = rownames(pKa_mat)
  colnames(mat_prot) = colnames(pKa_mat)
  is_prot = pKa_mat > pH
  mat_prot[is_prot] = 1.0
  
  prot_prop = colSums(mat_prot)
  prot_prop = prot_prop / nrow(mat_prot)
  stopifnot(prot_prop >= 0.0)
  stopifnot(prot_prop <= 1.0)
  
  print(paste0(
    "Observed protonation fraction for merged wildtype at pH=",
    pH,
    ":"
  ))
  print(prot_prop)
  #return(prot_prop)
}

# Report negative patch fractions for E17 and E31
do_get_prot_prop_merge_neg_patch <- function() {
  # Commented results are as of 10/23/2017
  # Note, could also combine res_ids into vector argument
  get_prot_prop_merge(pH = 7, "GLU.31") # 0.9316
  get_prot_prop_merge(pH = 5, "GLU.31") # 0.9988
  get_prot_prop_merge(pH = 7, "GLU.17") # 0.0704
  get_prot_prop_merge(pH = 5, "GLU.17") # 0.8882
}
