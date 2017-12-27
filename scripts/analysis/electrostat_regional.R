###################################################################
# electrostat_regional.R
#
# Analysis relating to macro-level regional stability
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
source(file.path(get_script_dir(), "electrostat_globals2.R"))

# Global list of loops. NONE is for all non-simulated regions
# which inlcudes loop 4 and the barrel.
LOOP_IDS = c(paste0("LOOP_", c(1, 2, 3, 5, 6, 7)), "NONE")

###################################################################
# Regional interaction filters
###################################################################

# Default regional interaction filter
DEF_REG_INT_FILT_ID = "none"

# Applies *no* filtering to design and meta parameters
regional_int_filt_none <- function(meta) {
  return(rep(TRUE, nrow(meta)))
}

# Selects interactions only between ionizable residues
regional_int_filt_ion <- function(meta) {
  return(meta$ION.ION)
}

# @return regional interaction filter from string identifier
get_reg_int_filt <- function(filt_int_id = DEF_REG_INT_FILT_ID) {
  return(list("none" = regional_int_filt_none,
              "ion" = regional_int_filt_ion)[[filt_int_id]])
}

###################################################################
# Regional region filters
###################################################################

# Filter by loop regions
regional_reg_filt_loop_pair <-
  function(meta, loop_id.a, loop_id.b) {
    a.b = ((meta$RESA.LOOP == loop_id.a) &
             (meta$RESB.LOOP == loop_id.b))
    b.a = ((meta$RESA.LOOP == loop_id.b) &
             (meta$RESB.LOOP == loop_id.a))
    sel = a.b | b.a
    return(sel)
  }

# Filter by residue name
regional_reg_filt_single_res <-
  function(meta, res.no) {
    sel = ((meta$RESA.NO == res.no) |
             (meta$RESB.NO == res.no))
    return(sel)
  }

# Filter by patch
regional_reg_filt_patch <- function(meta, patch_res_nos, loop_ids) {
  sel.loop = ((meta$RESA.LOOP %in% loop_ids) |
                (meta$RESB.LOOP %in% loop_ids))
  sel.patch = ((meta$RESA.NO %in% patch_res_nos) |
                 (meta$RESB.NO %in% patch_res_nos))
  sel = sel.loop & sel.patch
}

# Default regional region filter
DEF_REG_REG_FILT = regional_reg_filt_loop_pair
# Default regional region filter args list
DEF_REG_REG_FILT_ARGS_LS = list(loop_id.a = "LOOP_6", loop_id.b = "LOOP_6")

###################################################################
# Specific path utilities
###################################################################

# Path for storing intermediate regional scores
get_regional_scores_dir <-
  function(base_output_dir = get_def_base_output_dir(),
           pH = 5,
           sim_id = "wt") {
    return(file.path(
      get_capt_dir(
        base_output_dir = base_output_dir,
        pH = pH,
        sim_id = sim_id
      ),
      "regsc"
    ))
  }

# @return path to regional scores rdata
get_regional_scores_rdata_path <-
  function(regional_scores_dir = get_regional_scores_dir(),
           filt_int_id = DEF_REG_INT_FILT_ID,
           filt_reg_id) {
    file.path(regional_scores_dir,
              paste0("regsc.", filt_int_id, ".", filt_reg_id, ".rdata"))
  }

###################################################################
# Regional scores util
###################################################################

# Obtain regional scores based on interaction and region filters
# @param design - input matrix of residue-residue interaction energies
#   where columns are residue-residue interaction and rows are
#   the observed samples
# @param meta - input meta data such as loop identifiers for each
#   residue-residue interaction
# @return numeric vector with row sums of the filtered design matrix
#   - has length = nrow(design)
#   - has names = rownames(design)
# May return NULL if no elements pass selection filter
get_regional_scores_design <-
  function(design = WT_PH5_ELEC_RCB$design,
           meta = WT_PH5_ELEC_RCB$meta,
           filt_int = get_reg_int_filt(),
           filt_reg = DEF_REG_REG_FILT,
           filt_reg_args_ls = DEF_REG_REG_FILT_ARGS_LS) {
    stopifnot(ncol(design) == nrow(meta))
    stopifnot(all(colnames(design) == rownames(meta)))
    
    # Interaction filter
    sel_int = filt_int(meta = meta)
    # Region filter
    filt_reg_args_ls$meta = meta
    sel_reg = do.call(what = filt_reg, args = filt_reg_args_ls)
    # Combine filters
    sel = sel_int & sel_reg
    
    # Early out if no elements meet selection criteria
    if (!any(sel)) {
      return(NULL)
    }
    
    # Apply filter
    design = design[, sel]
    
    # Replace NAs with zero
    design[is.na(design)] = 0.0
    
    # Obtain scores
    reg_scores = c()
    if (sum(sel) > 1) {
      reg_scores = rowSums(x = design, na.rm = FALSE)
    } else {
      reg_scores = design
    }
    stopifnot(names(reg_scores) == rownames(design))
    
    return(reg_scores)
  }

###################################################################
# Munge region scores
###################################################################

# Reside-residue pairs
# @WARNING - all NAs replaced with 0.0
get_regional_scores_rp <- function(design = WT_PH5_ELEC_RCB$design,
                                   meta = WT_PH5_ELEC_RCB$meta,
                                   filt_int = get_reg_int_filt(),
                                   rdata_path = get_regional_scores_rdata_path(filt_reg_id = "rp"),
                                   overwrite = SHOULD_OVERWRITE) {
  # Check for cached data
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Apply interaction filter to design matrix and meta frame
  sel_int = filt_int(meta = meta)
  stopifnot(any(sel_int))
  out.design = design[, sel_int]
  # Replace NAs with zero
  out.design[is.na(out.design)] = 0.0
  out.meta = meta[sel_int,]
  
  # Save data set
  out_ls = list(meta = out.meta,
                design = out.design)
  save_rdata(out_ls, rdata_path)
  return(out_ls)
}

# Loop-loop pairs
get_regional_scores_lp <- function(design = WT_PH5_ELEC_RCB$design,
                                   meta = WT_PH5_ELEC_RCB$meta,
                                   filt_int = get_reg_int_filt(),
                                   rdata_path = get_regional_scores_rdata_path(filt_reg_id = "lp"),
                                   overwrite = SHOULD_OVERWRITE) {
  # Check for cached data
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Determine iteration elements
  loop_ids = LOOP_IDS
  n_loops = length(loop_ids)
  n_elems = ((n_loops * n_loops - n_loops) / 2) + n_loops - 1
  
  # Allocate output meta
  none_vec = rep("NONE", n_elems)
  out.meta = data.frame(LOOP.A = none_vec,
                        LOOP.B = none_vec,
                        stringsAsFactors = FALSE)
  
  # Allocate output design
  out.design = matrix(data = 0.0,
                      nrow = nrow(design),
                      ncol = n_elems)
  rownames(out.design) = rownames(design)
  colnames(out.design) = 1:n_elems
  
  # Iterate over regions
  filt_reg_args_ls = list(loop_id.a = "NONE",
                          loop_id.b = "NONE")
  k = 1
  for (i in 1:(n_loops - 1)) {
    loop_id.a = loop_ids[i]
    filt_reg_args_ls$loop_id.a = loop_id.a
    for (j in i:n_loops) {
      loop_id.b = loop_ids[j]
      filt_reg_args_ls$loop_id.b = loop_id.b
      reg_name = paste0(loop_id.a, ".", loop_id.b)
      rownames(out.meta)[k] = reg_name
      colnames(out.design)[k] = reg_name
      out.meta$LOOP.A[k] = loop_id.a
      out.meta$LOOP.B[k] = loop_id.b
      print(paste0("--------Processing ", reg_name))
      
      tmp = get_regional_scores_design(
        design = design,
        meta = meta,
        filt_int = filt_int,
        filt_reg = regional_reg_filt_loop_pair,
        filt_reg_args_ls = filt_reg_args_ls
      )
      
      # Add score if selection was not empty (not NULL)
      if (!is.null(tmp)) {
        out.design[, k] = tmp
        k = k + 1
      } else {
        print(paste0("Warning: No scores available for ", reg_name))
      }
    } # end iteration over second loop
  } # end iteration over first loop
  
  # Shrink outputs to account for any NULL elements
  if (k <= n_elems) {
    out.design = out.design[, 1:(k - 1)]
    out.meta = out.meta[1:(k - 1), ]
  }
  
  # Save data set
  out_ls = list(meta = out.meta,
                design = out.design)
  save_rdata(out_ls, rdata_path)
  return(out_ls)
}

# Single residue scores
get_regional_scores_srs <- function(design = WT_PH5_ELEC_RCB$design,
                                    meta = WT_PH5_ELEC_RCB$meta,
                                    filt_int = get_reg_int_filt(),
                                    loop_keep_ids = LOOP_IDS,
                                    hack_only_ion = FALSE,
                                    rdata_path = get_regional_scores_rdata_path(filt_reg_id = "srs"),
                                    overwrite = SHOULD_OVERWRITE) {
  # Check for cached data
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping as rdata already exists: ", rdata_path))
    return(load_rdata(rdata_path))
  }
  
  # Determine iteration elements
  res_nos = unique(c(meta$RESA.NO, meta$RESB.NO))
  n_elems = length(res_nos)
  
  # Allocate output meta
  zero_vec = as.integer(rep(0, n_elems))
  none_vec = rep("NONE", n_elems)
  false_vec = rep(FALSE, n_elems)
  out.meta = data.frame(
    RES.NO = zero_vec,
    RES.NAME = none_vec,
    RES.LOOP = none_vec,
    RES.ION = false_vec,
    RES.POL = false_vec,
    stringsAsFactors = FALSE
  )
  rownames(out.meta) = res_nos
  
  # Allocate output design
  out.design = matrix(data = 0.0,
                      nrow = nrow(design),
                      ncol = n_elems)
  rownames(out.design) = rownames(design)
  colnames(out.design) = res_nos
  
  meta_col_sfx = c("NO", "NAME", "LOOP", "ION", "POL")
  IX_COL_NO = 1
  IX_COL_NAME = 2
  IX_COL_LOOP = 3
  IX_COL_ION = 4
  IX_COL_POL = 5
  meta_col_a = paste0("RESA.", meta_col_sfx)
  meta_col_b = paste0("RESB.", meta_col_sfx)
  
  # Pre-alloc some loop vars
  meta_ix = 0
  meta_col = meta_col_a
  
  # Iterate over regions
  filt_reg_args_ls = list(res.no = "NONE")
  keep = rep(TRUE, n_elems)
  for (k in 1:n_elems) {
    res.no = res_nos[k]
    filt_reg_args_ls$res.no = res.no
    
    # Determine row, col for meta info
    sel = meta$RESA.NO == res.no
    meta_col = meta_col_a
    if (!any(sel)) {
      sel = meta$RESB.NO == res.no
      stopifnot(any(sel))
      meta_col = meta_col_b
    }
    meta_ix = which(sel)[1]
    
    # @HACK - skip if loop is not being tracked
    res_loop = meta[meta_ix, meta_col[IX_COL_LOOP]]
    if (!(res_loop %in% loop_keep_ids)) {
      keep[k] = FALSE
      next
    }
    
    # @HACK - skip if only tracking ionizable residues
    is_ion = meta[meta_ix, meta_col[IX_COL_ION]]
    if (hack_only_ion && !is_ion) {
      keep[k] = FALSE
      next
    }
    
    # Copy meta data
    out.meta$RES.NO[k] = meta[meta_ix, meta_col[IX_COL_NO]]
    out.meta$RES.NAME[k] = meta[meta_ix, meta_col[IX_COL_NAME]]
    out.meta$RES.LOOP[k] = res_loop
    out.meta$RES.ION[k] = is_ion
    out.meta$RES.POL[k] = meta[meta_ix, meta_col[IX_COL_POL]]
    
    print(paste0("--------Processing ", res.no))
    tmp = get_regional_scores_design(
      design = design,
      meta = meta,
      filt_int = filt_int,
      filt_reg = regional_reg_filt_single_res,
      filt_reg_args_ls = filt_reg_args_ls
    )
    
    # Only append non-NULL elements
    if (!is.null(tmp)) {
      out.design[, k] = tmp
    } else {
      print(paste0("Warning: No scores available for ", res.no))
      keep[k] = FALSE
    }
  }
  
  # Apply hack filters
  out.meta = out.meta[keep,]
  out.design = out.design[, keep]
  
  # Save data set
  out_ls = list(meta = out.meta,
                design = out.design)
  save_rdata(out_ls, rdata_path)
  return(out_ls)
}

# Patch-loop pairs
get_regional_scores_patch <-
  function(design = WT_PH5_ELEC_RCB$design,
           meta = WT_PH5_ELEC_RCB$meta,
           filt_int = get_reg_int_filt(),
           rdata_path = get_regional_scores_rdata_path(filt_reg_id = "patch"),
           overwrite = SHOULD_OVERWRITE) {
    # Check for cached data
    if (!overwrite && file.exists(rdata_path)) {
      print(paste0("Skipping as rdata already exists: ", rdata_path))
      return(load_rdata(rdata_path))
    }
    
    # Determine iteration elements
    # Patches:
    patch_res_nos_sets = list(c(15, 17, 31, 52), c(111, 113, 211, 216))
    patch_res_nos_names = sapply(X = patch_res_nos_sets,
                                 FUN = paste,
                                 collapse = ".")
    patch_res_nos_desc = c("NEG", "POS")
    stopifnot(length(patch_res_nos_sets) == length(patch_res_nos_names))
    stopifnot(length(patch_res_nos_sets) == length(patch_res_nos_desc))
    # Loops:
    l.ids = paste0("LOOP_", c(1, 2, 3, 5, 6, 7))
    loop_ids_sets = list(l.ids[1], l.ids[2], l.ids[3], l.ids[4], l.ids[5], l.ids[6], l.ids)
    loop_ids_names = c(l.ids, "LOOP_ALL")
    stopifnot(length(loop_ids_sets) == length(loop_ids_names))
    
    n_elems = length(patch_res_nos_sets) * length(loop_ids_sets)
    
    # Allocate output meta
    none_vec = rep("NONE", n_elems)
    out.meta = data.frame(
      PATCH_ID = none_vec,
      PATCH_DESC = none_vec,
      LOOP_ID = none_vec,
      stringsAsFactors = FALSE
    )
    
    # Allocate output design
    out.design = matrix(data = 0.0,
                        nrow = nrow(design),
                        ncol = n_elems)
    rownames(out.design) = rownames(design)
    colnames(out.design) = 1:n_elems
    
    # Iterate over regions
    filt_reg_args_ls = list(patch_res_nos = 0, loop_ids = "NONE")
    k = 1
    for (i in 1:length(patch_res_nos_sets)) {
      patch_res_nos = patch_res_nos_sets[[i]]
      filt_reg_args_ls$patch_res_nos = patch_res_nos
      p.name = patch_res_nos_names[[i]]
      p.desc = patch_res_nos_desc[[i]]
      
      for (j in 1:length(loop_ids_sets)) {
        loop_ids = loop_ids_sets[[j]]
        filt_reg_args_ls$loop_ids = loop_ids
        l.name = loop_ids_names[[j]]
        
        reg_name = paste0(p.name, "-", l.name)
        rownames(out.meta)[k] = reg_name
        colnames(out.design)[k] = reg_name
        
        out.meta$PATCH_ID[k] = p.name
        out.meta$PATCH_DESC[k] = p.desc
        out.meta$LOOP_ID[k] = l.name
        
        print(paste0("--------Processing ", reg_name))
        tmp = get_regional_scores_design(
          design = design,
          meta = meta,
          filt_int = filt_int,
          filt_reg = regional_reg_filt_patch,
          filt_reg_args_ls = filt_reg_args_ls
        )
        
        # Add score if selection was not empty (not NULL)
        if (!is.null(tmp)) {
          out.design[, k] = tmp
          k = k + 1
        } else {
          print(paste0("Warning: No scores available for ", reg_name))
        }
      } # end iteration over loop sets
    } # end iteration over patches
    
    # Shrink outputs to account for any NULL elements
    if (k <= n_elems) {
      out.design = out.design[, 1:(k - 1)]
      out.meta = out.meta[1:(k - 1), ]
    }
    
    # Save data set
    out_ls = list(meta = out.meta,
                  design = out.design)
    save_rdata(out_ls, rdata_path)
    return(out_ls)
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
NUM_BOOTSTRAP_TRIALS = 20000

# Default number of parallel estimates
NUM_BOOTSTRAP_CPUS = detectCores()

# Default significance level
BOOTSTRAP_ALPHA = 0.05

##############################
# Paths

# @return path to bootstrap summary data
get_reg_boot_sum_path <-
  function(regsc_rdata_path = get_regional_scores_rdata_path(filt_reg_id =
                                                               "lp"),
           stat_type = "med_diff",
           alpha = BOOTSTRAP_ALPHA) {
    return(insert_attrib_in_rdata(regsc_rdata_path, paste0("boot.", stat_type, ".", alpha)))
  }

# @return path to bootstrap replicate data path
get_reg_boot_rep_path <-
  function(reg_boot_sum_path = get_reg_boot_sum_path()) {
    return(insert_attrib_in_rdata(reg_boot_sum_path, "reps"))
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
reg_med_diff_stat <- function(obs, i, n.open) {
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
reg_mea_diff_stat <- function(obs, i, n.open) {
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

# Determines if difference is consistently greater in open vs closed state
# @param meta - data.frame contain $lo_id, $hi_id, and $fav_id columns
# @param lo_id - name of column containing lower confidence interval
# @param hi_id - name of column containing higher confidence interval
# @param fav_id - name of column to write favorability label
# @return meta data.frame with target favorability column
#   filled with "OPEN", "CLOSE", or "NONE" indicating which
#   state is favored
reg_diff_fav <- function(meta, lo_id, hi_id, fav_id) {
  # Check if open state is more favored
  fav_open = (meta[, hi_id] < 0)
  meta[fav_open, fav_id] = "OPEN"
  
  # Check if close state is more favored
  fav_close = (meta[, lo_id] > 0)
  meta[fav_close, fav_id] = "CLOSE"
  return(meta)
}

##############################
# Core

# Computes bootstrap statistics on regional data.frame
# @param reg_scores_mat - Aggregate residue scores matrix where
#   each column is a sample and each row is a region score
# @param meta - Data frame with meta info for each row reg_scores_mat
# @param oc - Open/close status of each sample column in reg_scores_mat
# @param stat - The statistic to compute
# @param fav - Method for determining favorability towards open v close
# @param alpha - The base significance level: (1-alpha) CI
# @param num_trials - Number of bootstrap replicates to generate
# @param num_cpus - Number of parallel CPUs to use
# @param sum_rdata_path - Path to output bootstrap summary data
# @param rep_rdata_path - Path to output bootstrap replicate data
# @param overwrite - If TRUE, any cached data stored on disk is
#   overwritten
# @return data.frame with bootstrap summary data
get_reg_boot <-
  function(reg_scores_mat = get_regional_scores_lp()$design,
           meta = get_regional_scores_lp()$meta,
           oc = WT_PH5_OC_RCB,
           stat = reg_med_diff_stat,
           fav = reg_diff_fav,
           alpha = BOOTSTRAP_ALPHA,
           num_trials = NUM_BOOTSTRAP_TRIALS,
           num_cpus = NUM_BOOTSTRAP_CPUS,
           sum_rdata_path = get_reg_boot_sum_path(alpha = alpha),
           rep_rdata_path = get_reg_boot_rep_path(sum_rdata_path),
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
    
    #########################################
    # Separate open and close samples
    
    n.open = sum(oc)
    n.close = length(oc) - n.open
    names.open = names(oc)[oc]
    names.close = names(oc)[!oc]
    names.strat = c(names.open, names.close)
    # Place open before close to allow stratified statistic computation
    obs = reg_scores_mat[names.strat, ]
    
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
      meta = fav(meta, lo_id, hi_id, fav_id)
      
      return(meta)
    }
    
    meta = do_ci(meta, conf, "Conf.Ci", "Lo.Ci", "Hi.Ci", "Fav.Ci")
    
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
    
    meta = do_ci(meta, conf, "Conf.Fcr", "Lo.Fcr", "Hi.Fcr", "Fav.Fcr")
    
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

reg_bar_feat_name_ident <-
  function(bt.ci) {
    return(rownames(bt.ci))
  }

reg_bar_feat_name_lp_desc <- function(bt.ci) {
  # Names are format in LOOP_<DIGIT in 1:6>.LOOP_<DIGIT in 1:6>
  # or LOOP_<DIGIT in 1:6>.NONE
  
  lp = rownames(bt.ci)
  
  loop2loop = !grepl("NONE", x = lp, fixed = TRUE)
  digit1 = substr(lp, start = 6, stop = 6)
  digit2 = substr(lp[loop2loop], start = 13, stop = 13)
  
  first = paste0("L", digit1)
  second.master = rep("N", length(first))
  second.loop = paste0("L", digit2)
  second.master[loop2loop] = second.loop
  
  feat_names = paste0(first, ".", second.master)
  return(feat_names)
} 

# Utility to obtain a human readable residue-residue interaction name
# Converts <resno>.<resno> to <res1LetterCode><resno>.<res1LetterCode><resno>
# Eg. "68.228" becomes "R68.R228"
reg_bar_feat_name_rp_desc <- function(bt.ci) {
  resa.code = AA3to1[bt.ci$RESA.NAME]
  resa.no = bt.ci$RESA.NO
  resb.code = AA3to1[bt.ci$RESB.NAME]
  resb.no = bt.ci$RESB.NO
  return(paste0(resa.code, resa.no, ".", resb.code, resb.no))
}

reg_bar_feat_name_srs_desc <- function(bt.ci) {
  res.no = bt.ci$RES.NO
  res.code = AA3to1[bt.ci$RES.NAME]
  return(paste0(res.code, res.no))
}

reg_bar_feat_name_patch_desc <- function(bt.ci) {
  return(paste0(bt.ci$PATCH_DESC, ".", bt.ci$LOOP_ID))
}

source(file.path(get_script_dir(), "theme_Publication.R"))

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
# @param col_err - The color for error bar line
# @param top_k - keeps only top-k ranked |TO| (observed statistic) samples, 0 is no filtering
# @return ggplot object if successful, NULL o/w
# See https://coolors.co/
# http://www.colorhexa.com/
plot_reg_fcr_bar <- function(bt.ci = do_reg_boot(
                                       pH.ls = 5,
                                       sim_id.ls = "wt",
                                       filt_int_id.ls = "none",
                                       stat_type.ls = "med_diff",
                                       regsc_type.ls = "patch"
                                     )[[1]]$boot,
                             feat_name_fnc = reg_bar_feat_name_patch_desc,
                             main = "Net Stability",
                             submain = paste0(round(100.0 * bt.ci$Conf.Fcr[1], 1)," % FCR"),
                             xlab = "OmpG Patch to Residue Interactions",
                             ylab = "Median Net Stability",
                             fill_lo = "cadetblue",
                             fill_hi = "greenyellow",
                             col_err = "darksalmon",
                             top_k = 25) {
  library(ggplot2)
  
  # @TODO - make optional
  # Keep only significant entries
  bt.ci = bt.ci[(bt.ci$Fav.Fcr != "NONE"),]
  if (nrow(bt.ci) < 1) {
    return(NULL)
  }
  
  # Keep only top-ranked "interesting" samples
  if ((top_k > 0) && (nrow(bt.ci) > top_k)) {
    abs.t0 = abs(bt.ci$T0)
    ord = order(abs.t0, decreasing = TRUE)
    bt.ci = bt.ci[ord[1:top_k],]
  }
  
  # Order by decreasing observed median
  ix_order = order(bt.ci$T0, decreasing = TRUE)
  bt.ci = bt.ci[ix_order, ]
  
  feat.names = feat_name_fnc(bt.ci)
  bt.ci$feature = feat.names # as.factor(feat.names)
  
  limits <-
    aes(ymax = Hi.Fcr, ymin = Lo.Fcr)
  
  dodge <- position_dodge(width = 0.9)
  
  #http://stackoverflow.com/questions/19957536/add-dynamic-subtitle-using-ggplot
  #http://felixfan.github.io/ggplot2-remove-grid-background-margin/
  
  p = ggplot(bt.ci,
             aes(
               fill = T0,
               y = T0,
               x = reorder(feature, T0, function(x) {
                 -x
               })
             )) +
    geom_bar(stat = "identity",
             position = dodge) +
    scale_fill_gradient(low = fill_lo, high = fill_hi) +
    labs(x = xlab,
         y = ylab) +
    ggtitle(bquote(atop(.(main), atop(italic(.(submain)), "")))) +
    theme_Publication(base_size=12) +
    theme(axis.text.x = element_text(angle = -90),
          legend.position = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    geom_errorbar(limits,
                  position = dodge,
                  width = 0.25,
                  color = col_err)
  return(p)
}

###################################################################
# Shell
###################################################################

##############################
# Munge

# Residue-residue pairs
do_regional_scores_rp <- function(pH.ls = c(5, 7),
                                  sim_id.ls = c("wt"),
                                  filt_int_id.ls = c("none", "ion")) {
  elec_ls = list("wt.5" = WT_PH5_ELEC_RCB,
                 "wt.7" = WT_PH7_ELEC_RCB)
  
  results = list()
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      regsc_scores_dir = get_regional_scores_dir(pH = pH, sim_id = sim_id)
      # Make sure directory exists
      dir.create(regsc_scores_dir,
                 showWarnings = FALSE,
                 recursive = TRUE)
      base_id = paste0(sim_id, ".", pH)
      
      design = elec_ls[[base_id]]$design
      meta = elec_ls[[base_id]]$meta
      
      for (filt_int_id in filt_int_id.ls) {
        rdata_path = get_regional_scores_rdata_path(
          regional_scores_dir = regsc_scores_dir,
          filt_int_id = filt_int_id,
          filt_reg_id = "rp"
        )
        
        final_id = paste0(base_id, ".", filt_int_id)
        
        print(paste0("Processing ", rdata_path))
        results[[final_id]] = list(
          regsc = get_regional_scores_rp(
            design = design,
            meta = meta,
            filt_int = get_reg_int_filt(filt_int_id = filt_int_id),
            rdata_path = rdata_path
          ),
          rdata_path = rdata_path
        )
        
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

# Loop-loop pairs
do_regional_scores_lp <- function(pH.ls = c(5, 7),
                                  sim_id.ls = c("wt"),
                                  filt_int_id.ls = c("none", "ion")) {
  elec_ls = list("wt.5" = WT_PH5_ELEC_RCB,
                 "wt.7" = WT_PH7_ELEC_RCB)
  
  results = list()
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      regsc_scores_dir = get_regional_scores_dir(pH = pH, sim_id = sim_id)
      # Make sure directory exists
      dir.create(regsc_scores_dir,
                 showWarnings = FALSE,
                 recursive = TRUE)
      base_id = paste0(sim_id, ".", pH)
      
      design = elec_ls[[base_id]]$design
      meta = elec_ls[[base_id]]$meta
      
      for (filt_int_id in filt_int_id.ls) {
        rdata_path = get_regional_scores_rdata_path(
          regional_scores_dir = regsc_scores_dir,
          filt_int_id = filt_int_id,
          filt_reg_id = "lp"
        )
        
        final_id = paste0(base_id, ".", filt_int_id)
        
        print(paste0("Processing ", rdata_path))
        results[[final_id]] = list(
          regsc = get_regional_scores_lp(
            design = design,
            meta = meta,
            filt_int = get_reg_int_filt(filt_int_id = filt_int_id),
            rdata_path = rdata_path
          ),
          rdata_path = rdata_path
        )
        
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

# Single residue scores
DEFAULT_REG_SRS_LOOP_KEEP_IDS = LOOP_IDS
do_regional_scores_srs <- function(pH.ls = c(5, 7),
                                   sim_id.ls = c("wt"),
                                   filt_int_id.ls = c("none", "ion"),
                                   loop_keep_ids = DEFAULT_REG_SRS_LOOP_KEEP_IDS,
                                   hack_only_ion = FALSE) {
  elec_ls = list("wt.5" = WT_PH5_ELEC_RCB,
                 "wt.7" = WT_PH7_ELEC_RCB)
  
  results = list()
  
  # Generate loop filter file encoding
  # e.g. if loop_keep_ids = c("LOOP_1", "LOOP_2", "NONE"), the generated
  # string will be "L12N"
  get_loop_filter_fattrib <- function(loop_keep_ids) {
    a = strsplit(x = loop_keep_ids,
                 split = "_",
                 fixed = TRUE)
    f = function(a.charv) {
      return(substr(
        x = a.charv[length(a.charv)],
        start = 1,
        stop = 1
      ))
    }
    a = as.character(unlist(lapply(X = a, FUN = f)))
    return(paste0("L", paste(a, collapse = "")))
  }
  
  loop_filter_fattrib = get_loop_filter_fattrib(loop_keep_ids)
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      regsc_scores_dir = get_regional_scores_dir(pH = pH, sim_id = sim_id)
      # Make sure directory exists
      dir.create(regsc_scores_dir,
                 showWarnings = FALSE,
                 recursive = TRUE)
      base_id = paste0(sim_id, ".", pH)
      
      design = elec_ls[[base_id]]$design
      meta = elec_ls[[base_id]]$meta
      
      for (filt_int_id in filt_int_id.ls) {
        # Force only ions to be reported if only ion-ion interactions are kept
        if (filt_int_id == "ion") {
          hack_only_ion = TRUE
        }
        
        rdata_path = get_regional_scores_rdata_path(
          regional_scores_dir = regsc_scores_dir,
          filt_int_id = filt_int_id,
          filt_reg_id = "srs"
        )
        
        rdata_path = insert_attrib_in_rdata(rdata_path, loop_filter_fattrib)
        
        if (hack_only_ion) {
          rdata_path = insert_attrib_in_rdata(rdata_path, "hkION")
        }
        
        final_id = paste0(base_id, ".", filt_int_id)
        
        print(paste0("Processing ", rdata_path))
        results[[final_id]] = list(
          regsc = get_regional_scores_srs(
            design = design,
            meta = meta,
            filt_int = get_reg_int_filt(filt_int_id),
            loop_keep_ids = loop_keep_ids,
            hack_only_ion = hack_only_ion,
            rdata_path = rdata_path
          ),
          rdata_path = rdata_path
        )
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

# Patches
do_regional_scores_patch <- function(pH.ls = c(5, 7),
                                     sim_id.ls = c("wt"),
                                     filt_int_id.ls = c("none", "ion")) {
  elec_ls = list("wt.5" = WT_PH5_ELEC_RCB,
                 "wt.7" = WT_PH7_ELEC_RCB)
  
  results = list()
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      regsc_scores_dir = get_regional_scores_dir(pH = pH, sim_id = sim_id)
      # Make sure directory exists
      dir.create(regsc_scores_dir,
                 showWarnings = FALSE,
                 recursive = TRUE)
      base_id = paste0(sim_id, ".", pH)
      
      design = elec_ls[[base_id]]$design
      meta = elec_ls[[base_id]]$meta
      
      for (filt_int_id in filt_int_id.ls) {
        rdata_path = get_regional_scores_rdata_path(
          regional_scores_dir = regsc_scores_dir,
          filt_int_id = filt_int_id,
          filt_reg_id = "patch"
        )
        
        final_id = paste0(base_id, ".", filt_int_id)
        
        print(paste0("Processing ", rdata_path))
        results[[final_id]] = list(
          regsc = get_regional_scores_patch(
            design = design,
            meta = meta,
            filt_int = get_reg_int_filt(filt_int_id = filt_int_id),
            rdata_path = rdata_path
          ),
          rdata_path = rdata_path
        )
        
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

##############################
# Boot

# Master bootstrap data generator
do_reg_boot <- function(pH.ls = c(5, 7),
                        sim_id.ls = c("wt"),
                        filt_int_id.ls = c("none", "ion"),
                        stat_type.ls = c("med_diff", "mea_diff"),
                        regsc_type.ls = c("srs", "patch"),
                        srs_loop_keep_ids = DEFAULT_REG_SRS_LOOP_KEEP_IDS,
                        srs_hack_only_ion = FALSE) {
  oc_ls = list(wt.5 = WT_PH5_OC_RCB,
               wt.7 = WT_PH7_OC_RCB)
  
  stat_ls = list(med_diff = reg_med_diff_stat,
                 mea_diff = reg_mea_diff_stat)
  
  
  regsc_fn_ls = list(
    rp = do_regional_scores_rp,
    lp = do_regional_scores_lp,
    srs = do_regional_scores_srs,
    patch = do_regional_scores_patch
  )
  
  results = list()
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      oc_id = paste0(sim_id, ".", pH)
      oc = oc_ls[[oc_id]]
      
      for (filt_int_id in filt_int_id.ls) {
        for (regsc_type in regsc_type.ls) {
          # Obtain regional scores
          # @HACK - handle extra parameters for single residues scoring
          if (regsc_type == "srs") {
            regsc_ls = regsc_fn_ls[[regsc_type]](
              pH.ls = pH,
              sim_id.ls = sim_id,
              filt_int_id.ls = filt_int_id,
              loop_keep_ids = srs_loop_keep_ids,
              hack_only_ion = srs_hack_only_ion
            )[[1]]
          } else {
            regsc_ls = regsc_fn_ls[[regsc_type]](
              pH.ls = pH,
              sim_id.ls = sim_id,
              filt_int_id.ls = filt_int_id
            )[[1]]
          }
          
          for (stat_type in stat_type.ls) {
            stat = stat_ls[[stat_type]]
            
            sum_rdata_path = get_reg_boot_sum_path(regsc_rdata_path = regsc_ls$rdata_path,
                                                   stat_type = stat_type)
            
            final_id = paste0(oc_id,
                              ".",
                              filt_int_id,
                              ".",
                              regsc_type,
                              ".",
                              stat_type)
            
            print(paste0("Processing ", sum_rdata_path))
            results[[final_id]] = list(
              boot = get_reg_boot(
                reg_scores_mat = regsc_ls$regsc$design,
                meta = regsc_ls$regsc$meta,
                oc = oc,
                stat = stat,
                sum_rdata_path = sum_rdata_path
              ),
              sum_rdata_path = sum_rdata_path
            )
            
          } # end iteration over stat type
        } # end iteration of region score type (regional region filter)
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
  
  return(results)
}

##############################
# Boot

# Master FCR bar plot generator
do_plot_reg_fcr_bar <- function(pH.ls = c(5, 7),
                                sim_id.ls = c("wt"),
                                filt_int_id.ls = c("none", "ion"),
                                stat_type.ls = c("mea_diff", "med_diff"),
                                regsc_type.ls = c("rp", "lp", "patch", "srs"),
                                srs_loop_keep_ids = DEFAULT_REG_SRS_LOOP_KEEP_IDS,
                                srs_hack_only_ion = FALSE,
                                plot_fill_lo = "cadetblue",
                                plot_fill_hi = "greenyellow",
                                plot_col_err = "darksalmon",
                                plot_top_k = 25,
                                collect_plots = FALSE) {
  library(ggplot2)
  
  regsc_to_feat_name_fnc = list(
    "rp" = reg_bar_feat_name_rp_desc,
    "lp" = reg_bar_feat_name_lp_desc,
    "srs" = reg_bar_feat_name_srs_desc,
    "patch" = reg_bar_feat_name_patch_desc
  )
  
  regsc_to_xlab = list(
    "rp" = "OmpG Reside-Residue Interaction",
    "lp" = "OmpG Loop Interaction",
    "srs" = "OmpG Residue",
    "patch" = "OmpG Patch-to-Residue Interaction"
  )
  
  stat_type_to_desc = list("mea_diff" = "Mean",
                           "med_diff" = "Median")
  
  out.ls = list()
  
  for (pH in pH.ls) {
    for (sim_id in sim_id.ls) {
      for (filt_int_id in filt_int_id.ls) {
        for (regsc_type in regsc_type.ls) {
          feat_name_fnc = regsc_to_feat_name_fnc[[regsc_type]]
          xlab = regsc_to_xlab[[regsc_type]]
          
          for (stat_type in stat_type.ls) {
            bt.res = do_reg_boot(
              pH.ls = pH,
              sim_id.ls = sim_id,
              filt_int_id.ls = filt_int_id,
              stat_type.ls = stat_type,
              regsc_type.ls = regsc_type,
              srs_loop_keep_ids = srs_loop_keep_ids,
              srs_hack_only_ion = srs_hack_only_ion
            )[[1]]
            
            pdf_path = get_ext_path_from_rdata(bt.res$sum_rdata_path, ".pdf")

            stat_type_desc = stat_type_to_desc[[stat_type]]
            
            print(paste0("Processing ", pdf_path))
            p = plot_reg_fcr_bar(
              bt.ci = bt.res$boot,
              feat_name_fnc = feat_name_fnc,
              main = "Net Stability",
              submain = paste0(round(100.0 * bt.res$boot$Conf.Fcr[1], 1),
                               "% FCR, pH ", pH),
              xlab = xlab,
              ylab = paste0(stat_type_desc, " Net Stability"),
              fill_lo = plot_fill_lo,
              fill_hi = plot_fill_hi,
              col_err = plot_col_err,
              top_k = plot_top_k
            )
            
            # Save plot if valid
            if (!is.null(p)) {
              ggplot2::ggsave(filename = pdf_path, plot = p)
              if (collect_plots) {
                out.ls[[paste0(sim_id,
                               ".",
                               pH,
                               ".",
                               get_fid_from_rdata(bt.res$sum_rdata_path))]] = p
              }
              
            } # end check if plot is valid
            
          } # end iteration over stat type
        } # end iteration of region score type (regional region filter)
      } # end iteration over interaction filter
    } # end iteration over sim_id
  } # end iteration over pH
 
  return(out.ls) 
}

# Output sample identifiers based on energy rank and open-close status
# This is meant to generate an easily parsible list to be fed into
# external utilities such as a python script to render each of the
# sample lists through pymol.
do_write_ener_oc_sample_ids <- function() {
  process_pH_and_barrel <-
    function(pH = 5,
             barrel = "merge",
             max_erank = 5000) {
      pH_id = paste0("pH", pH)
      dat = SIM_WT[[pH_id]][[barrel]]
      oc = dat$oc[1:max_erank]
      # Split energy ranks by open v close
      ranks_o = which(oc)
      ranks_o_df = data.frame(NAME = names(ranks_o), RANK = ranks_o)
      rownames(ranks_o_df) = c()
      ranks_c = which(!oc)
      ranks_c_df = data.frame(NAME = names(ranks_c), RANK = ranks_c)
      rownames(ranks_c_df) = c()
      # Determine output paths
      target_dir = get_regional_scores_dir(pH = pH)
      target_ofid = paste0("sids.", barrel, ".open.csv")
      target_cfid = paste0("sids.", barrel, ".close.csv")
      target_opath = file.path(target_dir, target_ofid)
      target_cpath = file.path(target_dir, target_cfid)
      write.csv(
        x = ranks_o_df,
        file = target_opath,
        quote = FALSE,
        row.names = FALSE
      )
      write.csv(
        x = ranks_c_df,
        file = target_cpath,
        quote = FALSE,
        row.names = FALSE
      )
    }
  
  process_pH_and_barrel(pH = 5, barrel = "merge")
  process_pH_and_barrel(pH = 7, barrel = "merge")
}

# Utility only generates images for JACS paper
do_plot_reg_fcr_bar_JACS <- function() {
  plots = do_plot_reg_fcr_bar(
    pH.ls = c(5, 7),
    sim_id.ls = c("wt"),
    filt_int_id.ls = c("none"),
    stat_type.ls = c("med_diff"),
    regsc_type.ls = c("lp", "srs"),
    srs_loop_keep_ids = c("NONE"),
    srs_hack_only_ion = TRUE,
    plot_fill_lo = "grey",
    plot_fill_hi = "grey",
    plot_col_err = "black",
    plot_top_k = 20,
    collect_plots = TRUE
  )
  
  # Make grid of images
  # http://stackoverflow.com/questions/1249548/side-by-side-plots-with-ggplot2
  DOCS_DIR = file.path(get_script_dir(), "..", "..", "docs")
  JACS_DIR = file.path(DOCS_DIR, "raw", "OmpG.JACS")
  OUT_PATH = file.path(JACS_DIR, "lp.srs.grid.pdf")
  dir.create(JACS_DIR, showWarnings = FALSE)
  
  pg = cowplot::plot_grid(plotlist = plots, labels = c("(a)", "(b)", "(c)", "(d)"))
  print(paste0("Saving ", OUT_PATH))
  cowplot::ggsave(filename = OUT_PATH, plot = pg)
}
