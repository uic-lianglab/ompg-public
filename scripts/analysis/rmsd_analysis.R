###################################################################
# rmsd_analysis.R
#
# Analysis relating to root-mean-square distance (RMSD)
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
# Libraries
###################################################################

library(Rpdb)
library(matrixStats)

###################################################################
# Globals
###################################################################

# Load globals
source(file.path(get_script_dir(), "electrostat_globals2.R"))

###################################################################
# More paths
###################################################################

# @return path to directory containing PDBs for topology analysis
get_castp_pdb_dir <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                              pH = DEF_PH,
                              sim_id = DEF_SIM_ID) {
  return(file.path(
    get_castp_dir(
      base_output_dir = base_output_dir,
      pH = pH,
      sim_id = sim_id
    ),
    "pdb"
  ))
}

# @return output dir for
get_rmsd_dir <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                         sim_prefix = DEF_SIM_PREFIX) {
  return(file.path(base_output_dir,
                   "rmsd",
                   sim_prefix))
}

###################################################################
# RMSD analysis
###################################################################

# @param loop_ids - if empty, all valid loop identifiers are
#   selected. Valid loop identifiers are any combination of 1, 2, 3,
#   5, 6, and 7. For instance, if loop_ids = c(6, 7), then only
#   residue identifiers for loops 6 and 7 will be returned.
# @param template_pdb_info - mapping from residue identifiers to
#   loop identifiers
# @return residue idenitifiers for target loop regions
get_target_residues <- function(loop_ids = c(),
                                template_pdb_info = DEF_TEMPLATE_PDB_INFO) {
  if (length(loop_ids) == 0) {
    return(which(
      grepl(
        pattern = "LOOP",
        x = template_pdb_info$loop_map,
        fixed = TRUE
      )
    ))
  } else {
    # Creates regular expression. Example, if loop_ids = c(1,6),
    # then pattern will equal "LOOP_1|LOOP_6"
    pattern = paste0("LOOP_", loop_ids, collapse = "|")
    return(which(
      grepl(
        pattern = pattern,
        x = template_pdb_info$loop_map,
        fixed = FALSE
      )
    ))
  }
}

# @return protein atom coordinates data.frame
load_pdb <- function(pdb_path, verbose = TRUE) {
  if (verbose) {
    print(paste0("Loading ", pdb_path))
  }
  prot = read.pdb(pdb_path)
  prot = prot$atoms
  # Strip heteroatoms
  prot = prot[prot$recname == "ATOM", ]
  prot = data.frame(
    elename = prot$elename,
    resname = prot$resname,
    resid = prot$resid,
    x = prot$x1,
    y = prot$x2,
    z = prot$x3
  )
  return(prot)
}

load_pdb_dir <- function(pdb_dir = get_castp_pdb_dir(),
                         barrel = DEF_BARREL,
                         eranks = WT_PH5_ENER_RCB) {
  print(paste0("Loading PDB directory: ", pdb_dir))
  # Filter by barrel
  barrel_filt = grepl(pattern = barrel, x = eranks$NAME)
  eranks = eranks$NAME[barrel_filt]
  
  prots = list()
  for (pid in eranks) {
    pdb_path = file.path(pdb_dir, paste0(pid, ".pdb"))
    if (file.exists(pdb_path)) {
      prots[[pid]] = load_pdb(pdb_path, verbose = FALSE)
    }
  }
  print(paste0("Loaded ", length(prots), " proteins."))
  return(prots)
}

# Computes root-mean-square distance between two proteins
# @param p1 - first protein
# @param p2 - second protein
# @param target_residues - residues to compute RMSDs for
# @return numeric vector with named elements:
#   N - RMSD of N atoms
#   CA - RMSD of CA atoms
#   C - RMSD of C atoms
#   O - RMSD of O atoms
#   ALL - RMSD of all atoms (N,CA,C,O)
rmsd_util <- function(p1, p2, target_residues) {
  p1 = p1[p1$resid %in% target_residues, ]
  p2 = p2[p2$resid %in% target_residues, ]
  bb_atoms = c('N', 'CA', 'C', 'O')
  p1 = p1[p1$elename %in% bb_atoms, ]
  p2 = p2[p2$elename %in% bb_atoms, ]
  stopifnot(all(p1$elename == p2$elename))
  stopifnot(all(p1$resid == p2$resid))
  
  dx2 = p1$x - p2$x
  dx2 = dx2 * dx2
  dy2 = p1$y - p2$y
  dy2 = dy2 * dy2
  dz2 = p1$z - p2$z
  dz2 = dz2 * dz2
  d_prof = data.frame(atom = p1$elename, dist2 = (dx2 + dy2 + dz2))
  
  d = rep(0.0, length(bb_atoms) + 1)
  names(d) = c(bb_atoms, "ALL")
  for (a in bb_atoms) {
    d[a] = sqrt(mean(d_prof$dist2[d_prof$atom == a]))
  }
  d["ALL"] = sqrt(mean(d_prof$dist2))
  return(d)
}

# Computes root-mean-square distances for low energy structures
# relative to the x-ray crystal structures
# @param sim - simulation structure as returned by load_sim()
# @param sim_prefix - name of simulation
# @param base_output_dir - directory to place rmsd results in
# @param pHs - which pH levels to process
# @param barrels - which barrel templates to process
# @param target_residues - which residues to include in RMSD
#   calculation
# @param tag - additional tag appended to output file name
# @param overwrite - if true, overwrites any cached data
# @return list of data.frames with root-mean-square distances
get_rmsd <- function(sim = SIM_WT,
                     sim_prefix = DEF_SIM_PREFIX,
                     base_output_dir = DEF_BASE_OUTPUT_DIR,
                     pHs = c(5, 7),
                     barrels = c("2iwv", "2iww"),
                     target_residues = get_target_residues(),
                     tag = "all",
                     overwrite = FALSE) {
  out_ls = list()
  
  for (barrel in barrels) {
    # Load template protein
    tprot = load_pdb(file.path(get_xray_dir(), paste0(barrel, ".pdb")))
    # Determine simulation identifier
    sim_id = paste0(sim_prefix, "_", barrel)
    for (pH in pHs) {
      print(paste0("----------------------"))
      print(paste0("Computing RMSD for ", barrel, " at pH = ", pH))
      
      # Determine output path to write rdata
      rdata_path = file.path(
        get_rmsd_dir(base_output_dir = base_output_dir, sim_prefix = sim_prefix),
        paste("rmsd", barrel, paste0("pH", pH), tag, "rdata", sep = ".")
      )
      
      rmsd_id = paste(sim_id, pH, tag, sep = ".")
      
      if (file.exists(rdata_path) && !overwrite) {
        print(paste0("Loading cached data at: ", rdata_path))
        out_ls[[rmsd_id]] = load_rdata(rdata_path)
        next
      }
      
      # Determine input protein directory
      pdb_dir = get_castp_pdb_dir(base_output_dir = base_output_dir,
                                  pH = pH,
                                  sim_id = sim_id)
      
      # Load proteins
      prots = load_pdb_dir(
        pdb_dir = pdb_dir,
        barrel = barrel,
        eranks = sim[[paste0("pH", pH)]]$merge$ener
      )
      
      n_prots = length(prots)
      # Pad by +4 to account for mean, meadian, min, max summaries
      zero_vec = rep(0.0, n_prots + 4)
      rmsd_df = data.frame(
        NAME = rep("", n_prots + 4),
        N = zero_vec,
        CA = zero_vec,
        C = zero_vec,
        O = zero_vec,
        ALL = zero_vec,
        stringsAsFactors = FALSE
      )
      
      # Process each protein
      for (i in 1:n_prots) {
        d = rmsd_util(tprot, prots[[i]], target_residues)
        rmsd_df[i, names(d)] = d
        rmsd_df$NAME[i] = names(prots)[i]
      } # end iteration over loaded proteins
      
      # Append summaries
      rmsd_cols = c("N", "CA", "C", "O", "ALL")
      rmsd_mat = as.matrix(rmsd_df[1:n_prots, rmsd_cols])
      # Mean
      rmsd_df$NAME[n_prots + 1] = "MEAN"
      rmsd_df[n_prots + 1, rmsd_cols] = colMeans(x = rmsd_mat)
      # Median
      rmsd_df$NAME[n_prots + 2] = "MEDIAN"
      rmsd_df[n_prots + 2, rmsd_cols] = colMedians(x = rmsd_mat)
      # Min
      rmsd_df$NAME[n_prots + 3] = "MIN"
      rmsd_df[n_prots + 3, rmsd_cols] = apply(X = rmsd_mat,
                                              MARGIN = 2,
                                              FUN = min)
      # Max
      rmsd_df$NAME[n_prots + 4] = "MAX"
      rmsd_df[n_prots + 4, rmsd_cols] = apply(X = rmsd_mat,
                                              MARGIN = 2,
                                              FUN = max)
      
      # Save result
      save_rdata(
        data = rmsd_df,
        rdata_path = rdata_path,
        b_create_subdirs = TRUE
      )
      save_csv(
        x = rmsd_df,
        csv_path = get_csv_path_from_rdata(rdata_path),
        row.names = FALSE,
        col.names = TRUE
      )
    } # end iteration over pH levels
  } # end iteration over barrels
  
  return(out_ls)
}

# Utility for computing RMSD for loop 6 only
get_rmsd_l6 <- function(sim = SIM_WT,
                        sim_prefix = DEF_SIM_PREFIX,
                        base_output_dir = DEF_BASE_OUTPUT_DIR,
                        pHs = c(5, 7),
                        barrels = c("2iwv", "2iww"),
                        overwrite = FALSE) {
  return(
    get_rmsd(
      sim = sim,
      sim_prefix = sim_prefix,
      base_output_dir = base_output_dir,
      pHs = pHs,
      barrels = barrels,
      target_residues = get_target_residues(6),
      tag = "l6",
      overwrite = overwrite
    )
  )
}
