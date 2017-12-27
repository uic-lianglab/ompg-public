###################################################################
# electrostat_globals.R
#
# Global data sets and data paths
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

# Files in this directory are separated by barrel template. All
# hydrogens along with simulated loop side chains (including Cb)
# have been relaxed.
# @return path to relax directory
get_relax_dir <- function() {
  file.path(get_root_dir(), "ompg", "output", "relax")
}

#@return path contain x-ray crystal structure PDBs
get_xray_dir <- function() {
  return(file.path(get_root_dir(),
                   'ompg',
                   'templates'))
}

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "load_data_utils.R"))

###################################################################
# Defaults
###################################################################

# Base output directory
DEF_BASE_OUTPUT_DIR = get_relax_dir()

# Base output directory method
get_def_base_output_dir <- get_relax_dir

# pH level
DEF_PH = 5

# Simulation prefix
DEF_SIM_PREFIX = "wt"

# Barrel identifier
DEF_BARREL = "2iwv"

# Simulation identifier (mutant_barrel)
DEF_SIM_ID = paste0(DEF_SIM_PREFIX, "_", DEF_BARREL)

# Probe radius in Angstroms for open/close analysis
# 2.75 Angstroms = K
# 1.75 Angstroms = Cl
DEF_PROBE = 2.75

# Only keep samples with energy rank less than or equal to
# this value. If -1, all samples are retained.
DEF_MAX_ENERGY_RANK = 5000

# Switch to toggle overwrite behavior. If TRUE, cached rdata
# will be regenerated
DEF_SHOULD_OVERWRITE = FALSE

# TRUE if design and meta info should be loaded for pairwise
# electrostatic interaction data
DEF_SHOULD_LOAD_ELEC = TRUE

# Template PDB information needed for pairwise electrostatic
# interaction data
DEF_TEMPLATE_PDB_INFO = load_pdb_info(file.path(get_xray_dir(),
                                                '2iwv.pdb'))

# An electrostatic interaction pair is discarded if the *total* number
# of observations is below this threshold
DEF_MIN_OBS_COUNT = 100

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "load_data_utils.R"))

###################################################################
# Amino Acid Code Mappings
###################################################################

# Amino Acid 1 letter to 3 letter code mapping
# Example usage: AA1to3["A"], will output "ALA"
AA1to3 = c(
  "ALA",
  "ARG",
  "ASN",
  "ASP",
  "CYS",
  "GLU",
  "GLN",
  "GLY",
  "HIS",
  "ILE",
  "LEU",
  "LYS",
  "MET",
  "PHE",
  "PRO",
  "SER",
  "THR",
  "TRP",
  "TYR",
  "VAL"
)

# Amino Acid 3 letter to 1 letter code mapping
# Example usage: AA3to1["ALA"], will output "A"
AA3to1 = c(
  "A",
  "R",
  "N",
  "D",
  "C",
  "E",
  "Q",
  "G",
  "H",
  "I",
  "L",
  "K",
  "M",
  "F",
  "P",
  "S",
  "T",
  "W",
  "Y",
  "V"
)

# Maps 1-letter code to 3-letter code
names(AA1to3) = AA3to1

# Maps 3-letter code to 1-letter code
names(AA3to1) = AA1to3

###################################################################
# Specific path utilities
###################################################################

# @return path to directory containing energy capture file
get_capt_dir <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                         pH = DEF_PH,
                         sim_id = DEF_SIM_ID) {
  file.path(base_output_dir, "capt", paste0("pH", pH), sim_id)
}

# @return path to energy capture file
get_capt_path <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                          pH = DEF_PH,
                          sim_id = DEF_SIM_ID) {
  file.path(get_capt_dir(
    base_output_dir = base_output_dir,
    pH = pH,
    sim_id = sim_id
  ),
  "ener.csv")
}

# @return path to CASTp base directory
get_castp_dir <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                          pH = DEF_PH,
                          sim_id = DEF_SIM_ID) {
  file.path(base_output_dir,
            "castp",
            paste0("pH", pH),
            sim_id)
}

# @return path to CASTp open/close file
get_castp_path <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                           pH = DEF_PH,
                           sim_id = DEF_SIM_ID,
                           probe = DEF_PROBE) {
  file.path(
    get_castp_dir(
      base_output_dir = base_output_dir,
      pH = pH,
      sim_id = sim_id
    ),
    paste0("oc", probe, ".csv")
  )
}

# @return path to electrostatic interaction meta-data
get_elec_path <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                          pH = DEF_PH,
                          sim_id = DEF_SIM_ID) {
  file.path(
    get_capt_dir(
      base_output_dir = base_output_dir,
      pH = pH,
      sim_id = sim_id
    ),
    "ints_electro.csv"
  )
}

###################################################################
# Simulation loader
###################################################################

# @param base_output_dir - directory containing all outputs
# @param sim_prefix - simulation identifier sans barrel type
# @param probe - CASTp probe radius in Angstroms
# @param max_energy_rank - Only keeps samples with energy rank less
#   than or equal to this value. If -1, all samples are retained.
# @param should_load_elec - TRUE if design and meta info should be
#   loaded for pairwise electrostatic interaction data
# @param template_pdb_info - only needed if should_load_elec is TRUE,
#   used for generating interaction meta info
# @param min_obs_count - only needed if should_load_elec is TRUE,
#   interaction pairs must be observed by at least this count for
#   retention in final data set.
# @param overwrite - TRUE if cached rdata should be regenerated
# @return list with the following members
#   $pH5 - list with members: $b2iww, $b2iwv, $merge each of which
#     is a list with members:
#       $ener - data.frame of energy values with columns
#         NAME, PH, POTENTIAL, ELECT, VDW, BOND, ANGLE, DIHED, IMPRP,
#         BOUNDARY, MISC
#       $oc - logical vector: TRUE is open, FALSE is closed state
#       $elec - list with members:
#         $design - matrix of size n_sample X n_interaction_pair
#           where each entry is electrostatic interaction for that
#           ij pair at that sample
#         $meta - data.frame where each row is meta information for
#           each interaction pair column of design matrix
#   $pH7 - list with same members as "pH5" but for pH 7
load_sim <- function(base_output_dir = DEF_BASE_OUTPUT_DIR,
                     sim_prefix = DEF_SIM_PREFIX,
                     probe = DEF_PROBE,
                     max_energy_rank = DEF_MAX_ENERGY_RANK,
                     should_load_elec = DEF_SHOULD_LOAD_ELEC,
                     template_pdb_info = DEF_TEMPLATE_PDB_INFO,
                     min_obs_count = DEF_MIN_OBS_COUNT,
                     overwrite = DEF_SHOULD_OVERWRITE) {
  # @param barrel - string barrel name
  # @param pH - pH at which data was captured
  # @return energy and open/close states for samples with target
  #   barrel type at parameter pH. Output is a list with members:
  #     $ener - data.frame of energy values for each sample
  #     $oc - logical vector with TRUE for open state, FALSE for close
  load_barrel <- function(barrel = DEF_BARREL,
                          pH = DEF_PH) {
    sim_id = paste0(sim_prefix, "_", barrel)
    # Determine path to captured energy data.frame
    ener_path = get_capt_path(base_output_dir = base_output_dir,
                              pH = pH,
                              sim_id = sim_id)
    # Load energy data.frame
    out = list()
    out[["ener"]] = ener_csv_2_rdata(csv_path = ener_path,
                                     k = max_energy_rank,
                                     overwrite = overwrite)
    # Determine path to open/close state data
    oc_path = get_castp_path(
      base_output_dir = base_output_dir,
      pH = pH,
      sim_id = sim_id,
      probe = probe
    )
    out[["oc"]] = oc_csv_2_rdata(
      csv_path = oc_path,
      ener_ids = out$ener$NAME,
      overwrite = overwrite
    )
    
    if (should_load_elec) {
      # Determine path to electrostatic interactions
      elec_path = get_elec_path(base_output_dir = base_output_dir,
                                pH = pH,
                                sim_id = sim_id)
      out[["elec"]] = elec_csv_2_rdata(
        csv_path = elec_path,
        template_pdb_info = template_pdb_info,
        ener_ids = out[["ener"]]$NAME,
        oc = out[["oc"]],
        min_obs_count = min_obs_count,
        overwrite = overwrite
      )
    }
    
    return(out)
  }
  
  # Merge barrel outputs into a single list sorted by energy
  merge_barrels <- function(lst_barrel_a, lst_barrel_b) {
    out = list()
    
    # Merge energy data.frame
    ener = rbind(lst_barrel_a$ener,
                 lst_barrel_b$ener)
    # Sort by potential energy
    erank = order(ener$POTENTIAL)
    ener = ener[erank, ]
    # Apply energy cutoff
    if (max_energy_rank > 0) {
      ener = ener[1:max_energy_rank,]
    }
    rownames(ener) = c()
    
    # Merge open/close state
    oc = c(lst_barrel_a$oc,
           lst_barrel_b$oc)
    # Sort by potential energy
    oc = oc[ener$NAME]
    stopifnot(names(oc) == ener$NAME)
    
    out[["ener"]] = ener
    out[["oc"]] = oc
    
    # Merge electrostatic interactions if present
    if (("elec" %in% names(lst_barrel_a)) &&
        ("elec" %in% names(lst_barrel_b))) {
      # Default to first barrel argument
      elec.design = lst_barrel_a$elec$design
      elec.meta = lst_barrel_a$elec$meta
      # Apply merged energy cutoff to interaction data
      samp.names.a = ener$NAME[ener$NAME %in% lst_barrel_a$ener$NAME]
      samp.names.b = ener$NAME[ener$NAME %in% lst_barrel_b$ener$NAME]
      # Merge if top-ranked elements are found in both barrel types:
      if ((length(samp.names.a) > 0) &&
          (length(samp.names.b) > 0)) {
        elec.design.a = lst_barrel_a$elec$design[samp.names.a,]
        elec.design.b = lst_barrel_b$elec$design[samp.names.b,]
        # Strip NA columns
        # http://stackoverflow.com/questions/15968494/how-to-delete-columns-that-contain-only-nas
        elec.design.a = elec.design.a[, colSums(is.na(elec.design.a)) != nrow(elec.design.a)]
        elec.design.b = elec.design.b[, colSums(is.na(elec.design.b)) != nrow(elec.design.b)]
        # Determine union and global order of columns
        int.names.all = union(colnames(elec.design.a), colnames(elec.design.b))
        int.order = order(as.numeric(int.names.all))
        int.names.all = int.names.all[int.order]
        # Insert missing columns
        # http://stackoverflow.com/questions/21968126/insert-nonexistent-columns-in-matrix-or-dataframe-in-given-order
        elec.design.a = elec.design.a[, match(int.names.all, colnames(elec.design.a))]
        colnames(elec.design.a) = int.names.all
        elec.design.b = elec.design.b[, match(int.names.all, colnames(elec.design.b))]
        colnames(elec.design.b) = int.names.all
        # Rbind design
        stopifnot(all(colnames(elec.design.a) == colnames(elec.design.b)))
        elec.design = rbind(elec.design.a, elec.design.b)
        # Sort rows by energy
        elec.design = elec.design[ener$NAME,]
        # Combine meta data (strip summary statistics as they are tricky to combine)
        meta.keep.cols = c(
          "RESA.NO",
          "RESB.NO",
          "RESA.NAME",
          "RESB.NAME",
          "RESA.LOOP",
          "RESB.LOOP",
          "N",
          "N.OPEN",
          "N.CLOSE",
          "ION.ION",
          "POL.POL",
          "LOOP.LOOP",
          "LOOP.ION",
          "LOOP.POL",
          "RESA.ION",
          "RESB.ION",
          "RESA.POL",
          "RESB.POL"
        )
        elec.meta = rbind(lst_barrel_a$elec$meta[, meta.keep.cols], lst_barrel_b$elec$meta[, meta.keep.cols])
        elec.meta = unique(elec.meta)
        elec.meta = elec.meta[int.names.all, ]
      } else {
        # Only one barrel type is represented in top energy samples.
        # Since default is first barrel (a), we only need to check
        # if it's empty, in which case, switch to second barrel (b).
        if (length(samp.names.b) > 0) {
          stopifnot(length(samp.names.b) == max_energy_rank)
          stopifnot(length(samp.names.a) == 0)
          elec.design = lst_barrel_b$elec$design
          elec.meta = lst_barrel_b$elec$meta
        }
      }
      stopifnot(nrow(elec.design) == max_energy_rank)
      stopifnot(ncol(elec.design) == nrow(elec.meta))
      stopifnot(all(colnames(elec.design) == rownames(elec.meta)))
      out[["elec"]] = list(design = elec.design, meta = elec.meta)
    }
    
    return(out)
  }
  
  out = list()
  pHs = c(5, 7)
  barrels = c("2iwv", "2iww")
  
  for (pH in pHs) {
    pH_id = paste0("pH", pH)
    out[[pH_id]] = list()
    
    for (barrel in barrels) {
      barrel_id = paste0("b", barrel)
      out[[pH_id]][[barrel_id]] = load_barrel(barrel = barrel,
                                              pH = pH)
      
    } # end iteration over barrels
    
    # Create merge set
    out[[pH_id]][["merge"]] = merge_barrels(out[[pH_id]]$b2iww,
                                            out[[pH_id]]$b2iwv)
    
  } # end iteration over pH
  
  return(out)
}
