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

# Files in this directory are merged 2iwv and 2iww templates,
# but c_beta carbon atom remains *fixed* during minimization
# DEPRECATED! Use relax_c_beta instead
# @return path to merge directory
get_merge_dir <- function() {
  file.path(get_root_dir(), "ompg", "output", "merge")
}

# Files in this directory are merged 2iwv and 2iww templates,
# but c_beta carbon is *not* fixed during minimization
# @return path to relax_c_beta directory
get_relax_c_beta_dir <- function() {
  file.path(get_root_dir(), "ompg", "output", "relax_c_beta")
}


# @return path to directory containing energy capture file
get_capt_dir <- function(base_output_dir = get_relax_c_beta_dir(),
                         pH = 5,
                         sim_id = "wt") {
  file.path(base_output_dir, "capt", paste0("pH", pH), sim_id)
}

# @return path to energy capture file
get_capt_path <- function(base_output_dir = get_relax_c_beta_dir(),
                          pH = 5,
                          sim_id = "wt") {
  file.path(get_capt_dir(
    base_output_dir = base_output_dir,
    pH = pH,
    sim_id = sim_id
  ),
  "ener.csv")
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

# @return path to electrostatic interaction meta-data
get_elec_path <- function(base_output_dir = get_relax_c_beta_dir(),
                          pH = 5,
                          sim_id = "wt") {
  file.path(
    get_capt_dir(
      base_output_dir = base_output_dir,
      pH = pH,
      sim_id = sim_id
    ),
    "ints_electro.csv"
  )
}

# Default base output directory
get_def_base_output_dir <- get_relax_c_beta_dir

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
WT_TEMPLATE_PDB_PATH = file.path(get_root_dir(), 'ompg', 'templates', '2iwv.pdb')
WT_PH5_ELEC_CSV_PATH_RCB = get_elec_path(get_relax_c_beta_dir(), 5, "wt")
WT_PH7_ELEC_CSV_PATH_RCB = get_elec_path(get_relax_c_beta_dir(), 7, "wt")

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

##############################
# Template reference PDBs

# WT
WT_TEMPLATE_PDB_INFO = load_pdb_info(WT_TEMPLATE_PDB_PATH)

##############################
# Electrostatic interactions

# An electrostatic interaction pair is discarded if the *total* number of
# observations is below this threshold
MIN_OBS_COUNT = 100

# WT RCB
WT_PH5_ELEC_RCB = elec_csv_2_rdata(
  WT_PH5_ELEC_CSV_PATH_RCB,
  WT_TEMPLATE_PDB_INFO,
  WT_PH5_ENER_RCB$NAME,
  WT_PH5_OC_RCB,
  MIN_OBS_COUNT,
  SHOULD_OVERWRITE
)

WT_PH7_ELEC_RCB = elec_csv_2_rdata(
  WT_PH7_ELEC_CSV_PATH_RCB,
  WT_TEMPLATE_PDB_INFO,
  WT_PH7_ENER_RCB$NAME,
  WT_PH7_OC_RCB,
  MIN_OBS_COUNT,
  SHOULD_OVERWRITE
)
