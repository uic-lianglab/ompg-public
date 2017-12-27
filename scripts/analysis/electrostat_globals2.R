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

###################################################################
# Data set utilities
###################################################################

source(file.path(get_script_dir(), "electrostat_merge_utilities.R"))

###################################################################
# Load binary data
###################################################################

# Switch to toggle overwrite behavior
SHOULD_OVERWRITE = FALSE

##############################
# Simulation data sets

# @Warning, template_pdb_info must be changed for mutant simulations
SIM_WT = load_sim(
  base_output_dir = get_relax_dir(),
  sim_prefix = "wt",
  probe = DEF_PROBE,
  max_energy_rank = DEF_MAX_ENERGY_RANK,
  should_load_elec = TRUE,
  template_pdb_info = DEF_TEMPLATE_PDB_INFO,
  min_obs_count = DEF_MIN_OBS_COUNT,
  overwrite = SHOULD_OVERWRITE
)

##############################
# Energy data sets

# WT RCB
WT_PH5_ENER_RCB = SIM_WT$pH5$merge$ener
WT_PH7_ENER_RCB = SIM_WT$pH7$merge$ener

##############################
# Open-close data sets

# WT RCB
WT_PH5_OC_RCB = SIM_WT$pH5$merge$oc
WT_PH7_OC_RCB = SIM_WT$pH7$merge$oc

##############################
# Electrostatic interactions

# WT RCB
WT_PH5_ELEC_RCB = SIM_WT$pH5$merge$elec
WT_PH7_ELEC_RCB = SIM_WT$pH7$merge$elec
