#!/bin/bash 

# Script exports multi-loop samples to CSV file for clustering

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Path to loop regions directory
LOOP_REGION_DIR=$DATA_DIR/loop_regions
# Path to base loop PDBs directory
BASE_LOOP_DIR=$DATA_DIR/output/multi_loop
# Path to base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop
# Path for storing logs
LOG_DIR=$DATA_DIR/logs/multi_loop
# Path to executable directory
EXE_DIR=$ROOT_DIR
# Name of executable
EXE_NAME=mDisgro

###################################################
# Arrays
###################################################

# Simulation identifiers
SIM_ID=("wt")

# The template pdb structures to use
# 2iwv -> open conformation crystal structure captured at pH 7
# 2iww -> closed conformation crystal structure captured at pH 5
TEMPLATE=("2iwv" "2iww")

###################################################
# Scalars
###################################################

# Maximum number of conformations to export for clustering
# A value <= 0 results in no limit
MAX_CONF=-1

# If 1, then only backbone atoms (within loop regions) are
# exported
BBONE_ONLY=1

# File specifying loop regions
LOOP_REGION_PATH="$LOOP_REGION_DIR/loops_1_2_3_5_6_7.txt"

###################################################
# Timestamp utils
###################################################

# From http://stackoverflow.com/questions/17066250/create-timestamp-variable-in-bash-script
timestamp() {
  date +"%T"
}

###################################################
# Queue jobs
###################################################

log_dir="$LOG_DIR/clust/a_csv_pdbs"
mkdir -p $log_dir

# Iterate over simulation identifiers(i.e. mutations and wild types)
for ix_sim_id in "${!SIM_ID[@]}"
do

    # Name of simulation
    sim_id="${SIM_ID[$ix_sim_id]}"
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): simulation: $sim_id"

    for ix_template in "${!TEMPLATE[@]}"
    do
        # Name of template
        template="${TEMPLATE[$ix_template]}"

        # Prefix for this (simulation, template) pair
        sim_template_prefix="$sim_id"_$template

        # Directory containing generated loop pdbs
        pdb_dir=$BASE_LOOP_DIR/$sim_template_prefix/pdb
        echo ""
        echo "$(timestamp): template: $template"

        # Output directory
        out_dir="$BASE_OUTPUT_DIR/clust/a_csv_pdbs/$sim_template_prefix"
        mkdir -p $out_dir

        # Output files
        out_csv_path="$out_dir/$sim_template_prefix.csv"
        out_ls_path="$out_dir/$sim_template_prefix.ls.txt"

        # Logs
        log_path="$log_dir/$sim_template_prefix.log.txt"

        # Build command
        cmd="./$EXE_NAME -exp_cl_csv_pdbs $pdb_dir $out_csv_path $out_ls_path $MAX_CONF $BBONE_ONLY -multiloops $LOOP_REGION_PATH"

        # Echo command
        echo "$(timestamp): Log path: $log_path"
        echo "$(timestamp): Running cmd:"
        echo "$cmd"
                    
        # Switch to executable dir and run command
        pushd ./
        cd $EXE_DIR
        $cmd > $log_path
        popd

    done # end iteration over templates
 
done #end iteration over simulation identifiers
       
echo "$(timestamp): Export finished."

