#!/bin/bash 

# Script cleans .propka_input files

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Path to base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop

###################################################
# Loop invariant arrays
###################################################

# Simulation identifiers
SIM_ID=("wt")

# The template pdb structures to use
# 2iwv -> open conformation crystal structure captured at pH 7
# 2iww -> closed conformation crystal structure captured at pH 5
TEMPLATE=("2iwv" "2iww")

###################################################
# Timestamp utils
###################################################

# From http://stackoverflow.com/questions/17066250/create-timestamp-variable-in-bash-script
timestamp() {
  date +"%T"
}

###################################################
# Clean propka folders
###################################################

# Iterate over simulation identifiers(i.e. mutations and wild types)
for ix_sim_id in "${!SIM_ID[@]}"
do

    # Name of simulation
    sim_id="${SIM_ID[$ix_sim_id]}"

    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): simulation: $sim_id"
    
    # Iterate over template PDBs
    for template in "${TEMPLATE[@]}"
    do
        # Prefix for this (simulation, template) pair
        sim_template_prefix="$sim_id"_$template
        
        # directory to write propka data
        propka_dir=$BASE_OUTPUT_DIR/$sim_template_prefix/propka
        
        # Clean directory if it exists
        if [ -d "$propka_dir" ]; then
            # Switch to output directory
            pushd .
            cd $propka_dir
            
            # Process each .propka_input file - iterate over individual files as
	        # otherwise can receive 'argument list too long error'
	        # http://stackoverflow.com/questions/11289551/argument-list-too-long-error-for-rm-cp-mv-commands
	        for file in $propka_dir/*.propka_input;
	        do
	            rm $file
	        done
	        popd
	        
        fi # end check if propka directory exists
        
    done # end iteration over templates

done # end iteration over sim ids

echo "Finished."
