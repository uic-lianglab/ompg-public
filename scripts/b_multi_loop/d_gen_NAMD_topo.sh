#!/bin/bash

# Script calls NAMD psfgen with patch_topology.tcl
# The patch_topology.tcl script will protonate the input
# pdb according to the parameter pH and associated
# PROPKA pKa file. The final output will be a .psf and
# .coor/.pdb suitable for hydrogen relaxation with NAMD.

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
# Path for storing logs
LOG_DIR=$DATA_DIR/logs/multi_loop
# Path to NAMD psfgen
PSFGEN_PATH=$ROOT_DIR/tools/namd/bin/nix/psfgen
# Path to patch topology tcl script
PATCH_TOPOLOGY_TCL_PATH=$ROOT_DIR/tools/namd/scripts/patch_topology.tcl

###################################################
# Loop invariant arrays
###################################################

# Simulation identifiers
SIM_ID=("wt")

# The template pdb structures to use
# 2iwv -> open conformation crystal structure captured at pH 7
# 2iww -> closed conformation crystal structure captured at pH 5
TEMPLATE=("2iwv" "2iww")

# pH levels to simulate
PH_LEVEL=("5" "7")

###################################################
# Timestamp utils
###################################################

# From http://stackoverflow.com/questions/17066250/create-timestamp-variable-in-bash-script
timestamp() {
  date +"%T"
}

###################################################
# JOBQUEUE UTILS
###################################################

# From https://pebblesinthesand.wordpress.com/2008/05/22/a-srcipt-for-running-processes-in-parallel-in-bash/

NUM=0
QUEUE=""
MAX_NPROC=18

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}

function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

###################################################
# Queue jobs
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
        
        # Directory containing generated loop pdbs
        pdb_dir=$BASE_OUTPUT_DIR/$sim_template_prefix/pdb
        # Directory containing propka data
        propka_dir=$BASE_OUTPUT_DIR/$sim_template_prefix/propka
        
        echo "+++++++++++++++++++++++++++++"
        echo "$(timestamp): Applying template: $template"
        echo "$(timestamp): pdb dir: $pdb_dir"
        echo "$(timestamp): propka dir: $propka_dir"
        echo ""
        
        # Iterate over pH levels
        for ix_pH in "${!PH_LEVELS[@]}"
        do
            pH="${PH_LEVEL[$ix_pH]}"
            # directory to write topology and coordinates data
            pH_dir="$BASE_OUTPUT_DIR/$mut_template_prefix/pH$pH/prot_raw"
            mkdir -p $pH_dir
            # Make log directory
            pH_logdir="$LOG_DIR/$mut_template_prefix/pH$pH/prot_raw"
            mkdir -p $pH_logdir
            
            echo "-----------------------------"
            echo "$(timestamp): pH: $pH"
            echo "$(timestamp): target dir: $pH_dir"
            echo ""
            
            out_coords_dir=$pH_dir/pdb
            mkdir -p $out_coords_dir
            
            out_psf_dir=$pH_dir/psf
            mkdir -p $out_psf_dir
            
            # Process pdb files with associated .pka files
            abs_pdb_dir=$(cd "$pdb_dir"; pwd)
            abs_propka_dir=$(cd "$propka_dir"; pwd)
            for propka_file_path in $abs_propka_dir/*.pka;
            do
                # Extract pdb name
                pdb_name=$(basename "$propka_file_path" .pka)
                pdb_file_path=$abs_pdb_dir/$pdb_name.pdb
                # Path to output protonated coordinates
                out_coords_path=$out_coords_dir/$pdb_name.pH$pH.pdb
                # Path to output protonated topology psf
                out_psf_path=$out_psf_dir/$pdb_name.pH$pH.psf
                # Build command line
                cmd="$PSFGEN_PATH $PATCH_TOPOLOGY_TCL_PATH"
                cmd="$cmd $pdb_file_path $propka_file_path $pH"
                cmd="$cmd $out_coords_path $out_psf_path"
                echo "$(timestamp): $cmd"
                
                log_path=$pH_logdir/$pdb_name.pH$pH.log.txt
                
                # Run command in background
                # http://www.cyberciti.biz/faq/redirecting-stderr-to-stdout/
                $cmd &> $log_path &
                # Capture process id of last process ran
                PID=$!
                # Store process id
                queue $PID
                
                # Check if we need to spin until at least one of the jobs finishes
                while [ $NUM -ge $MAX_NPROC ]; do
                    checkqueue
                    # Sleep for *s -> seconds, *m -> minutes, *h -> hours
                    sleep 1s
                done # end queue spinning
            done # end iteration over pdb files
            
        done # end iteration over pH
    done # end iteration over templates
done # end iteration over sim ids

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
