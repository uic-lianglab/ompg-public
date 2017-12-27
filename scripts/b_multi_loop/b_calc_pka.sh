#!/bin/bash 
# Script calls propka to calculate pKa on multiloop samples

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
# Path to PROPKA script
PROPKA_PATH=$ROOT_DIR/tools/propka-3.1/propka31.py

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
MAX_NPROC=17

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
        # Directory to write propka data
        propka_dir=$BASE_OUTPUT_DIR/$sim_template_prefix/propka
        # Make propka directory
        mkdir -p $propka_dir
        # Make log directory
        propka_logdir=$LOG_DIR/$sim_template_prefix/propka
        mkdir -p $propka_logdir
        
        echo "+++++++++++++++++++++++++++++"
        echo "$(timestamp): Applying template: $template"
        echo "$(timestamp): pdb dir: $pdb_dir"
        echo "$(timestamp): propka dir: $propka_dir"
        echo ""
        
        # Switch to output directory
        pushd .
        cd $propka_dir
        
        # Process each file in target PDB directory
        file_count=0
        for file in $pdb_dir/*.pdb;
        do
            # Update file count
            ((file_count++))
            
            # Skip if pKa already computed
            base_name_no_ext="$(basename $file .pdb)"
            pka_file="$base_name_no_ext".pka
            if [ -f "$pka_file" ]; then
                echo "$file_count: $(timestamp): Skipping $file as pKa already computed."
                continue
            fi
            
            # Build command line
            cmd="python $PROPKA_PATH $file"
            echo "$file_count: $(timestamp): $cmd"
            pdb=$(basename "$file")
            log_path="$propka_logdir/$pdb.propka.log.txt"
            
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
            done
        done
        popd
    done # End iteration over templates

done # End iteration over sim ids

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
