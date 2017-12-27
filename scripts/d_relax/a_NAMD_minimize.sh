#!/bin/bash

# Will parse psf mapping and feed (.pdb, .psf) pairs
# to python script which calls NAMD minimization routine
# Note: the specific atoms that are held fixed are
# defined in NAMD_minimize.py fix_loop_bb()

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
BASE_OUTPUT_DIR=$DATA_DIR/output
# Base directory for minimization files
BASE_MIN_DIR=$BASE_OUTPUT_DIR/relax/min
# Parent folder to parse for psf files
BASE_TOPO_INPUT_DIR=$BASE_OUTPUT_DIR/multi_loop/topo
# Path to py script for atom potential minimization
PY_SCRIPT_PATH=$ROOT_DIR/scripts/NAMD_minimize.py

# Name of file containing psf mappings
PSF_MAP_NAME=psf_templates.csv

# Path for storing logs
BASE_LOG_DIR=$DATA_DIR/logs
BASE_MIN_LOG_DIR=$BASE_LOG_DIR/relax/min

###################################################
# Completed data sets to ignore
###################################################

# Set of previously completed loops which are ignored by script
# To ignore an entire loop, add element "loop_<#>_<mut>"
# - example, to ignore loop 1 wildtype, add element: "loop_1_wt"
# To ignore only a template within a loop, add element "loop_<#>_<mut>/<template>"
# - example, to ignore loop 1 wildtype at 2iwv template, add element: "loop_1_wt/2iwv"
# To ignore a pH folder within merge, add element "pH<#>/<sim_id>"
# - example,to ignore wildtype at pH 5, add element: "pH5/wt"
COMPLETED=()

# Check if any array element is a substring of parameter
# http://stackoverflow.com/questions/14366390/bash-if-condition-check-if-element-is-present-in-array
# http://stackoverflow.com/questions/7109720/behavior-of-return-statement-in-bash-functions
array_contains_substring_of () { 
    local master=$1; shift
    for substring; do
        if [[ "$master" =~ "$substring" ]]; then
            # Note that echo is actually the return value!
            echo 1
            return
        fi
    done
    echo 0
    return
}

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
MAX_NPROC=37

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

pushd ./

cd $BASE_TOPO_INPUT_DIR

# Iterate over subdirectories containing psf files
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Skip folder if no mapping
    # http://www.shellhacks.com/en/HowTo-Check-If-a-File-Exists
    psf_map_path="$d/$PSF_MAP_NAME"
    if [ ! -f "$psf_map_path" ]; then continue; fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Full path to input directory
    infulld=$(cd "$d"; pwd)
        
    # Output directory
    min_dir=$BASE_MIN_DIR/$d
    mkdir -p $min_dir
    min_dir=$(cd "$min_dir"; pwd)
    
    # Make log directory
    log_dir="$BASE_MIN_LOG_DIR/$d"
    mkdir -p $log_dir
    
    # Process list file, extract paths
    line_count=0
    while read line
    do
        # Update line count
        ((line_count++))
        
        # Below code splits a CSV tuple into an array
        # - first element contains the input psf name
        # - second element contains the mapped template psf to actually use
        # http://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
        # http://www.tldp.org/LDP/abs/html/internalvariables.html
        # IFS stands for input field separator
        IFS=',' read -r -a psf_map <<< "$line"
        
        base_name_no_ext="$(basename ${psf_map[0]} .psf)"
        
        in_coord_path="$infulld/pdb/$base_name_no_ext.pdb"
        in_psf_path="$infulld/psf/${psf_map[1]}"
        
        # Build command line
        cmd="python $PY_SCRIPT_PATH $min_dir $in_coord_path $in_psf_path"
        echo "$line_count: $(timestamp): $cmd"
        
        log_path="$log_dir/$base_name_no_ext.txt"
        
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
            sleep 2s
        done # end queue spinning
        
    done < $psf_map_path
done # end iteration over target folders

popd

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
