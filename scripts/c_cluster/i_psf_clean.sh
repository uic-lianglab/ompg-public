#!/bin/bash

# PSF files take up too much disk space. This script removes duplicate .psf
# files within psf folders and writes a mapping from original file name
# to unique template file name in parent folder.

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
# Path to py script for removing duplicates
REMOVE_DUPS_SCRIPT_PATH=$ROOT_DIR/scripts/remove_duplicates.py
# Parent folder to parse for psf files
BASE_TOPO_INPUT_DIR=$BASE_OUTPUT_DIR/multi_loop/topo

# Path for storing logs
BASE_LOG_DIR=$DATA_DIR/logs
BASE_PSF_CLEAN_LOG_DIR=$BASE_LOG_DIR/multi_loop/psf_clean

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
MAX_NPROC=25

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
    
    # Count number of target files in child directory
    has_target=`find $d -maxdepth 1 -type f -name '*.psf' | wc -l`
    # If count > 0, then this directory contains files of interest
    # So, skip all folder with no files of interest
    if [ $has_target == 0 ]; then 
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Full path to input directory
    psf_dir=$(cd "$d"; pwd)
    # Path to output mapping
    psf_templates_path="$psf_dir/../psf_templates.csv"
    # Build command line
    cmd="python $REMOVE_DUPS_SCRIPT_PATH $psf_dir $psf_templates_path 1"
    echo "$(timestamp): $cmd"
    
    # Make log directory
    log_dir="$BASE_PSF_CLEAN_LOG_DIR/$d"
    mkdir -p $log_dir
    log_path="$log_dir/psf_clean.txt"
    
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
        sleep 10s
    done # end queue spinning
done # end iteration over target folders

popd

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
