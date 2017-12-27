#!/bin/bash 

# Script calls python cluster resampler for each (.c2s, .ex)
# in target directory. The resampler only generates a list of
# indices into the raw fragment library that are to be selected
# for the final "resampled" library

###################################################
# Script paths
###################################################

# Directory containing this script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Script for cluster resampling
RS_SCRIPT_NAME=cluster_id_resampler.py
RS_SCRIPT_PATH=$ROOT_DIR/scripts/$RS_SCRIPT_NAME
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Output directories
BASE_OUTPUT_DIR=$DATA_DIR/output
BASE_FRAG_LIB_INPUT_DIR=$BASE_OUTPUT_DIR/frag_libs/f_clust_libs
BASE_FRAG_LIB_OUTPUT_DIR=$BASE_OUTPUT_DIR/frag_libs/g_resample_ids
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_FRAG_LIB_LOG_DIR=$BASE_LOG_DIR/frag_libs/g_resample_ids

###################################################
# Completed data sets to ignore
###################################################

# Set of previously completed loops which are ignored by script
# To ignore an entire loop, add element "loop_<#>_<mut>"
# - example, to ignore loop 1 wildtype, add element: "loop_1_wt"
# To ignore only a template within a loop, add element "loop_<#>_<mut>/<template>"
# - example, to ignore loop 1 wildtype at 2iwv template, add element: "loop_1_wt/2iwv"
COMPLETED_LOOPS=()

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
MAX_NPROC=15

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
# Loop invariant arguments to resampling tool
###################################################

# Minimum number of samples from each cluster
MIN_COVERAGE=10
# Minimum total number of samples in final fragment library
MIN_OUTPUT_COUNT=10000

###################################################
# Process c2s files in directory
###################################################

pushd ./

cd $BASE_FRAG_LIB_INPUT_DIR

# Iterate over subdirectories containing c2s.csv mappings
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED_LOOPS[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of c2s.csv files in child directory
    has_lib=`ls -1 $d/*.c2s.csv 2>/dev/null | wc -l`
    # If count > 0, then this directory is a frag lib
    # So, skip all folder with no c2s.csv mappings
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Directory to write cluster results
    outdir="$BASE_FRAG_LIB_OUTPUT_DIR/$d"
    echo "$(timestamp): Creating $outdir"
    mkdir -p "$outdir"
    # Capture absolute path to output directory
    outdir=$(cd "$outdir"; pwd)
    
    # Directory to write logs
    logdir="$BASE_FRAG_LIB_LOG_DIR/$d"
    mkdir -p $logdir
    
    # Full path to input directory
    infulld=$(cd "$d"; pwd)
    
    # Process each library in directory
    for f in $d/*.c2s.csv
    do
        # I/O data path arguments to clustering tool
        in_c2s_path="$infulld"/"$(basename $f)"
        # Path to exemplar indices
        in_ex_path="$infulld"/"$(basename $f .c2s.csv)".ex.csv
        # Output path for resample indices
        out_rs_path="$outdir"/"$(basename $f .c2s.csv)".rs.csv
        
        logpath="$logdir"/"$(basename $f .c2s.csv)".rs.txt
        
        # Build command string
        cmd="python $RS_SCRIPT_PATH"
        cmd="$cmd -c2s $in_c2s_path"
        cmd="$cmd -ex $in_ex_path"
        cmd="$cmd -rs $out_rs_path"
        cmd="$cmd -mcov $MIN_COVERAGE"
        cmd="$cmd -mcout $MIN_OUTPUT_COUNT"
        
        # Switch to script dir
        pushd ./
        cd $SCRIPT_DIR
        echo $cmd
        # Run command in background
        $cmd > $logpath &
        # Capture process id of last process ran
        PID=$!
        # Store process id
        queue $PID
        popd
        
        # Check if we need to spin until at least one of the jobs finishes
        while [ $NUM -ge $MAX_NPROC ]; do
            checkqueue
            # Sleep for x seconds
            sleep 5s
        done # end spin until queue frees up
        
    done # end iteration over cluster mappings in current directory
    
done # end iteration over all child directories containing cluster mappings

popd

###################################################
# Wait on last set of jobs
###################################################

echo "$(timestamp): Waiting on final set of jobs to finish..."
wait

echo "$(timestamp): Finished"
