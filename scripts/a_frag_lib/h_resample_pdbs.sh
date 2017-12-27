#!/bin/bash 

# Script calls python fragment resampler for each (.rs, .pdb) in corresponding
# directories. The resampler generates a new resampled fragment library based
# on the resampling indices.

###################################################
# Script paths
###################################################

# Directory containing this script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Name of resampling script
RS_SCRIPT_NAME=frag_lib_pdb_resampler.py
# Path to resampling script
RS_SCRIPT_PATH=$SCRIPT_DIR/$RS_SCRIPT_NAME
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Output directories
BASE_OUTPUT_DIR=$DATA_DIR/output
# Path to raw fragment libraries directory
BASE_RAW_FRAG_LIB_INPUT_DIR=$BASE_OUTPUT_DIR/frag_libs/b_cat_libs
# Path to resampled indices directory
BASE_RS_IX_INPUT_DIR=$BASE_OUTPUT_DIR/frag_libs/g_resample_ids
# Path to write resampled fragment libraries
BASE_RS_FRAG_LIB_OUTPUT_DIR=$BASE_OUTPUT_DIR/frag_libs/h_resample_pdbs
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_FRAG_LIB_LOG_DIR=$BASE_LOG_DIR/frag_libs/h_resample_pdbs

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
# Process raw fragment libraries
###################################################

pushd ./

cd $BASE_RAW_FRAG_LIB_INPUT_DIR

# Iterate over subdirectories containing pdb files
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED_LOOPS[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of pdb files in child directory
    has_lib=`ls -1 $d/*.pdb 2>/dev/null | wc -l`
    # If count > 0, then this directory is a frag lib
    # So, skip all folder with no pdb files
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Directory to write cluster results
    outdir="$BASE_RS_FRAG_LIB_OUTPUT_DIR/$d"
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
    for f in $d/*.pdb
    do
        base_name_no_ext="$(basename $f .pdb)"
        echo "Processing raw library for $base_name_no_ext"
        
        # I/O data path arguments resampler
        in_raw_frag_path="$infulld"/"$base_name_no_ext".pdb
        # Path to resampling indices
        in_rs_ix_path="$BASE_RS_IX_INPUT_DIR"/"$d"/"$base_name_no_ext".rs.csv
        # Output path for resampled fragment library
        out_rs_frag_path="$outdir"/"$base_name_no_ext".rs.pdb
    
        logpath="$logdir"/"$base_name_no_ext".rs.pdb.txt
        
        # Build command string
        cmd="python $RS_SCRIPT_PATH"
        cmd="$cmd -rs $in_rs_ix_path"
        cmd="$cmd -ipdb $in_raw_frag_path"
        cmd="$cmd -opdb $out_rs_frag_path"
        
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
            # Sleep for x minutes
            sleep 1m
        done # end spin until queue frees up
        
    done # end iteration over raw fragment libraries in current directory
    
done # end iteration over all raw fragment libraries

popd

###################################################
# Wait on last set of jobs
###################################################

echo "$(timestamp): Waiting on final set of jobs to finish..."
wait

echo "$(timestamp): Finished"
