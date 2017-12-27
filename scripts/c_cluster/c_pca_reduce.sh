#!/bin/bash 

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop/clust
# Input subdirectory
BASE_CLUST_INPUT_DIR=$BASE_OUTPUT_DIR/b_bin_pdbs
# Output subdirectory
BASE_CLUST_OUTPUT_DIR=$BASE_OUTPUT_DIR/c_pca_pdbs
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_CLUST_LOG_DIR=$BASE_LOG_DIR/multi_loop/clust/c_pca_pdbs
# Executables
EXE_DIR=$ROOT_DIR
EXE_NAME=scoby

###################################################
# Scalars
###################################################

# The percentage of variance explained by the retained PCA coordinatins in (0.0, 1.0]
PERC_KEEP=0.95
# If sample observations are in rows (not as efficient), set this to 1, else set to 0
OBS_IN_ROWS=0
# The number of columns each sample occupies
N_DIMS=3

###################################################
# Completed data sets to ignore
###################################################

# Set of previously completed loops which are ignored by script
# To ignore an entire loop, add element "loop_<#>_<mut>"
# - example, to ignore loop 1 wildtype, add element: "loop_1_wt"
# To ignore only a template within a loop, add element "loop_<#>_<mut>/<template>"
# - example, to ignore loop 1 wildtype at 2iwv template, add element: "loop_1_wt/2iwv"
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
MAX_NPROC=8

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

cd $BASE_CLUST_INPUT_DIR

# Iterate over subdirectories containing bin files
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of bin files in child directory
    has_lib=`ls -1 $d/*.bin 2>/dev/null | wc -l`
    # If bin count > 0, then this directory is a lib
    # So, skip all folder with no binaries
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "$(timestamp): Processing $d"
    
    # Directory to write PCA reduced libraries
    outdir="$BASE_CLUST_OUTPUT_DIR/$d"
    echo "$(timestamp): Creating $outdir"
    mkdir -p "$outdir"
    # Capture absolute path to output directory
    outdir=$(cd "$outdir"; pwd)
    
    # Directory to write logs
    logdir="$BASE_CLUST_LOG_DIR/$d"
    mkdir -p $logdir
    
    # Full path to input bin directory
    binfulld=$(cd "$d"; pwd)
    
    # Process each library in directory
    for lib in $d/*.bin
    do
        # Get absolute input and output paths
        binpath="$binfulld/$(basename $lib)"
        pcapath="$outdir"/"$(basename $lib .bin)".pca.bin
        
        logpath="$logdir"/"$(basename $lib .bin)".pca.txt
        
        # Build command string
        cmd="./$EXE_NAME pca_reduce $PERC_KEEP $OBS_IN_ROWS $N_DIMS $binpath $pcapath"
        
        # Switch to exe dir
        pushd ./
        cd $EXE_DIR
        echo $cmd
        # Run command in background
        $cmd &> $logpath &
        # Capture process id of last process ran
        PID=$!
        # Store process id
        queue $PID
        popd
        
        # Check if we need to spin until at least one of the jobs finishes
        while [ $NUM -ge $MAX_NPROC ]; do
            checkqueue
            # Sleep for x seconds
            sleep 10s
        done # end spin until queue frees up
        
    done # end iteration over bin files in current directory
done # end iteration over all child directories containing bin files

popd

###################################################
# Wait on last set of clustering jobs
###################################################

echo "$(timestamp): Waiting on final set of jobs to finish..."
wait

echo "$(timestamp): Finished"

