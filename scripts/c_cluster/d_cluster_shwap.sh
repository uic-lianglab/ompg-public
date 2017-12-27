#!/bin/bash 

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop/clust
# Input subdirectory
BASE_CLUST_INPUT_DIR=$BASE_OUTPUT_DIR/c_pca_pdbs
# Output subdirectory
BASE_CLUST_OUTPUT_DIR=$BASE_OUTPUT_DIR/d_clust_pdbs
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_CLUST_LOG_DIR=$BASE_LOG_DIR/multi_loop/clust/d_clust_pdbs
# Executables
EXE_DIR=$ROOT_DIR
EXE_NAME=scoby

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
MAX_NPROC=10

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
# Loop invariant arguments to clustering tool
###################################################

# Number of columns defining a single sample
N_DIMS=1
# The first INIT_SIZE elements are clustered using affinity propagation
# to find the initial exemplar and cluster assignments
INIT_SIZE=10000
# The type of heuristic to use for attempting to associate an unlabeled sample
# to an existing exemplar
# 0 -> nascent universal average distance to exemplar 
# 1 -> nascent universal median average distance to exemplar
EX_ASSOC_HEUR=1
# The association heuristic cutoff distance is scaled by this factor
# -> values less than 1.0 make it harder to associate an unlabeled sample to an exemplar
# -> values greater than 1.0 make it easier to associate an unlabeled sample to an exemplar
EX_ASSOC_HEUR_SCALE=0.4
# Threshold for triggering a weighted, local affinity propagation on samples
# which failed to be associated to an exemplar
MAX_UNLABELED_THRESH=2500
# Threshold for triggering a weighted, hierarchical affinity propagation run on only
# the exemplars. This is to attempt to merge exemplars and reduce the total count in order
# for the local affinity propagation runs to be performant.
MAX_EXEMPLARS_THRESH=10000
# Each exemplar has a window length defined as the number of samples that have been processed since
# the last time a sample was associated to that exemplar. If the window length exceeds the decay value
# then the exemplar is removed. Furthermore, the multiplicity weights associated to that exemplar are
# reduced as a function of this decay length.
DECAY_WINDOW=15000
# Specifies which quantile similarity to assign as the no-prior preference
Q=0.5
# Maximum number of iterations for each affinity propagation run
MAXITS=1000
# A single run of affinity propagation terminates if the examplars have not changed for convits
# iterations
CONVITS=100
# Damping factor: should be a value in the range [0.5, 1); higher values correspond
# to heavy damping which may be needed if oscillations occur
LAM=0.90
# If false, adds a small amount of noise to similarities to prevent degenerate cases
NONOISE=0
# Quantile preferences are decayed by this amount if any single affinity propagation
# trial fails to converge. This will usually allow convergence on the subsequent run.
P_CONV_DECAY_FACT=1.05
# Maximum number of decay steps to take when stock affinity propagation fails to
# converge. The program will fail if convergence is not reached by this number of steps.
P_CONV_MAX_ATTEMPTS=5
# Maximum number of times to swizzle the input samples in an attempt to
# obtain an ensemble that converges during the initial affinity propagation run.
INIT_SWIZ_MAX_ATTEMPTS=5

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
    # So, skip all folder with no CSVs
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "$(timestamp): Processing $d"
    
    # Directory to write cluster results
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
        in_matrix_path="$binfulld/$(basename $lib)"
        out_samples_to_clusters_path="$outdir"/"$(basename $lib .coords.pca.bin)".s2c.csv
        out_clusters_to_samples_path="$outdir"/"$(basename $lib .coords.pca.bin)".c2s.csv
        out_exemplars_path="$outdir"/"$(basename $lib .coords.pca.bin)".ex.csv
        
        logpath="$logdir"/"$(basename $lib .coords.bin)".shwap.txt
        
        # Build command string
        cmd="./$EXE_NAME shwap $in_matrix_path $N_DIMS"
        cmd="$cmd $out_samples_to_clusters_path $out_clusters_to_samples_path"
        cmd="$cmd $out_exemplars_path $INIT_SIZE $EX_ASSOC_HEUR $EX_ASSOC_HEUR_SCALE"
        cmd="$cmd $MAX_UNLABELED_THRESH $MAX_EXEMPLARS_THRESH $DECAY_WINDOW $Q $MAXITS"
        cmd="$cmd $CONVITS $LAM $NONOISE $P_CONV_DECAY_FACT $P_CONV_MAX_ATTEMPTS $INIT_SWIZ_MAX_ATTEMPTS"
        
        # Switch to exe dir
        pushd ./
        cd $EXE_DIR
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
        
    done # end iteration over bin files in current directory
done # end iteration over all child directories containing bin files

popd

###################################################
# Wait on last set of clustering jobs
###################################################

# Wait on last set of jobs
echo "$(timestamp): Waiting on final set of jobs to finish..."
wait

echo "$(timestamp): Finished"
