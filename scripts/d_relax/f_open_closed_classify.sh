#!/bin/bash

# Script classifies PDB as being in "open" or "close" state based on
# previously generated CASTp data

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
# Path to base parent directory
BASE_PARENT_DIR=$BASE_OUTPUT_DIR/relax
# Base output directory for CASTp
BASE_CASTP_OUTPUT_DIR=$BASE_PARENT_DIR/castp
# Path to python script
PY_SCRIPT_PATH=$ROOT_DIR/scripts/pdb_open_closed_classifier_template_aware.py
# Path for storing logs
BASE_LOG_DIR=$DATA_DIR/logs
BASE_CASTP_LOG_DIR=$BASE_LOG_DIR/relax/castp

###################################################
# Additional globals
###################################################

# Note: 1.4 angstroms is default probe_radius (for H20)
# Potassium = 2.75 angstroms (wolfram alpha)
# Chlorine = 1.75 angstroms (wolfram alpha)
# According to Chen et al, PNAS 2008:
# widest point in 2iww (open) is 8 angstroms
# widest point in 2iwv (closed) is 1.4 angstroms
PROBE_RADIUS=("2.75")

# Minimum number of mouth atoms above the bilayer
# for a pocket to classified as spanning the outer bilayer
MIN_MOUTH_ATOMS_ABOVE_BILAYER="2"

# Minimum number of mouth atoms below the bilayer
# for a pocket to be classified as spanning the inner bilayer
MIN_MOUTH_ATOMS_BELOW_BILAYER="2"

# For a pocket with mouths on both sides of the bilayer,
# minimum solvent accessible surface area
# necessary for the protein to be considered 'open'
MIN_SA_AREA="50.0" # small value for now

# Only process the top-k energy ranked PDBs
# - Assumes existence of "pdb.<TOP_K>" subfolder which has been pre-filtered
# If less than or equal to 0, then:
# - Assumes existence of "pdb" subfolder with no energy filtering
TOP_K=-1

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
# Derived variables
###################################################

# Determine name of PDB subdirectory
PDB_DNAME="pdb"
if [ "$TOP_K" -gt  0 ]; then
    PDB_DNAME="pdb.$TOP_K"
fi

###################################################
# Queue jobs
###################################################

pushd ./

cd $BASE_CASTP_OUTPUT_DIR

# Iterate over directories containing CASTp subdirs
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Check for subfolder with CASTp formatted PDBs
    # http://stackoverflow.com/questions/59838/check-if-a-directory-exists-in-a-shell-script
    if [ ! -d "$d/$PDB_DNAME" ]; then
        continue
    fi
    
    # Full path to PDB directory
    pdb_dir="$d/$PDB_DNAME"
    pdb_dir=$(cd "$pdb_dir"; pwd)
    
    # Iterate over probe radius
    for radius in "${PROBE_RADIUS[@]}"
    do
        # Check if probe radius subdirectory is missing
        # http://stackoverflow.com/questions/59838/check-if-a-directory-exists-in-a-shell-script
        if [ ! -d "$d/$radius" ]; then
            continue
        fi
        
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo "$(timestamp): Processing probe radius: $d/$radius"
        
        # Parent directory containing CASTp data
        parent_dir=$(cd "$d"; pwd)
        # Directory containing CASTp topology data for probe radius
        castp_topo_dir="$parent_dir/$radius"
        # Path to output open/closed classification
        oc_path="$parent_dir/oc$radius.csv"
        
        # Build command line
        cmd="python $PY_SCRIPT_PATH -omc_ic $castp_topo_dir -omc_ip $pdb_dir -omc_op $oc_path"
        cmd="$cmd -omc_min_o $MIN_MOUTH_ATOMS_ABOVE_BILAYER"
        cmd="$cmd -omc_min_i $MIN_MOUTH_ATOMS_BELOW_BILAYER"
        cmd="$cmd -omc_min_sa $MIN_SA_AREA"
        
        echo "$(timestamp): $cmd"
        
        log_dir="$BASE_CASTP_LOG_DIR/$d"
        mkdir -p $log_dir
        log_path="$log_dir/oc$radius.txt"
        
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
        
    done # end iteration over probe directories
done # end iteration over target directories

popd

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
