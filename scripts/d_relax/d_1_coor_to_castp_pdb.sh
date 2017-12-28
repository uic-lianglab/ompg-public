#!/bin/bash

# Converts NAMD side chain minimized loop coor to newcast pdb format

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
# Base directory for minimized coordinates
BASE_MIN_DIR=$BASE_PARENT_DIR/min
# Name of subdirectory with minimized coordinate files
MIN_PDB_DIR_NAME=fixed_post

# Path to base directoy for pretzel multiloop samples. These files are in PDB
# format and not NAMD format; therefore, they can serve as a template for
# converting from one format to the other
BASE_LOOP_DIR=$BASE_OUTPUT_DIR/multi_loop
# Make sure path is global
BASE_LOOP_DIR=$(cd "$BASE_LOOP_DIR"; pwd)

# Base output directory for CASTp
BASE_CASTP_OUTPUT_DIR=$BASE_PARENT_DIR/castp
# Path to python script
PY_SCRIPT_PATH=$ROOT_DIR/scripts/coor_to_castp_pdb.py
# Path for storing logs
BASE_LOG_DIR=$DATA_DIR/logs
BASE_CASTP_LOG_DIR=$BASE_LOG_DIR/relax/castp

###################################################
# Completed data sets to ignore
###################################################

# Set of previously completed loops which are ignored by script
# To ignore an entire loop, add element "loop_<#>_<mut>"
# - example, to ignore loop 1 wildtype, add element: "loop_1_wt"
# To ignore only a template within a loop, add element "loop_<#>_<mut>/<template>"
# - example, to ignore loop 1 wildtype at 2iwv template, add element: "loop_1_wt/2iwv"
# To ignore a pH folder, add element "pH<#>/<sim_id>"
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
MAX_NPROC=20

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

cd $BASE_MIN_DIR

# Iterate over subdirectories containing minimized coords
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # http://stackoverflow.com/questions/59838/check-if-a-directory-exists-in-a-shell-script
    # Minimized PDBs exist in MIN_PDB_DIR_NAME subdirectory
    if [ ! -d "$d/$MIN_PDB_DIR_NAME" ]; then
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Full path to minimized coords directory
    min_coords_dir="$d/$MIN_PDB_DIR_NAME"
    min_coords_dir=$(cd "$min_coords_dir"; pwd)
    
    # Full path to CASTp formatted PDBs directory
    castp_pdbs_dir="$BASE_CASTP_OUTPUT_DIR/$d/pdb"
    mkdir -p $castp_pdbs_dir
    castp_pdbs_dir=$(cd "$castp_pdbs_dir"; pwd)
    
    # Make log directory
    log_dir="$BASE_CASTP_LOG_DIR/$d/pdb"
    mkdir -p $log_dir
    
    # Process minimized coordinates
    count=0
    for f in $d/$MIN_PDB_DIR_NAME/*.coor
    do
        # Update count
        ((count++))
        
        # Output PDB
        fid="$(basename $f .coor)"
        pdb_path="$castp_pdbs_dir/$fid.pdb"
        
        # Input coordinates
        coor_path="$min_coords_dir/$fid.coor"
        
        # Determine path to template PDB. Here we use the raw multi-loop file
        # prior to relaxation. We assume the minimized .coor file name follows
        # this format:
        #   <sim_id>.<template_id:2iwv|2iww><misc parameters and timestamp...>.pH<integer>.dup<integer>.min.fixed.coor
        # Example .coor file name:
        #   loop_6_full_neu.2iww.rot1.ofa0.85.ofna0.95.ofm0.75.nds32.nsc16.mxbb0.mxsc15.20160801133526.99_2iww_9.pH5.dup0.min.fixed.coor
        # Here we see that:
        #   sim_id = loop_6_full_neu
        #   template_id = 2iww
        # We assume the corresponding input PDB is named:
        #   <sim_id>.<template_id:2iwv|2iww><misc parameters and timestamp...>.pdb
        # and is located at:
        #   <BASE_LOOP_DIR>/<sim_id>_<template_id:2iww|2iww>/pdb/
        # For the example input file, the corresponding template pdb is named:
        #   loop_6_full_neu.2iww.rot1.ofa0.85.ofna0.95.ofm0.75.nds32.nsc16.mxbb0.mxsc15.20160801133526.99_2iww_9.pdb
        # and is located at:
        #   <BASE_LOOP_DIR>/loop_6_full_neu_2iww/pdb/
        # Hence, to find the corresponding template (do not confuse with x-ray templates 2iwv and 2iwv as our input
        # .coor may be mutated and hence requires the original mutated source file), we are going to use bash specific
        # string manipulations to strip and capture various file extensions. Please see the following references:
        # http://www.tldp.org/LDP/LG/issue18/bash.html
        # http://stackoverflow.com/questions/2664740/extract-file-basename-without-path-and-extension-in-bash/
        
        # Equivalent to basename - remove directory info
        template_path=${coor_path##*/}
        # Strip outermost extension '.coor'
        template_path=${template_path%.*}
        # Strip next outermost extension '.fixed'
        template_path=${template_path%.*}
        # Strip next outermost extension '.min'
        template_path=${template_path%.*}
        # Strip next outermost extension '.dup<integer>'
        template_path=${template_path%.*}
        # Strip next outermost extension '.pH<integer>'
        template_path=${template_path%.*}
        # Append .pdb extension
        template_path="$template_path.pdb"
        # Get simulation identifier
        sim_id=${template_path%%.*}
        
        # Obtain xray input - default is 2iwv
        # http://stackoverflow.com/questions/229551/string-contains-in-bash
        xray="2iwv"
        if [[ $fid =~ .*2iww.* ]]
        then
            xray="2iww"
        fi
        
        # Final template path with directory information
        template_path="$BASE_LOOP_DIR/$sim_id"_"$xray/pdb/$template_path"
        
        # Build command line
        cmd="python $PY_SCRIPT_PATH $pdb_path $coor_path $template_path"
        echo "$count: $(timestamp): $cmd"
        
        log_path="$log_dir/$fid.txt"
        
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
        
    done # end iteration over files in current directory
    
done # end iteration over target folder

popd

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
