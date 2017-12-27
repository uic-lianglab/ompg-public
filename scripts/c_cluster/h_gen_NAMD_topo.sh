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
# Path to NAMD base directory
NAMD_BASE_DIR=$ROOT_DIR/tools/namd
# Path to NAMD psfgen
PSFGEN_PATH=$NAMD_BASE_DIR/bin/nix/psfgen
# Path to patch topology tcl script
PATCH_TOPOLOGY_TCL_PATH=$NAMD_BASE_DIR/scripts/patch_topology.tcl

# Path to base output directory
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop
BASE_LS_INPUT_DIR=$BASE_OUTPUT_DIR/clust/f_resample_ls

# Path for storing logs
BASE_LOG_DIR=$DATA_DIR/logs
BASE_TOPO_LOG_DIR=$BASE_LOG_DIR/multi_loop/topo

###################################################
# Loop invariant arrays
###################################################

# pH levels to simulate
PH_LEVELS=("5" "7")

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
# Process list files
###################################################

pushd ./

cd $BASE_LS_INPUT_DIR

# Iterate over subdirectories containing ls files
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of list files in child directory
    has_lib=`ls -1 $d/*.ls.txt 2>/dev/null | wc -l`
    # If count > 0, then this directory is a lib
    # So, skip all folder with no lib files
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Process each pH setting
    for pH in "${PH_LEVELS[@]}"
    do
        pH_dir=$BASE_OUTPUT_DIR/topo/pH$pH/$d
        
        out_coords_dir=$pH_dir/pdb
        mkdir -p $out_coords_dir
        
        out_psf_dir=$pH_dir/psf
        mkdir -p $out_psf_dir
        
        pH_log_dir=$BASE_TOPO_LOG_DIR/pH$pH/$d
        mkdir -p $pH_log_dir
        
        # Process each library in directory
        for f in $d/*.ls.txt
        do
            # Create an associative array for tracking
            # http://www.artificialworlds.net/blog/2012/10/17/bash-associative-array-examples/
            # http://stackoverflow.com/questions/13219634/easiest-way-to-check-for-an-index-or-a-key-in-an-array
            unset path_count
            declare -A path_count
            
            # Process list file, extract paths
            linecount=0
            while read line
            do
                # Check if duplicate line
                if [ ${path_count[$line]+_} ]
                then
                    # http://stackoverflow.com/questions/33150161/bash-shell-increment-associative-array-value/33150162
                    (('path_count[$line]'++))
                else
                    path_count[$line]=0
                fi
                
                ((linecount++))
                
                # Check if file is a duplicate
                # Note: am not doing anything special other than creating
                # a unique file name. It's too much hassle to try and save
                # time by not processing duplicates. FWIW, found about 7%
                # duplicates in a 100k resample file of ~590 clusters.
                dup_id=${path_count[$line]}
                # Input PDB
                in_pdb_path=$line
                in_pdb_dir="$(dirname $line)"
                base_name_no_ext="$(basename $line .pdb)"
                # Input pKas
                in_pka_dir="$in_pdb_dir/../propka"
                in_pka_path="$in_pka_dir/$base_name_no_ext".pka
                # Output protonated coordinates
                out_coords_path="$out_coords_dir/$base_name_no_ext.pH$pH.dup$dup_id.pdb"
                # Output topology
                out_psf_path="$out_psf_dir/$base_name_no_ext.pH$pH.dup$dup_id.psf"
                
                # Build command line
                cmd="$PSFGEN_PATH $PATCH_TOPOLOGY_TCL_PATH"
                cmd="$cmd $in_pdb_path $in_pka_path $pH"
                cmd="$cmd $out_coords_path $out_psf_path"
                echo "$linecount: $(timestamp): $cmd"
                
                log_path="$pH_log_dir/$base_name_no_ext.pH$pH.dup$dup_id.log.txt"
                
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
                
            done < $f
            
        done # end iteration over current directory
    
    done # end iteration over pH
done # end iteration over libs   

popd

###################################################
# Wait on last set of jobs
###################################################

echo "Waiting on final set of jobs to finish..."
wait

echo "Finished."
