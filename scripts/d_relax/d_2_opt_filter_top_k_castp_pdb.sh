#!/bin/bash

# Script creates a directory containing only the top-k energy ranked CASTp
# PDBs. This can make the processing time of subsequent steps faster if
# we only need to process a subset of the energy ranks.

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
# Path to base merge directory
BASE_PARENT_DIR=$BASE_OUTPUT_DIR/relax
# Base output directory for CASTp
BASE_CASTP_DIR=$BASE_PARENT_DIR/castp
# Base output directory for captured NAMD energies
BASE_CAPT_DIR=$BASE_PARENT_DIR/capt
# Name of spreadsheet containing energy rankings
ENER_CSV_FNAME="ener.csv"

###################################################
# Additional globals
###################################################

# The max energy ranks to keep in cast PDB folder
TOP_K=5000

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
# Queue jobs
###################################################

pushd ./

cd $BASE_CASTP_DIR

# Iterate over subdirectories containing PDB subdirs
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
    if [ ! -d "$d/pdb" ]; then
        continue
    fi
    
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "$(timestamp): Processing $d"
    
    # Full path to input PDB directory
    pdb_src_dir="$d/pdb"
    pdb_src_dir=$(cd "$pdb_src_dir"; pwd)
    
    # Full path to output PDB directory
    pdb_dst_dir="$d/pdb.$TOP_K"
    mkdir -p $pdb_dst_dir
    pdb_dst_dir=$(cd "$pdb_dst_dir"; pwd)
    
    # Full path to target energy spreadsheet
    ener_csv_dir="$BASE_CAPT_DIR/$d"
    ener_csv_dir=$(cd "$ener_csv_dir"; pwd)
    ener_csv_path="$ener_csv_dir/$ENER_CSV_FNAME"
    
    # Process energy file, extract ranked PDB basenames
    is_header=1
    erank=0
    while read line
    do
        # Skip header row
        if [ $is_header == 1 ]; then
            is_header=0
            continue
        fi
        
        # Update energy rank
        ((erank++))
        
        # Exit loop if we've processed top-k ranked already
        if [ "$erank" -gt "$TOP_K" ]; then
            break
        fi
        
        # Below code splits a CSV tuple into an array
        # http://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
        # http://www.tldp.org/LDP/abs/html/internalvariables.html
        # IFS stands for input field separator
        IFS=',' read -r -a ener <<< "$line"
        
        # Get file name of PDB (first column of spreadsheet)
        pdb_name="${ener[0]}.pdb"
        
        # Input PDB path
        pdb_src="$pdb_src_dir/$pdb_name"
        pdb_dst="$pdb_dst_dir/$pdb_name"
        
        # Create command line
        cmd="cp $pdb_src $pdb_dst"
        echo "$erank: $(timestamp): $cmd"
        
        # Copy to target
        cp $pdb_src $pdb_dst
         
    done < $ener_csv_path
    
done # end iteration over target directories

popd

echo "Finished."
