#!/bin/bash 

# Scripts moves all fragment libraries for output generating
# folders to final location expected by DiSGro.
#
# WARNING: ALL FRAGMENT LIBRARIES MUST HAVE UNIQUE FILENAMES!
#   IF NOT, NAME CLASHES WILL RESULT IN OVERWRITING.

###################################################
# Script paths
###################################################

# Directory containing this script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Path with resampled fragment libraries
BASE_RS_FRAG_LIB_INPUT_DIR=$DATA_DIR/output/frag_libs/h_resample_pdbs
# All fragment libraries are moved to this folder
# (DiSGro expects all libraries to reside in same folder)
FINAL_FRAG_LIB_DIR=$DATA_DIR/frag_libs

# Make sure output folder exists
mkdir -p $FINAL_FRAG_LIB_DIR

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
# Process fragment libraries
###################################################

pushd ./

cd $BASE_RS_FRAG_LIB_INPUT_DIR

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
        
    # Process each library in directory
    for f in $d/*.pdb
    do
        frag_lib_name="$(basename $f)"
        outpath="$FINAL_FRAG_LIB_DIR"/"$frag_lib_name"
        echo "$(timestamp): Copying $f to $outpath"
        
        cp $f $outpath
        
    done # end iteration over fragment libraries in current directory
    
done # end iteration over all fragment libraries

popd

echo "$(timestamp): Finished."
