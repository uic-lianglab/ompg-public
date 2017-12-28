#!/bin/bash 

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Output directories
BASE_OUTPUT_DIR=$DATA_DIR/output/multi_loop/clust
# Input subdirectory
BASE_CLUST_INPUT_DIR=$BASE_OUTPUT_DIR/a_csv_pdbs
# Output subdirectory
BASE_CLUST_OUTPUT_DIR=$BASE_OUTPUT_DIR/b_bin_pdbs
# Executable
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
# Convert files
###################################################

pushd ./

cd $BASE_CLUST_INPUT_DIR

# Iterate over subdirectories containing CSV files
for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of CSV files in child directory
    has_lib=`ls -1 $d/*.csv 2>/dev/null | wc -l`
    # If CSV count > 0, then this directory has a library
    # So, skip all folder with no CSVs
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "$(timestamp): Processing $d"
    
    outdir="$BASE_CLUST_OUTPUT_DIR/$d"
    
    echo "$(timestamp): Creating $outdir"
    mkdir -p "$outdir"
    
    # Full path to input csv directory
    csvfulld=$(cd "$d"; pwd)
    
    # Process each library in directory
    for lib in $d/*.csv
    do
        # Get absolute path to this library
        csvpath="$csvfulld/$(basename $lib)"
        
        # Run pretzel converter
        pushd ./
        cd $EXE_DIR
        binpath="$outdir"/"$(basename $lib .csv)".bin
        echo "$(timestamp): Converting $csvpath to $binpath"
        # 0 -> no header row
        ./$EXE_NAME convert_matrix_csv_to_bin "$csvpath" "$binpath" 0
        popd
    done # end iteration over CSV files in current directory
done # end iteration over all child directories containing CSVs

popd

echo "$(timestamp): Finished cluster prep conversion from CSV to BIN"
