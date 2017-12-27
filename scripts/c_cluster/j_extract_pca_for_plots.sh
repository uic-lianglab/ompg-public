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
BASE_CLUST_OUTPUT_DIR=$BASE_OUTPUT_DIR/g_pca_for_plots
# Log directories
BASE_LOG_DIR=$DATA_DIR/logs
BASE_CLUST_LOG_DIR=$BASE_LOG_DIR/multi_loop/clust/g_pca_for_plots
# Executables
EXE_DIR=$ROOT_DIR
EXE_NAME=scoby

###################################################
# Scalars
###################################################

# The number of PCA dimensions to keep
# Each sample is assumed to occupy a single column
N_KEEP_ROWS="2"
# Number of samples to keep, if <= 0, all samples kept
N_KEEP_COLS="-1"

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
# Process datasets
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
        
        #########################################
        # CROP FULL PCA MATRIX
        #
        
        in_matrix_path="$binfulld/$(basename $lib)"
        
        # Output cropped bin path
        out_pca_plot_bin_path="$outdir"/"$(basename $lib .pca.bin)".pca.plot.bin
        
        logpath="$logdir"/"$(basename $lib .pca.bin)".crop.txt
        
        # Build command string
        cmd="./$EXE_NAME crop_fast_matrix_bin $N_KEEP_ROWS $N_KEEP_COLS $in_matrix_path $out_pca_plot_bin_path"
        
        # Switch to exe dir
        pushd ./
        cd $EXE_DIR
        echo $cmd
        # Run command
        $cmd &> $logpath
        popd
        
        #########################################
        # CONVERT TO CSV
        #
        
        in_matrix_path="$out_pca_plot_bin_path"
        
        # Output cropped csv path
        out_pca_plot_csv_path="$outdir"/"$(basename $lib .pca.bin)".pca.plot.csv
        
        logpath="$logdir"/"$(basename $lib .pca.bin)".tocsv.txt
        
        cmd="./$EXE_NAME convert_fast_matrix_bin_to_csv $in_matrix_path $out_pca_plot_csv_path"
        
        # Switch to exe dir
        pushd ./
        cd $EXE_DIR
        echo $cmd
        # Run command
        $cmd &> $logpath
        popd
        
    done # end iteration over bin files in current directory
done # end iteration over all child directories containing bin files

popd

echo "$(timestamp): Finished"
