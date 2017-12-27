#!/bin/bash 

# Script merges/concatenates all raw fragment libraries into single file

###################################################
# Script paths
###################################################

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
# Root project directory
ROOT_DIR=$SCRIPT_DIR/../..
# Base data directory
DATA_DIR=$ROOT_DIR/ompg
# Output directories
BASE_PDB_OUTPUT_DIR=$DATA_DIR/output
BASE_FRAG_LIB_PDB_INPUT_DIR=$BASE_PDB_OUTPUT_DIR/frag_libs/a_raw_libs
BASE_FRAG_LIB_PDB_OUTPUT_DIR=$BASE_PDB_OUTPUT_DIR/frag_libs/b_cat_libs

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
# Concatenate files
###################################################

pushd ./

cd $BASE_FRAG_LIB_PDB_INPUT_DIR

for d in $(find . -type d)
do
    # Skip this folder if already processed
    # http://unix.stackexchange.com/questions/131568/is-a-directory-error-when-trying-to-pass-directory-name-into-function
    is_complete="$(array_contains_substring_of "$d" "${COMPLETED_LOOPS[@]}")"
    if [ $is_complete == 1 ]; then
        echo "Skipping $d as previously completed."
        continue
    fi
    
    # Count number of PDB files in child directory
    has_lib=`ls -1 $d/*.pdb 2>/dev/null | wc -l`
    # If PDB count > 0, then this directory is a frag lib
    # So, skip all folder with no PDBs
    if [ $has_lib == 0 ]; then 
        continue
    fi
    
    echo "$(timestamp): Processing $d"
    
    outdir="$BASE_FRAG_LIB_PDB_OUTPUT_DIR/$d"
    echo "$(timestamp): Creating $outdir"
    mkdir -p "$outdir"
    
    # Determine name of fragment library - derive it
    # from its current child directory paths
    
    # Strip "./" from the directory listing
    libname=$d
    # http://stackoverflow.com/questions/2172352/in-bash-how-can-i-check-if-a-string-begins-with-some-value
    if [[ "$libname" == \.\/* ]]; then
        # http://stackoverflow.com/questions/11469989/how-can-i-strip-first-x-characters-from-string-in-shellscript-using-sed
        libname="${libname:2}"
    fi
    
    # Replace any "/" character with a "."
    # http://stackoverflow.com/questions/5928156/replace-a-space-with-a-period-in-bash
    # http://tldp.org/LDP/abs/html/string-manipulation.html
    libname=${libname//\//\.}
    
    # Add 'frag.lib.pdb' suffix
    libname="$libname.frag.lib.pdb"
    libpath="$outdir/$libname"
    
    echo "$(timestamp): Appending fragment library: $libpath"
    
    # Append files
    cat $d/*.pdb >> $libpath
    
done

popd

echo "$(timestamp): Finished fragment library concatenation"
