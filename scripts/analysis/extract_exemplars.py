#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes:
# - .ex.csv - column list of examplar indices (0-based indices)
# - .s2c.csv: file containing sample to cluster identifier mappings
# - .ls.txt - column list which maps an identifier (0-based index) to its
#             file path
# Script will copy all [unrelaxed] exemplars to target directory

###########################################
# Imports
###########################################

# For directory creation error handling
import errno

# For parsing user supplied arguments
import argparse

# For parsing list of files in a directory
import os

from shutil import copyfile

###########################################
# Default parameters
###########################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Number of digits in cluster id prefix for destination PDB file names
ZFILL_PAD = 3

###########################################
# Directories for default data paths

# Default data set to process
DEF_MUT_ID="wt"
DEF_BARREL_ID="2iwv"
DEF_MUT_BARREL_ID=DEF_MUT_ID + "_" + DEF_BARREL_ID
DEF_PH_ID="pH5"

# Root project directory
DEF_ROOT_DIR = os.path.join(SCRIPT_DIR, "..", "..")

# Base output directory
DEF_OUTPUT_DIR = os.path.join(DEF_ROOT_DIR, "ompg", "output")

# Base cluster directory
DEF_CLUST_DIR = os.path.join(DEF_OUTPUT_DIR, "multi_loop", "clust")

# Directory containing energy relaxed data sets
DEF_RELAX_DIR = os.path.join(DEF_OUTPUT_DIR, "relax")

# Directory containing 0-run MD runs
DEF_CAPT_BASE_DIR = os.path.join(DEF_RELAX_DIR, "capt")

###########################################
# Default input arguments

# Default exemplar file path
DEFAULT_EX_PATH = os.path.join(DEF_CLUST_DIR, "d_clust_pdbs",
                               DEF_MUT_BARREL_ID,
                               DEF_MUT_BARREL_ID + ".pca.bin.ex.csv")

# Mapping from cluster identifier (0-based) to original file name
DEFAULT_LS_PATH = os.path.join(DEF_CLUST_DIR, "a_csv_pdbs",
                               DEF_MUT_BARREL_ID,
                               DEF_MUT_BARREL_ID + ".ls.txt")

###########################################
# Default output arguments

# Directory to copy exemplar PDBs to
DEF_OUT_EX_PDB_DIR = os.path.join(DEF_CAPT_BASE_DIR, DEF_PH_ID, DEF_MUT_ID,
                                  "regsc", "ex.pdb." + DEF_BARREL_ID)

###########################################
# Utility functions
###########################################

def mkdir_p(path):
    """http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python"""
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# @param ex_path - path to exemplar identifiers
# @return list of exemplar identifiers
def load_ex_ids(ex_path):
    """Reads list of exemplar identifiers."""
    # Process exemplar ids
    print "Loading " + ex_path
    exs = []
    with open(ex_path) as f_ex:
        for line in f_ex:
            int_id = int(line)
            exs.append(int_id)
    return exs
    
# Parses single column of values
# @param path - path to list file
# @return list of strings, each string contains a single line from list file
def load_ls(ls_path):
    """Read list of sample file names."""
    print "Loading " + ls_path
    ls = []
    with open(ls_path) as f:
        ls = f.readlines()
    # Strip newlines
    ls = map(lambda s: s.strip(), ls)
    return ls

###########################################
# Extract exemplars
###########################################

def extract_exemplars(ex_path, ls_path, o_dir):
    """Copy exemplars to separate folder."""
    # Load data    
    ex = load_ex_ids(ex_path)
    ls = load_ls(ls_path)
    # Create output folder
    mkdir_p(o_dir)
    
    # cid: cluster identifier
    # fid: index into ls (map to file name)
    for cid, fid in enumerate(ex):
        src = ls[fid]
        out_pdb_name = str(cid).zfill(ZFILL_PAD) + ".pdb"
        dst = os.path.join(o_dir, out_pdb_name)
        print str(cid) + ": Copying from: " + src + "\n\tto: " + dst
        copyfile(src, dst)

###########################################
# Main
###########################################

# Main script entry point
def __main__():
    print "======================= Examplar Extractor ======================="
    # Initialize command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-ex', '--path_to_exemplar_ids',
        default=DEFAULT_EX_PATH,
        help='Path to exemplar identifiers file')
    parser.add_argument('-ls', '--path_to_ls_file',
        default=DEFAULT_LS_PATH,
        help='Path to input identifier to file path mapping')
    parser.add_argument('-od', '--path_to_out_dir',
        default=DEF_OUT_EX_PDB_DIR,
        help='Path to copy exemplar PDBs')

    # Parse command line
    args = parser.parse_args()

    # Print command line
    print '\t-ex = ' + args.path_to_exemplar_ids
    print '\t-ls = ' + args.path_to_ls_file
    print '\t-od = ' + args.path_to_out_dir

    # Load data
    extract_exemplars(args.path_to_exemplar_ids,
                      args.path_to_ls_file,
                      args.path_to_out_dir)

    print "Exemplar extraction finished."

# Run main if we are the active script
if __name__ == '__main__':
    __main__()
