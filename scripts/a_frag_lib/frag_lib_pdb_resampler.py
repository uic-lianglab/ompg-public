#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes:
# - .rs.csv - column list of selected fragment ids (0-based indices)
# - .lib.pdb - a fragment library containing pdb atoms for each fragment
# Script will load sampled indices from .rs.csv and then select the
# fragments at those indices to generate a new .lib.pdb library

###########################################
# Imports
###########################################

# For parsing user supplied arguments
import argparse

# For parsing list of files in a directory
import os

###########################################
# Default parameters
###########################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Default path to directory containing .rs files
DEFAULT_RS_DIR = \
    os.path.join(SCRIPT_DIR, "../ompg/output/frag_libs/g_resample_ids/loop_6_wt/2iwv")

# Default path to directory containing .pdb files
DEFAULT_PDB_DIR = \
    os.path.join(SCRIPT_DIR, "../ompg/output/frag_libs/b_cat_libs/loop_6_wt/2iwv")

# Default data set to process
DEFAULT_DATASET_NAME = "loop_6_wt.2iwv.frag.lib"

# Default resample indices file
DEFAULT_RS_FILE = \
    os.path.join(DEFAULT_RS_DIR, DEFAULT_DATASET_NAME + ".rs.csv")
    
# Default raw pdb atom coordinates (fragment library) file
DEFAULT_IN_PDB_FILE = \
    os.path.join(DEFAULT_PDB_DIR, DEFAULT_DATASET_NAME + ".pdb")

# Default out resampled file path
DEFAULT_OUT_PDB_FILE = \
    os.path.join(DEFAULT_RS_DIR, DEFAULT_DATASET_NAME + ".rs.pdb")

###########################################
# Utility functions
###########################################

# Parses resample indices file
# @param rs_path - path to resample indices file
# @return list of integer indices denoting which fragments to keep
def load_rs_indices(rs_path):
    rs = []
    with open(rs_path) as f:
        for line in f:
            i = int(line.strip())
            rs.append(i)
    return rs

# Parses pdb fragment libraries
# @param pdb_path - path to pdb atom coordinates of fragments (single file)
# @return list of strings, each string contains atom coordinates of an entire fragment
def load_frag_lib(pdb_path):
    frag_lib = []
    with open(pdb_path) as f:
        frag = ""
        for line in f:
            if "ENDMDL" in line:
                # strip any newline character
                frag += line.strip()
                frag_lib.append(frag)
                frag = ""
            else:
                # keep newline character
                frag += line
    return frag_lib

# Utility for loading resample indices and fragment library
# @param rs_path - path to resample indices file
# @param pdb_path - path to pdb atom coordinates of fragments (single file)
# @return tuple - first element is resample indices (list of integers),
#   second element is fragment library (list of strings)
def load_data(rs_path, pdb_path):
    print "Loading rs indices: " + rs_path
    rs = load_rs_indices(rs_path)
    print "Loaded " + str(len(rs)) + " rs indices"
    print "Loading fragment library: " + pdb_path
    frag_lib = load_frag_lib(pdb_path)
    print "Loaded " + str(len(frag_lib)) + " fragments"
    return rs, frag_lib

# Selects the indices specified by rs from param frag_lib
# @param rs - list of integer indices into param frag_lib
# @param frag_lib - list of strings, each string contains atom coordinates of an entire fragment
# @return list of strings, the resampled fragment library as specified by indices within rs
def resample_frag_lib(rs, frag_lib):
    print "Resampling fragment library"
    rs_frag_lib = []
    for i in rs:
        rs_frag_lib.append(frag_lib[i])
    return rs_frag_lib

# Writes fragment library to disk
# @param frag_lib - list of strings, each string contains atom coordinates of an entire fragment
# @param pdb_path - output pdb path to write fragment library's atom coordinates
def write_frag_lib(frag_lib, pdb_path):
    print "Writing " + str(len(frag_lib)) + " fragments to: " + pdb_path
    # Write fragments to disk
    with open(pdb_path, 'w') as f:
        for idx, frag in enumerate(frag_lib):
            if ((idx+1) < len(frag_lib)):
                f.write(frag + '\n')
            else:
                f.write(frag)

###########################################
# Main
###########################################

# Main script entry point
def __main__():
    print "======================= Frag Lib Resampler ======================="

    # Initialize command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-rs', '--path_to_in_rs_file',
        default=DEFAULT_RS_FILE,
        help='Path to input resampling indices file')
    parser.add_argument('-ipdb', '--path_to_in_pdb_file',
        default=DEFAULT_IN_PDB_FILE,
        help='Path to input fragment library pdb file')
    parser.add_argument('-opdb', '--path_to_out_pdb_file',
        default=DEFAULT_OUT_PDB_FILE,
        help='Path to output fragment')

    # Parse command line
    args = parser.parse_args()
    
    # Print command line
    print '\t-rs = ' + args.path_to_in_rs_file
    print '\t-ipdb = ' + args.path_to_in_pdb_file
    print '\t-opdb = ' + args.path_to_out_pdb_file

    # Load data
    rs, frag_lib = load_data(args.path_to_in_rs_file, args.path_to_in_pdb_file)

    # Resample
    rs_frag_lib = resample_frag_lib(rs, frag_lib)

    # Write re-sampled data
    write_frag_lib(rs_frag_lib, args.path_to_out_pdb_file)

    print "Resampling finished."

# Run main if we are the active script
if __name__ == '__main__':
    __main__()
