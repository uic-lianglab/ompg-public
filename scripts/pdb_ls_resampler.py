#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes:
# - .rs.csv - column list of selected pdb ids (0-based indices)
# - .ls.txt - column list which maps a pdb identifier to its file path
# Script will load sampled indices from .rs.csv and then select the
# pdb files at those indices to generate a new pdb set.

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

# Root project directory
DEFAULT_ROOT_DIR = os.path.join(SCRIPT_DIR, "..")

# Base cluster directory
DEFAULT_CLUST_DIR = os.path.join(DEFAULT_ROOT_DIR, "ompg", "output", "multi_loop",
                                 "clust") 

# Default path to directory containing .rs files
DEFAULT_RS_DIR = os.path.join(DEFAULT_CLUST_DIR, "e_resample_ids")

# Default data set to process
DEFAULT_DATASET_NAME = "wt_2iwv"

# Default resample indices file
DEFAULT_RS_FILE_PATH = \
    os.path.join(DEFAULT_RS_DIR, DEFAULT_DATASET_NAME,
                 DEFAULT_DATASET_NAME + "pca.bin.rs.csv")

# Default path to input file containing paths to .pdb files
DEFAULT_IN_LS_PATH = \
    os.path.join(DEFAULT_CLUST_DIR, "a_csv_pdbs", DEFAULT_DATASET_NAME,
                 DEFAULT_DATASET_NAME + ".ls.txt")

# Default path to output file containing cluster sampled paths to .pdb files
DEFAULT_OUT_LS_PATH = \
    os.path.join(DEFAULT_CLUST_DIR, "f_resample_ls", DEFAULT_DATASET_NAME,
                 DEFAULT_DATASET_NAME + ".rs.ls.txt")

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

# Parses single column of values
# @param path - path to list file
# @return list of strings, each string contains a single line from list file
def load_ls(ls_path):
    ls = []
    with open(ls_path) as f:
        ls = f.readlines()
    return ls

# Utility for loading resample indices and list values
# @param rs_path - path to resample indices file
# @param ls_path - path to list (ls) file
# @return tuple - first element is resample indices (list of integers),
#   second element is list file (list of strings)
def load_data(rs_path, ls_path):
    print "Loading rs indices: " + rs_path
    rs = load_rs_indices(rs_path)
    print "Loaded " + str(len(rs)) + " rs indices"
    print "Loading list file: " + ls_path
    ls = load_ls(ls_path)
    print "Loaded " + str(len(ls)) + " list values"
    return rs, ls

# Selects the indices specified by rs from ls
# @param rs - list of integer indices into param ls
# @param ls - list of strings
# @return list of strings, the resampled values as specified by indices within rs
def resample_ls(rs, ls):
    print "Resampling"
    rs_result = []
    for i in rs:
        rs_result.append(ls[i])
    return rs_result

# Writes resampled values to disk
# @param rs_result - list of resampled strings
# @param path - output path to write resampled values
def write_resampled_ls(rs_result, path):
    print "Writing " + str(len(rs_result)) + " resampled values to: " + path
    # Write fragments to disk
    with open(path, 'w') as f:
        f.writelines(rs_result)

###########################################
# Main
###########################################

# Main script entry point
def __main__():
    print "======================= File List Resampler ======================="
    # Initialize command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-rs', '--path_to_in_rs_file',
        default=DEFAULT_RS_FILE_PATH,
        help='Path to input resampling indices file')
    parser.add_argument('-ils', '--path_to_in_ls_file',
        default=DEFAULT_IN_LS_PATH,
        help='Path to input value list file to resample from')
    parser.add_argument('-ols', '--path_to_out_ls_file',
        default=DEFAULT_OUT_LS_PATH,
        help='Path to output resampled values from input list file')

    # Parse command line
    args = parser.parse_args()
    
    # Print command line
    print '\t-rs = ' + args.path_to_in_rs_file
    print '\t-ils = ' + args.path_to_in_ls_file
    print '\t-ols = ' + args.path_to_out_ls_file

    # Load data
    rs, ls = load_data(args.path_to_in_rs_file, args.path_to_in_ls_file)

    # Resample
    rs_result = resample_ls(rs, ls)

    # Write re-sampled data
    write_resampled_ls(rs_result, args.path_to_out_ls_file)

    print "Resampling finished."

# Run main if we are the active script
if __name__ == '__main__':
    __main__()
