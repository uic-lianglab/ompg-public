#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes exported <sample name>, <energy rank> lists and copies each
# sample PDB with the energy rank prefixed in the name. The copied PDBs are
# stored in a subdirectory within the same folder containing the input list.
# The name of the output subfolder is <input_list_name.PDBs>.

###########################################
# Imports
###########################################

import errno
import os
import sys

from shutil import copyfile

###########################################
# Default parameters
###########################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Command line arguments are in this order within sys.argv list
ARG_SID_PATH = 1
CMD_LINE_SIZE = 2

# Number of digits in energy rank prefix for destination PDB file names
ZFILL_PAD = 6

# Project root directory
ROOT_DIR =  os.path.join(SCRIPT_DIR, "..", "..")
ROOT_DIR = os.path.abspath(ROOT_DIR)

# Assumed base directories
OUTPUT_DIR = os.path.join(ROOT_DIR, "ompg", "output")
RELAX_DIR = os.path.join(OUTPUT_DIR, "relax")
CAPT_BASE_DIR = os.path.join(RELAX_DIR, "capt")
CASTP_BASE_DIR = os.path.join(RELAX_DIR, "castp")

# Defaults for testing
DEF_SID_PATH = os.path.join(CAPT_BASE_DIR, "pH5", "wt", "regsc",
                            "sids.merge.close.csv")

#############################################################################
# Methods
#############################################################################

def mkdir_p(path):
    """http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python"""
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def to_pdb_path(sid):
    """Parses sample identifier and determines path to PDB file."""
    IX_MUT = 0
    IX_BARREL = 1
    dots = sid.split(".")
    mut_barrel = dots[IX_MUT] + "_" + dots[IX_BARREL]
    pH = "pH" + sid.split("pH")[1].split(".")[0]
    pdb_name = sid + ".pdb"
    return os.path.join(CASTP_BASE_DIR, pH, mut_barrel, "pdb", pdb_name)
    

def extract_samples(sid_path):
    """Copies each sample in sample list file."""
    target_dir = os.path.splitext(sid_path)[0] + ".PDBs"
    print "Creating output PDB directory: " + target_dir
    mkdir_p(target_dir)

    # Column identifiers for sample file
    COL_SID = 0
    COL_ERANK = 1

    # Process each sample identifier
    with open(sid_path) as f:
        # Eat header row
        next(f)
        for line in f:
            rec = line.strip()
            rec = rec.split(",")
            sid = rec[COL_SID]
            erank = rec[COL_ERANK]
            src = to_pdb_path(sid)
            out_pdb_name = erank.zfill(ZFILL_PAD) + "." + sid + ".pdb"
            dst = os.path.join(target_dir, out_pdb_name)
            print erank + ": Copying from: " + src + "\n\tto: " + dst
            copyfile(src, dst)
    
#############################################################################
# Main
#############################################################################


def __main__():

    sid_path = DEF_SID_PATH

    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count == CMD_LINE_SIZE):
        sid_path = sys.argv[ARG_SID_PATH]
        
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "EXTRACT SAMPLES"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "-sid_path: " + sid_path
    
    extract_samples(sid_path)

    print "Extract samples script finished."


if __name__ == '__main__':
    __main__()
