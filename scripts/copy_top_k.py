#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script will extract top-k open and top-k closed PDBs.

Usage:

python <path_to_script>
"""


#############################################################################
# Imports
#############################################################################

import math
import os
import sys
from shutil import copyfile
from collections import namedtuple


#############################################################################
# Global paths
#############################################################################

# Path to script directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Path to base data directory
BASE_DATA_DIR = os.path.join(SCRIPT_DIR, "..", "ompg", "output", "relax_c_beta")


#############################################################################
# Configuration
#############################################################################

# The individual mutations
MUTATION_TYPES = ["wt"]

# The pH levels simulated
PH_LEVELS = ["5", "7"]

# Probe size used in open/close classification
PROBE_RADII = [2.75]

# Store the current configuration
OcConfig = namedtuple('OcConfig', "mut, pH, probe")

# Leave as empty string if all .coor are converted to PDB
MAX_OC_PDB_SUFFIX = ".5000"

# The top 'k' to copy
K = 10

# Make sure file names are padded with zeros so that they are all the same
# length and therefore sorted properly in file browsers
Z_FILL = 1 + int(math.log(K, 10))

#############################################################################
# Offsets
#############################################################################

# Indexed entries within a score row as stored on disk
IX_SCORE_FILE_COORD_NAME = 0
IX_SCORE_FILE_PH = 1
IX_SCORE_FILE_POTENTIAL = 2

# Indexed entries within a score tuple as stored in memory
IX_SCORE_COORD_NAME = 0
IX_SCORE_POTENTIAL = 1
IX_SCORE_LINE = 2

# Indexed entries within an open/closed record (disk or memory)
IX_OC_COORD_NAME = 0
IX_OC_STATUS = 1

# Status codes for open and closed configurations
OC_STATUS_CLOSE = 0
OC_STATUS_OPEN = 1


#############################################################################
# Utilities
#############################################################################


# http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Will make intermediate directories
def make_directory(dir_path):
    """Utility for making a directory if not existing."""
    abs_dir_path = os.path.abspath(dir_path)
    if not os.path.exists(abs_dir_path):
        os.makedirs(abs_dir_path)


def get_oc_path(cfg):
    """Obtain path to open/close information."""
    return os.path.join(
        BASE_DATA_DIR,
        "castp",
        "pH" + str(cfg.pH),
        str(cfg.mut),
        "oc" + str(cfg.probe) + ".csv")


def get_score_dir(cfg):
    """Obtain directory to energy score rankings."""
    return os.path.join(
        BASE_DATA_DIR,
        "capt",
        "pH" + str(cfg.pH),
        str(cfg.mut))


def get_score_path(cfg):
    """Obtain path to energy score rankings."""
    return os.path.join(
        get_score_dir(cfg),
        "ener.csv")


def get_in_pdb_path(cfg, coord_id):
    """Obtain path to pdb file."""
    return os.path.join(
        BASE_DATA_DIR,
        "castp",
        "pH" + str(cfg.pH),
        str(cfg.mut),
        "pdb" + MAX_OC_PDB_SUFFIX,
        coord_id + ".pdb")


def get_top_k_dirs(cfg):
    """Obtain folders to write each category of top k file"""
    base_top_k_dir = os.path.join(get_score_dir(cfg), "top_k" + "." + str(cfg.probe))
    global_dir = os.path.join(base_top_k_dir, "global")
    open_dir = os.path.join(base_top_k_dir, "open")
    close_dir = os.path.join(base_top_k_dir, "close")
    return {"global":global_dir, "open":open_dir, "close":close_dir}
        

#############################################################################
# Parsing
#############################################################################


# @return list of 3-tuples
#   tuple[0] -> base pdb name (string)
#   tuple[1] -> potential (float)
#   tuple[2] -> line in score.csv (string)
def read_score(path):
    """Read energy profile."""
    out = []
    with open(path) as f:
        # Skip header row
        f.next()
        for line in f:
            rec = line.split(',')
            out.append((rec[IX_SCORE_FILE_COORD_NAME],
                        float(rec[IX_SCORE_FILE_POTENTIAL]),
                        line))
    return out


# @return dictionary from coord identifier to open/close status
#   0 -> closed, 1 -> open
def read_oc(path):
    """Read open/close status."""
    out = {}
    with open(path) as f:
        for line in f:
            rec = line.split(',')
            coord_id = rec[IX_OC_COORD_NAME]
            out[coord_id] = int(rec[IX_OC_STATUS])
    return out


#############################################################################
# Global data sets
#############################################################################

# Lazy-loaded on demand mapping from (mutant, pH) tuple to energy rank records
SCORE_MAP = {}

# Lazy-loaded on demand mapping from (mutant, pH) tuple to oc status dict
OC_MAP = {}

# Maps and lists are passed by reference
def lazy_load(map_obj, key, loader, pather, cfg):
    """Lazy load a key if not present."""
    if key not in map_obj:
        map_obj[key] = loader(pather(cfg))
    return map_obj[key]


def get_score(cfg):
    """Obtain energy ranks for (mutation, pH)."""
    key = (cfg.mut, cfg.pH)
    return lazy_load(SCORE_MAP, key, read_score, get_score_path, cfg)


def get_oc(cfg):
    """Obtain open/close status for (mutation, pH, probe)."""
    key = (cfg.mut, cfg.pH, cfg.probe)
    return lazy_load(OC_MAP, key, read_oc, get_oc_path, cfg)


def get_top_k(cfg):
    """Obtain coord ids for top k global, open, and close."""
    score = get_score(cfg)
    oc = get_oc(cfg)
    top_k_open = []
    top_k_close = []
    top_k_global = []
    for score_tup in score:
        coord_id = score_tup[IX_SCORE_COORD_NAME]
        oc_status = oc[coord_id]
        if (len(top_k_global) < K):
            top_k_global.append(coord_id)
        if (oc_status == OC_STATUS_OPEN):
            if (len(top_k_open) < K):
                top_k_open.append(coord_id)
        else:
            if (len(top_k_close) < K):
                top_k_close.append(coord_id)
        if ((len(top_k_close) == K) and (len(top_k_open) == K) and (len(top_k_global) == K)):
            break
    return {"global":top_k_global, "open":top_k_open, "close":top_k_close}


def copy_top_k(cfg):
    """Make copies of top k in proper directories"""
    top_k = get_top_k(cfg)
    dirs = get_top_k_dirs(cfg)
    for key in top_k:
        make_directory(dirs[key])
        k = 0
        for coord_id in top_k[key]:
            in_pdb = get_in_pdb_path(cfg, coord_id)
            #out_pdb = os.path.join(dirs[key], str(k).zfill(Z_FILL) + "-" + coord_id + ".pdb")
            out_pdb = os.path.join(dirs[key], str(k).zfill(Z_FILL) + ".pdb")
            print("Copying " + in_pdb + " to " + out_pdb)
            k = k + 1
            copyfile(in_pdb, out_pdb)


#############################################################################
# Main
#############################################################################


def __main__():
    for mut in MUTATION_TYPES:
        for pH in PH_LEVELS:
            for probe in PROBE_RADII:
                cfg = OcConfig(mut=mut,
                               pH=pH,
                               probe=probe)
                copy_top_k(cfg)


if __name__ == '__main__':
    __main__()