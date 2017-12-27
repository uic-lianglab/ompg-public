# Usage:
# python <path_to_this_script> <input_pdb_dir> <output_png_dir>
# Script process an input folder of PDBs and renders their occlusion
# using Pymol.

#############################################################################
# Imports
#############################################################################

import os
import sys
import platform
import subprocess

#############################################################################
# Globals
#############################################################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Expected command line arguments
ARG_IN_PDB_DIR = 1
ARG_OUT_PNG_DIR = 2
CMD_LINE_SIZE = 3

# Determine executable name
EXEC_NAME = "pymol"
# On windows, pymol has problems initializing the python open gl module,
# so, we must call pymolwin instead.
# http://stackoverflow.com/questions/1387222/reliably-detect-windows-in-python
if (platform.system() == 'Windows'):
    EXEC_NAME = "pymolwin"

# Change here to switch between surface and cartoon rendering
PYMOL_RENDER_SCRIPT_TYPE = "surface" # change this to "cartoon" to switch rendering type
PYMOL_RENDER_SCRIPT_NAME = "pymol_render_occlusion_" + PYMOL_RENDER_SCRIPT_TYPE + ".py"

# Path to script for rendering in PYMOL
PYMOL_RENDER_SCRIPT_PATH = os.path.join(SCRIPT_DIR, PYMOL_RENDER_SCRIPT_NAME)


# DEFAULT ARGUMENTS
DEF_ROOT_DIR =  os.path.join(SCRIPT_DIR, "..")
DEF_ROOT_DIR = os.path.abspath(DEF_ROOT_DIR)
DEF_OUTPUT_DIR = os.path.join(DEF_ROOT_DIR, "ompg", "output")
DEF_RELAX_DIR = os.path.join(DEF_OUTPUT_DIR, "relax.old")
DEF_CAPT_BASE_DIR = os.path.join(DEF_RELAX_DIR, "capt")
DEF_IN_PDB_DIR = os.path.join(DEF_CAPT_BASE_DIR, "pH5", "wt", "regsc",
                              "sids.merge.close.PDBs")
DEF_OUT_PNG_DIR = os.path.join(DEF_CAPT_BASE_DIR, "pH5", "wt", "regsc",
                               "sids.merge.close.PNGs")

#############################################################################
# Methods
#############################################################################


# http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Will make intermediate directories
def make_directory(dir_path):
    """Utility for making a directory if not existing."""
    abs_dir_path = os.path.abspath(dir_path)
    if not os.path.exists(abs_dir_path):
        os.makedirs(abs_dir_path)


def get_fid(fname):
    """Get basename from file path with extension"""
    return os.path.splitext(fname)[0]


def batch_render_occlusions(in_pdb_dir, out_png_dir):
    # Make sure output directory exists
    make_directory(out_png_dir)

    abs_pdb_dir = os.path.abspath(in_pdb_dir)
    abs_png_dir = os.path.abspath(out_png_dir)

    print "Processing PDBs in: " + in_pdb_dir
    for fname in os.listdir(in_pdb_dir):
        pdb_path = os.path.join(abs_pdb_dir, fname)
        fid = get_fid(fname)
        png_name = fid + ".png"
        png_path = os.path.join(abs_png_dir, png_name)
        # Call pymol
        proc_stdout = subprocess.check_output([EXEC_NAME, '-qrc', PYMOL_RENDER_SCRIPT_PATH, '--', pdb_path, png_path])
        print proc_stdout


#############################################################################
# Main
#############################################################################


def __main__():
    print "=========== BATCH RENDER OCCLUSIONS ==========="
    # Check if user arguments exist
    in_pdb_dir = DEF_IN_PDB_DIR
    out_png_dir = DEF_OUT_PNG_DIR
    arg_count = len(sys.argv)
    if (arg_count == CMD_LINE_SIZE):
        in_pdb_dir = sys.argv[ARG_IN_PDB_DIR]
        out_png_dir = sys.argv[ARG_OUT_PNG_DIR]

    print "\t-in_pdb_dir: " + in_pdb_dir
    print "\t-out_png_dir: " + out_png_dir

    # Munge interaction records
    print "Using EXEC_NAME: " + EXEC_NAME
    batch_render_occlusions(in_pdb_dir,
                            out_png_dir)

    print "Finished."

if __name__ == '__main__':
    __main__()