#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script will minimize potentials for generated hydrogen atoms as well as simulated loop side chain atoms.

Usage:

python <path_to_script> <output_base_dir> <input_coord_path> <input_psf_path>
"""

#############################################################################
# Imports
#############################################################################

import os
import subprocess
import sys

#############################################################################
# Globals
#############################################################################

# Path to script directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Path to root project directory
ROOT_DIR = os.path.join(SCRIPT_DIR, '..')

# Path to base NAMD directory
NAMD_DIR = os.path.join(ROOT_DIR, 'tools', 'namd')

# The force field to use
FORCE_FIELD = "charmm36"
FORCE_FIELD_PATH = os.path.join(NAMD_DIR, 'forcefield', FORCE_FIELD, 'params.prm')

# Path to NAMD binary
NAMD_PATH = os.path.join(NAMD_DIR, 'bin', 'nix', 'namd2')

# Start residue identifiers for simulated loop regions
# Note: pair (10000, 10001) is not real loop region but makes
# code simpler than checking if we run past end of loop
# Note: These arrays must be sorted in ascending order!
# @TODO - read these values from file
LOOP_START = [18, 54, 97, 177, 217, 259, 10000]
# End residue identifiers for simulated loop regions
LOOP_END = [29, 65, 106, 188, 234, 267, 10001]

# Command line arguments are in this order within sys.argv list
ARG_OUTPUT_BASE_DIR = 1
ARG_INPUT_COORD_PATH = 2
ARG_INPUT_PSF_PATH = 3
CMD_LINE_SIZE = 4

#############################################################################
# Methods
#############################################################################

# http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
def make_directory(directory):
    """Utility for making a directory if not existing."""
    if not os.path.exists(directory):
        os.makedirs(directory)


def remove_bak_file(fpath_no_bak_ext):
    """Remove back-up version of file."""
    bak_fpath = fpath_no_bak_ext + ".BAK"
    if (os.path.exists(bak_fpath)):
        os.remove(bak_fpath)


def remove_misc_file(dir, id, ext):
    """Utility for removing an unneeded file."""
    fpath = os.path.join(dir, id + ext)
    os.remove(fpath)
    remove_bak_file(fpath)


def get_sim_dir_name():
    """@return name of simulation directory."""
    return 'sim'


def get_sim_dir(output_base_dir):
    """@return path to directory for minimization simulation data."""
    return os.path.join(output_base_dir, get_sim_dir_name())


def get_tcl_dir(output_base_dir):
    """@return path for storing tcl scripts used in minimization."""
    return os.path.join(output_base_dir, 'tcl')


def get_fixed_pre_dir(output_base_dir):
    """@return path for storing pdb file with loop backbone regions fixed."""
    return os.path.join(output_base_dir, 'fixed_pre')


def get_fixed_post_dir(output_base_dir):
    """@return path for storing coordinate files with fixed non-loops regions."""
    return os.path.join(output_base_dir, 'fixed_post')


def get_minimize_tcl(coord_path, struct_path, output_prefix):
    """@return tcl script for minimizing atom positions with NAMD."""
    # Note: minimize <steps> must be multiple of stepspercycle else NAMD errors out
    tcl = """coordinates %s
structure %s
paraTypeCharmm on
parameters %s
temperature 300.0
gbis on
alphaCutoff 15.0
ionConcentration 1.0
GBISDelta 0.8
GBISBeta 0.0
GBISGamma 2.90912
exclude scaled1-4
1-4scaling 1.0
cutoff 16
switching off
pairlistdist 24
PME no
fixedAtoms on
outputName %s
binaryoutput no
stepspercycle 3
minimize 66
""" % (coord_path, struct_path, FORCE_FIELD_PATH, output_prefix)
    return tcl


def fix_non_loop_regions(min_coord_path, fixed_coord_path):
    """Post-process a minimized coordinate file by fixing non-loop regions."""
    min_lines = []
    with open(min_coord_path, "r") as f_min:
        min_lines = f_min.readlines()
    curr_loop_ix = 0
    curr_loop_start = LOOP_START[0]
    curr_loop_end = LOOP_END[0]
    fixed_lines = []
    # Process each record in coordinate file
    # Note: min_h_line should contain a newline character already
    for min_line in min_lines:
        fixed_line = min_line
        # Only process ATOM records
        if min_line.startswith('ATOM'):
            # Residue identifier is from indices 22 through 25
            res_id = int(min_line[22:26])
            # Default occupancy is unsimulated (fixed)
            occ = "1.00"
            # Check if current loop region needs to be updated
            if (res_id > curr_loop_end):
                curr_loop_ix = curr_loop_ix + 1
                curr_loop_start = LOOP_START[curr_loop_ix]
                curr_loop_end = LOOP_END[curr_loop_ix]
            # Check if in loop region -> set to non-fixed atom
            if ((res_id >= curr_loop_start) and (res_id <= curr_loop_end)):
                occ = "0.00"
            # Occupancy is from indices 56 through 59
            # (Note: technically occupancy is from 54 to 59,
            #  but we only need last 4 characters)
            # http://stackoverflow.com/questions/1228299/change-one-character-in-a-string-in-python/1228332#1228332
            fixed_line = min_line[:56] + occ + min_line[60:]
        fixed_lines.append(fixed_line)

    # Write to disk
    with open(fixed_coord_path, "w") as f_fixed:
        f_fixed.writelines(fixed_lines)


def fix_bb(coord_path, coord_file, fixed_bb_dir):
    """Fix all backbone atoms."""
    raw_prot_lines = []
    with open(coord_path, "r") as f_raw_prot:
        raw_prot_lines = f_raw_prot.readlines()
    fixed_lines = []
    # Process each record in coordinate file
    # Note: raw_prot_line should contain a newline character already
    for raw_prot_line in raw_prot_lines:
        fixed_line = raw_prot_line
        # Only process ATOM records
        if raw_prot_line.startswith('ATOM'):
            # Atom name is from indices 12-15
            atom_name = raw_prot_line[12:16]
            atom_name = atom_name.strip()
            # Default occupancy is unsimulated (fixed)
            occ = "1.00"
            if atom_name not in ["N", "CA", "C", "O"]:
                occ = "0.00"
            # Occupancy is from indices 56 through 59
            # (Note: technically occupancy is from 54 to 59,
            #  but we only need last 4 characters)
            # http://stackoverflow.com/questions/1228299/change-one-character-in-a-string-in-python/1228332#1228332
            fixed_line = raw_prot_line[:56] + occ + raw_prot_line[60:]
        fixed_lines.append(fixed_line)

    # Write to disk
    out_path = os.path.join(fixed_bb_dir, coord_file)
    with open(out_path, "w") as f_fixed:
        f_fixed.writelines(fixed_lines)

    return out_path


def minimize(output_base_dir, input_coord_path, input_psf_path):
    """Minimize side chain atom positions using NAMD for (.pdb, .psf) pair.

    Will also unfix loop side chains prior to minimzation as well as fix
    non-loop regions post minimization. Fixing of loop atoms is done by
    changing the occupancy column of the PDB file which NAMD will interpret
    as a fixed atom.
    """

    # Create directories
    sim_dir = get_sim_dir(output_base_dir)
    make_directory(sim_dir)
    tcl_dir = get_tcl_dir(output_base_dir)
    make_directory(tcl_dir)
    fixed_pre_dir = get_fixed_pre_dir(output_base_dir)
    make_directory(fixed_pre_dir)
    fixed_post_dir = get_fixed_post_dir(output_base_dir)
    make_directory(fixed_post_dir)

    # Determine base file name for minimized data
    raw_coord_file = os.path.basename(input_coord_path)
    raw_coord_id = os.path.splitext(raw_coord_file)[0]
    min_coord_id = raw_coord_id + ".min"

    # Skip minimization if fixed version already exists
    fixed_coord_file = min_coord_id + ".fixed.coor"
    fixed_coord_path = os.path.join(fixed_post_dir, fixed_coord_file)
    if (os.path.exists(fixed_coord_path)):
        print "Skipping " + input_coord_path
        return

    # Fix backbone atoms within loop regions (overwrite coord path)
    input_coord_path = fix_bb(input_coord_path, raw_coord_file, fixed_pre_dir)

    # Note: am using relative output paths within tcl script because:
    # 1.) NAMD switches to this parameter tcl script directory and
    #   sets it as the working directory
    # 2.) NAMD truncates output prefixes that are too long - which can
    #   cause files to be overwritten as they end up with the same name.
    output_prefix = os.path.join('..', get_sim_dir_name(), min_coord_id)
    # Generate tcl configuration script
    tcl_body = get_minimize_tcl(input_coord_path, input_psf_path,
                                output_prefix)
    tcl_file = min_coord_id + ".tcl"
    tcl_path = os.path.join(tcl_dir, tcl_file)
    with open(tcl_path, "w") as script:
        script.write(tcl_body)
    # Call NAMD to minimize hydrogens
    # Note: consider tossing stdout due to disk space concerns when logging
    namd_stdout = subprocess.check_output([NAMD_PATH, tcl_path])
    print namd_stdout

    # Avoid future calculations of non-loop to non-loop interaction
    # energies by fixing these atoms.
    min_coord_path = os.path.join(sim_dir, min_coord_id + ".coor")
    fix_non_loop_regions(min_coord_path, fixed_coord_path)

    # We are low on disk space, so make sure to remove unneeded files
    remove_misc_file(sim_dir, min_coord_id, ".vel")
    remove_misc_file(sim_dir, min_coord_id, ".xsc")
    remove_misc_file(sim_dir, min_coord_id, ".coor")
    os.remove(tcl_path)
    # Remove generated pre-minimization file
    os.remove(input_coord_path)


#############################################################################
# Main
#############################################################################


def __main__():
    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count != CMD_LINE_SIZE):
        sys.exit("Expected " + str(CMD_LINE_SIZE) + " arguments but received " + str(arg_count) + ".")

    output_base_dir = sys.argv[ARG_OUTPUT_BASE_DIR]
    input_coord_path = sys.argv[ARG_INPUT_COORD_PATH]
    input_psf_path = sys.argv[ARG_INPUT_PSF_PATH]

    # Minimize clashes
    minimize(output_base_dir,
             input_coord_path,
             input_psf_path)


if __name__ == '__main__':
    __main__()
