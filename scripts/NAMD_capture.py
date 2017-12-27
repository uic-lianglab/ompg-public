#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script will score and capture interactions for a coordinate file using NAMD.

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
NAMD_PATH = os.path.join(NAMD_DIR, 'bin', 'nix', 'namd2_0_step_ir')

# Filter interactions which are not within loop regions
# @TODO - read these values from file
LOOPS = [xrange(18, 29+1), xrange(54, 65+1), xrange(97, 106+1),
         xrange(177, 188+1), xrange(217, 234+1), xrange(259, 267+1)]

# Indexed entries within an interaction record row
IX_IR_RES_A = 0
IX_IR_RES_B = 1
IX_IR_ELECTRO = 2
IX_IR_GBIS = 3

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


def get_energy_dir(output_base_dir):
    """@return path for storing energies."""
    return os.path.join(output_base_dir, 'ener')


def get_interaction_dir(output_base_dir):
    """@return path for storing interactions."""
    return os.path.join(output_base_dir, 'ints')


def get_tcl_dir(output_base_dir):
    """@return path for storing tcl scripts."""
    return os.path.join(output_base_dir, 'tcl')


def get_score_tcl(coord_path, struct_path):
    """@return tcl script for 0-step scoring run with NAMD."""
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
cutoff 16.0
switching off
pairlistdist 24
PME no
fixedAtoms on
outputName not_used
binaryoutput no
run 0
""" % (coord_path, struct_path, FORCE_FIELD_PATH)
    return tcl


def involves_loop(record):
    """Check if residue pair involves interacting loop region."""
    for loop in LOOPS:
        if (int(record[IX_IR_RES_A]) in loop):
            return True
        elif (int(record[IX_IR_RES_B]) in loop):
            return True
    return False


def write_energies(tokens_stdout, energy_dir, coord_id):
    """Write energies to disk."""
    ix_etitle = tokens_stdout.index("ETITLE:")
    ix_energy = tokens_stdout.index("ENERGY:", ix_etitle)
    # ix_e_ts = ix_energy + 1
    ix_e_bond = ix_energy + 2
    ix_e_angle = ix_energy + 3
    ix_e_dihed = ix_energy + 4
    ix_e_imprp = ix_energy + 5
    ix_e_elect = ix_energy + 6
    ix_e_vdw = ix_energy + 7
    ix_e_boundary = ix_energy + 8
    ix_e_misc = ix_energy + 9
    # ix_e_kinetic = ix_energy + 10
    # ix_e_total = ix_energy + 11
    # ix_e_temp = ix_energy + 12
    ix_e_potential = ix_energy + 13
    # ix_e_total3 = ix_energy + 14
    # ix_e_tempavg = ix_energy + 15

    csv_profile = tokens_stdout[ix_e_potential] + "," + \
        tokens_stdout[ix_e_elect] + "," + \
        tokens_stdout[ix_e_vdw] + "," + \
        tokens_stdout[ix_e_bond] + "," + \
        tokens_stdout[ix_e_angle] + "," + \
        tokens_stdout[ix_e_dihed] + "," + \
        tokens_stdout[ix_e_imprp] + "," + \
        tokens_stdout[ix_e_boundary] + "," + \
        tokens_stdout[ix_e_misc]
        
    # Write energy profile to disk
    energy_file = coord_id + ".ener.csv"
    energy_path = os.path.join(energy_dir, energy_file)
    with open(energy_path, 'w') as ener:
        ener.write(csv_profile + '\n')


def write_interactions(tokens_stdout, interaction_dir, coord_id):
    """Write interactions involving loop regions to disk."""
    # Extract interactions
    ix_int_begin = tokens_stdout.index("BEGIN_INTERACTION_RECORDS")
    ix_int_end = tokens_stdout.index("END_INTERACTION_RECORDS",
                                     ix_int_begin)

    # Write interactions to disk
    interaction_file = coord_id + ".ints.csv"
    interaction_path = os.path.join(interaction_dir, interaction_file)
    with open(interaction_path, 'w') as ints:
        # Eat begin token and header row
        for i in xrange(ix_int_begin+2, ix_int_end):
            # Process each interaction tuple
            record = tokens_stdout[i].split(',')
            # Check if loop regions are involved
            if (involves_loop(record)):
                ints.write(tokens_stdout[i] + '\n')


def capture(output_base_dir, input_coord_path, input_psf_path):
    """Capture energy and pairwise interactions"""
    # Create folders if not existent.
    energy_dir = get_energy_dir(output_base_dir)
    make_directory(energy_dir)
    interaction_dir = get_interaction_dir(output_base_dir)
    make_directory(interaction_dir)
    tcl_dir = get_tcl_dir(output_base_dir)
    make_directory(tcl_dir)
    
    # Input coordinates file name and identifier
    coord_file = os.path.basename(input_coord_path)
    coord_id = os.path.splitext(coord_file)[0]
    
    # Generate tcl configuration script
    tcl_body = get_score_tcl(input_coord_path, input_psf_path)
    tcl_file = coord_id + ".tcl"
    tcl_path = os.path.join(tcl_dir, tcl_file)
    with open(tcl_path, "w") as script:
        script.write(tcl_body)

    # Call NAMD for 0-step scoring run
    proc_stdout = subprocess.check_output([NAMD_PATH, tcl_path])

    # Low on disk space, so remove any unneeded files
    os.remove(tcl_path)

    # Tokenize output
    tokens_stdout = proc_stdout.split()

    # Parse output and write to disk
    write_energies(tokens_stdout, energy_dir, coord_id)
    write_interactions(tokens_stdout, interaction_dir, coord_id)


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

    # Capture energy and interactions
    capture(output_base_dir,
            input_coord_path,
            input_psf_path)

if __name__ == '__main__':
    __main__()
