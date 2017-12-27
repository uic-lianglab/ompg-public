#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Convert NAMD ascii coordinate output to pdb format for use with CASTp.

Will strip hydrogens and rename atoms according to template pdb.

Usage:

python <path_to_script> <out_pdb_path> <in_coor_path> <template_pdb_path>
"""

#############################################################################
# Imports
#############################################################################

import os
import sys

#############################################################################
# Globals
#############################################################################

# Path to script directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Command line arguments are in this order within sys.argv list
ARG_OUT_PDB_PATH = 1
ARG_IN_COOR_PATH = 2
ARG_TEMPLATE_PDB_PATH = 3
CMD_LINE_SIZE = 4

# Element indices within PDB ATOM RECORD
# http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
# python slice coords are in []
IX_PDB_RECORD_NAME = 0  # char 1-6 [0:6]
IX_PDB_ATOM_SERIAL_NO = 1  # char 7-11 [6:11]
IX_PDB_ATOM_NAME = 2  # char 13-16 [12:16]
IX_PDB_ALT_LOC = 3  # char 17 [16]
IX_PDB_RES_NAME = 4  # char 18-20 [17:20]
IX_PDB_CHAIN_ID = 5  # char 22 [21]
IX_PDB_RES_SEQ_NO = 6  # char 23-26 [22:26]
IX_PDB_RES_INS = 7  # char 27 [26]
IX_PDB_ATOM_X = 8  # char 31-38 [30:38]
IX_PDB_ATOM_Y = 9  # char 39-46 [38:46]
IX_PDB_ATOM_Z = 10  # char 47-54 [46:54]
IX_PDB_OCC = 11  # char 55-60 [54:60]
IX_PDB_TEMP_FACT = 12  # char 61-66 [60:66]
IX_PDB_SEG_IDENT = 13  # char 73-76 [72:76]
IX_PDB_ELEM_SYM = 14  # char 77-78 [76:78]
IX_PDB_ATOM_CHARGE = 15  # char 79-80 [78:80]


#############################################################################
# Methods
#############################################################################


def get_atom_record(line):
    """Parse PDB atom record into a list."""
    # WARNING: Does not check if actually an ATOM record!
    rec = [line[0:6],  # record name
           line[6:11],  # atom serial no
           line[12:16],  # atom name
           line[16],  # alt location
           line[17:20],  # residue name
           line[21],  # chain id
           line[22:26],  # residue sequence no
           line[26],  # residue insertion code
           line[30:38],  # atom x-coord
           line[38:46],  # atom y-coord
           line[46:54],  # atom z-coord
           line[54:60],  # atom occupancy
           line[60:66],  # temp factor
           line[72:76],  # segment identifier
           line[76:78],  # element symbol
           line[78:80]]  # charge
    return rec


def is_hydrogen(atom_rec):
    """@return TRUE if record is for hydrogen atom, FALSE o/w."""
    return atom_rec[IX_PDB_ATOM_NAME].strip().lower().startswith("h")


def load_template(pdb_path):
    """Load template pdb into record tuples."""
    atom_map = {}
    with open(pdb_path) as fpdb:
        for line in fpdb:
            # Only process ATOM records
            if line.startswith('ATOM'):
                rec = get_atom_record(line)
                # Catch special cases:
                if (("ILE" in rec[IX_PDB_RES_NAME]) and
                   ("CD" in rec[IX_PDB_ATOM_NAME])):
                    # Map coord ILE->CD to ILE->CD1
                    key = (rec[IX_PDB_RES_SEQ_NO], " CD ")
                if (("ILE" in rec[IX_PDB_RES_NAME]) and
                   ("CD" in rec[IX_PDB_ATOM_NAME])):
                    pass
                else:
                    key = (rec[IX_PDB_RES_SEQ_NO], rec[IX_PDB_ATOM_NAME])

                atom_map[key] = rec
    return atom_map


def rec_to_str(rec):
    """Convert pdb atom record list to a line."""
    return rec[IX_PDB_RECORD_NAME] + \
        rec[IX_PDB_ATOM_SERIAL_NO] + \
        " " + \
        rec[IX_PDB_ATOM_NAME] + \
        rec[IX_PDB_ALT_LOC] + \
        rec[IX_PDB_RES_NAME] + \
        " " + \
        rec[IX_PDB_CHAIN_ID] + \
        rec[IX_PDB_RES_SEQ_NO] + \
        rec[IX_PDB_RES_INS] + \
        "   " + \
        rec[IX_PDB_ATOM_X] + \
        rec[IX_PDB_ATOM_Y] + \
        rec[IX_PDB_ATOM_Z] + \
        rec[IX_PDB_OCC] + \
        rec[IX_PDB_TEMP_FACT] + \
        "      " + \
        rec[IX_PDB_SEG_IDENT] + \
        rec[IX_PDB_ELEM_SYM] + \
        rec[IX_PDB_ATOM_CHARGE]


def process_with_template(coor_path, template_pdb):
    """Strip hydrogens, renumber atoms, and rename residues."""
    coor_lines = []
    with open(coor_path, "r") as fcoor:
        coor_lines = fcoor.readlines()

    castp_lines = []

    for coor_line in coor_lines:
        if coor_line.startswith('ATOM'):
            coor_rec = get_atom_record(coor_line)
            # Strip hydrogens
            if is_hydrogen(coor_rec):
                continue
            # Check if template atom record exists
            key = (coor_rec[IX_PDB_RES_SEQ_NO], coor_rec[IX_PDB_ATOM_NAME])
            if key in template_pdb:
                template_rec = template_pdb[key]
                # Renumber atom
                coor_rec[IX_PDB_ATOM_SERIAL_NO] = \
                    template_rec[IX_PDB_ATOM_SERIAL_NO]
                # Rename residue
                coor_rec[IX_PDB_RES_NAME] = \
                    template_rec[IX_PDB_RES_NAME]
                # Append new (possibly out of order) atom record
                castp_line = rec_to_str(coor_rec)
                atom_no = int(template_rec[IX_PDB_ATOM_SERIAL_NO].strip())
                castp_tup = (atom_no, castp_line)
                castp_lines.append(castp_tup)
            else:
                print "Warning: missing template atom record for:"
                print coor_line

    # Sort by atom number
    castp_lines.sort(key=lambda tup: tup[0])
    # Make into list of lines (remove tuple)
    castp_lines = [tup[1] for tup in castp_lines]
    return castp_lines


def coor_to_castp_pdb(out_pdb_path, in_coor_path, template_pdb_path):
    """Convert coor NAMD ascii format to pdb format used by CASTp.

    Strips hydrogens and renames atoms and residues according to template pdbs.
    """
    # Load templates
    template = load_template(template_pdb_path)

    lines = process_with_template(in_coor_path, template)
    
    with open(out_pdb_path, "w") as fpdb:
        fpdb.writelines(lines)


#############################################################################
# Main
#############################################################################


def __main__():
    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count != CMD_LINE_SIZE):
        sys.exit("Expected " + str(CMD_LINE_SIZE) + " arguments but received " + str(arg_count) + ".")

    out_pdb_path = sys.argv[ARG_OUT_PDB_PATH]
    in_coor_path = sys.argv[ARG_IN_COOR_PATH]
    template_pdb_path = sys.argv[ARG_TEMPLATE_PDB_PATH]

    print "============ Coor to CASTp PDB ============"
    print "\t-out_pdb_path: " + out_pdb_path
    print "\t-in_coor_path: " + in_coor_path
    print "\t-template_pdb_path: " + template_pdb_path

    # Convert pdb format
    coor_to_castp_pdb(out_pdb_path,
                      in_coor_path,
                      template_pdb_path)

    print "Finished conversion."

if __name__ == '__main__':
    __main__()
