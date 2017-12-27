#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script will combine all .ener.csv files in a folder

Usage:

python <path_to_script> <in_ener_dir> <out_combined_csv_path>
"""

#############################################################################
# Imports
#############################################################################

import os
import sys


#############################################################################
# Globals
#############################################################################


# Command line arguments are in this order within sys.argv list
ARG_IN_ENER_DIR = 1
ARG_OUT_COMB_CSV_PATH = 2
CMD_LINE_SIZE = 3


# Column indices for output energy file
IX_E_NAME = 0
IX_E_PH = 1
IX_E_POTENTIAL = 2
IX_E_ELECT = 3
IX_E_VDW = 4
IX_E_BOND = 5
IX_E_ANGLE = 6
IX_E_DIHED = 7
IX_E_IMPRP = 8
IX_E_BOUNDARY = 9
IX_E_MISC = 10

#############################################################################
# Methods
#############################################################################


def get_fid(fname):
    # extract <fid>.ener.csv
    return os.path.splitext(os.path.splitext(fname)[0])[0]


def cat_energy(in_ener_dir, out_comb_csv_path):
    scores = []
    # Load records
    for fname in os.listdir(in_ener_dir):
        if fname.endswith(".ener.csv"):
            # Extract identifier
            fid = get_fid(fname)
            ix_pH_start = fid.index('.pH') + 3
            ix_pH_end = fid.index('.', ix_pH_start)
            # Extract pH
            pH = fid[ix_pH_start : ix_pH_end]
            # Full path to input file
            fpath = os.path.join(in_ener_dir, fname)
            rec=[fid, pH]
            # Extract energies
            with open(fpath, "r") as f:
                for line in f:
                    if ',' in line:
                        rec.extend(line.strip().split(','))
            scores.append(rec)

    # Early out if no records found
    if not scores:
        return

    # Sort by potential
    scores.sort(key=lambda rec: float(rec[IX_E_POTENTIAL]))

    # Write records
    with open(out_comb_csv_path, "w") as f:
        # Write header
        f.write('NAME,PH,POTENTIAL,ELECT,VDW,BOND,ANGLE,DIHED,IMPRP,BOUNDARY,MISC\n')
        for rec in scores:
            f.write(','.join(map(str, rec)) + '\n')


#############################################################################
# Main
#############################################################################


def __main__():
    print "=========== CAT NAMD ENERGY ==========="
    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count != CMD_LINE_SIZE):
        sys.exit("Expected " + str(CMD_LINE_SIZE) + " arguments but received " + str(arg_count) + ".")

    in_ener_dir = sys.argv[ARG_IN_ENER_DIR]
    out_comb_csv_path = sys.argv[ARG_OUT_COMB_CSV_PATH]

    print "\t-in_ener_dir: " + in_ener_dir
    print "\t-out_comb_csv_path: " + out_comb_csv_path

    # Munge energy files
    cat_energy(in_ener_dir,
               out_comb_csv_path)

    print "Finished."

if __name__ == '__main__':
    __main__()
