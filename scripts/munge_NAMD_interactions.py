#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script will parse interaction records and merge into a master CSV

Usage:

python <path_to_script> <in_ints_dir> <in_ener_csv_path> <out_combined_csv_path> <max_energy_rank>
"""

#############################################################################
# Imports
#############################################################################

import os
import sys

from collections import defaultdict

#############################################################################
# Globals
#############################################################################


# Command line arguments are in this order within sys.argv list
ARG_IN_INTS_DIR = 1
ARG_IN_ENER_CSV_PATH = 2
ARG_OUT_COMB_CSV_PREFIX = 3
ARG_OUT_MAX_ENERGY_RANK = 4
CMD_LINE_SIZE = 5

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

# Indexed entries within an interaction record row
IX_IR_RES_A = 0
IX_IR_RES_B = 1
IX_IR_ELECTRO = 2
IX_IR_GBIS = 3

#############################################################################
# Methods
#############################################################################


def get_fid(fname):
    # extract <fid>.ener.csv
    return os.path.splitext(os.path.splitext(fname)[0])[0]


def load_sorted_energy_ids(ener_csv_path):
    # Extract energy ranking
    ener_ids=[]
    # Lines in file assumed to be sorted according to an energy criteria
    # (e.g. lowest potential energy)
    with open(ener_csv_path, "r") as f:
        # Skip header row
        f.next()
        for line in f:
            if ',' in line:
                rec = line.strip().split(',')
                ener_ids.append(rec[IX_E_NAME])
    return ener_ids


def write_records(records, pair_keys, ener_ids, output_path):
    """Write interaction records to disk."""
    # Ensures that all interactions are homogeneous:
    # By iterating over a set of interaction pair keys, if a sample is missing
    # a cutoff interaction, it will be written as 'NA'
    print "Writing to: " + output_path
    with open(output_path, "w") as f:
        # Write header
        header = "RESA,RESB"
        sample_names = ','.join(map(str, ener_ids))
        header = header + ',' + sample_names
        f.write(header + '\n')
        # Write interaction energies
        for res_pair in pair_keys:
            # Skip self interactions (not interesting)
            if (res_pair[IX_IR_RES_A] == res_pair[IX_IR_RES_B]):
                continue
            f.write(res_pair[IX_IR_RES_A] + ',' + res_pair[IX_IR_RES_B])
            for rec in records:
                # defaultdict will output 'NA' if this pair wasn't recorded for
                # this sample due to cutoff distance setting
                f.write(',' + rec[res_pair])
            f.write('\n')


def munge_interactions(in_ints_dir, in_ener_csv_path, out_comb_csv_prefix, max_energy_rank):
    # Load energy identifiers
    # Interaction records will be listed in same column order
    ener_ids = load_sorted_energy_ids(in_ener_csv_path)

    # Crop energy identifiers to maximum rank allowed
    if (max_energy_rank <= 0) or (max_energy_rank > len(ener_ids)):
        max_energy_rank = len(ener_ids)
    ener_ids = ener_ids[:max_energy_rank]

    # Create hash table of energy identifiers for quick membership testing
    ener_ids_set = frozenset(ener_ids)

    # Set of unique recorded interaction pairs.
    # Used to output a consistent set of interaction pairs as, due to cutoff
    # distances, not all interaction pairs are recorded for each each loop.
    pair_keys = set()

    # Output interaction sets for all processed tuples
    # Map is from energy identifier -> set of interactions
    # Electrostatic potential = Coulomb + GBIS
    record_electrostatic = {}
    # GBIS potential = effect of atom screening
    record_GBIS = {}

    # Load records
    print "Processing records in: " + in_ints_dir
    for fname in os.listdir(in_ints_dir):
        if fname.endswith(".ints.csv"):
            # Extract identifier
            fid = get_fid(fname)
            # Skip if sample is not in energy set
            if fid not in ener_ids_set:
                continue
            # Full path to input file
            fpath = os.path.join(in_ints_dir, fname)
            # Parse residue-residue interactions for this sample.
            # Create seperate dictionaries for total electrostatic and gbis
            # respectively. Also, if a residue pair was not encountered,
            # this will write 'NA' to the final output file.
            recs_es = defaultdict(lambda: 'NA')
            recs_gbis = defaultdict(lambda: 'NA')

            # Extract interactions
            with open(fpath, "r") as f:
                for line in f:
                    if ',' in line:
                        rec = line.strip().split(',')
                        rec_key = (rec[IX_IR_RES_A], rec[IX_IR_RES_B])
                        # Register interaction pair
                        pair_keys.add(rec_key)
                        recs_es[rec_key] = rec[IX_IR_ELECTRO]
                        recs_gbis[rec_key] = rec[IX_IR_GBIS]

            # Store interactions for this sample
            record_electrostatic[fid] = recs_es
            record_GBIS[fid] = recs_gbis

    # Convert universal set of interaction pairs into a sorted list
    pair_keys = list(pair_keys)
    pair_keys.sort(key=lambda tup: (int(tup[0]), int(tup[1])))

    # Convert interaction sets map to sorted list
    record_electrostatic_sorted = []
    record_GBIS_sorted = []
    for ener_id in ener_ids:
        record_electrostatic_sorted.append(record_electrostatic[ener_id])
        record_GBIS_sorted.append(record_GBIS[ener_id])

    if record_electrostatic_sorted:
        write_records(record_electrostatic_sorted,
                      pair_keys,
                      ener_ids,
                      out_comb_csv_prefix + "_electro." + str(max_energy_rank) + ".csv")

    if record_GBIS_sorted:
        write_records(record_GBIS_sorted,
                      pair_keys,
                      ener_ids,
                      out_comb_csv_prefix + "_gbis." + str(max_energy_rank) + ".csv")


#############################################################################
# Main
#############################################################################


def __main__():
    print "=========== MUNGE NAMD INTERACTIONS ==========="
    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count != CMD_LINE_SIZE):
        sys.exit("Expected " + str(CMD_LINE_SIZE) + " arguments but received " + str(arg_count) + ".")

    in_ints_dir = sys.argv[ARG_IN_INTS_DIR]
    in_ener_csv_path = sys.argv[ARG_IN_ENER_CSV_PATH]
    out_comb_csv_prefix = sys.argv[ARG_OUT_COMB_CSV_PREFIX]
    max_energy_rank = int(sys.argv[ARG_OUT_MAX_ENERGY_RANK])

    print "\t-in_ints_dir: " + in_ints_dir
    print "\t-in_ener_csv_path: " + in_ener_csv_path
    print "\t-out_comb_csv_prefix: " + out_comb_csv_prefix
    print "\t-max_energy_rank: "+ str(max_energy_rank)

    # Munge interaction records
    munge_interactions(in_ints_dir,
                       in_ener_csv_path,
                       out_comb_csv_prefix,
                       max_energy_rank)

    print "Finished."

if __name__ == '__main__':
    __main__()
