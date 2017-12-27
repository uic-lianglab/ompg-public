#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes:
# - .c2s.csv: file containing cluster identifier to sample mappings
# - .ex.csv: file containing exemplar identifiers for each cluster
# Script resamples from each cluster according to
# - min_coverage: the minimum number of samples to obtain from each cluster
# - min_output_count: the minimum total number of samples to write to disk
# The final written count samples is:
#   max(min_coverage * num_clusters, min_output_count)
# Script prioritizes sampling without replacement; however once a cluster has
#   had all its unique sample ids selected, then sampling with replacement will ensue
# After all samples have been selected, the resulting set is shuffled once to
#   locally decorrelate samples
# Also, script will output a silhouette style plot to visualize how many samples
#   were selected from each cluster

###########################################
# Imports
###########################################

# For parsing user supplied arguments
import argparse

# For parsing list of files in a directory
import os

# For generating random numbers
import random

###########################################
# Default parameters
###########################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Default path to directory containing .c2s and .ex files
DEFAULT_DATA_DIR = \
    os.path.join(SCRIPT_DIR, "../ompg/output/frag_libs/f_clust_libs/loop_6_wt/2iwv")

# Default data set to process
DEFAULT_DATASET_NAME = "loop_6_wt.2iwv.frag.lib"

# Default cluster to sample mapping file path
DEFAULT_C2S_FILE = \
    os.path.join(DEFAULT_DATA_DIR, DEFAULT_DATASET_NAME + ".c2s.csv")

# Default exemplar file path
DEFAULT_EX_FILE = \
    os.path.join(DEFAULT_DATA_DIR, DEFAULT_DATASET_NAME + ".ex.csv")

# Default out resampled file path
DEFAULT_RS_FILE = \
    os.path.join(DEFAULT_DATA_DIR, DEFAULT_DATASET_NAME + ".rs.csv")

# Minimum number of samples from each cluster
DEFAULT_MIN_COVERAGE = 10

# Minimum number of samples to output
DEFAULT_MIN_OUTPUT_COUNT = 10000

###########################################
# Utility functions
###########################################

# Reads cluster to samples mapping and exemplar indices
# into internal list structures
# @param c2s_path - path to cluster to samples mapping
# @param ex_path - path to exemplar identifiers
# @return tuple with first item a c2s list of lists and
#   second item a list of exemplar identifiers
def load_data(c2s_path, ex_path):
    # Process cluster to samples mapping
    print "Loading " + c2s_path
    c2s = []
    with open(c2s_path) as f_c2s:
        for line in f_c2s:
            str_ids = line.split(",")
            int_ids = [int(i) for i in str_ids]
            c2s.append(int_ids)

    # Process exemplar ids
    print "Loading " + ex_path
    exs = []
    with open(ex_path) as f_ex:
        for line in f_ex:
            int_id = int(line)
            exs.append(int_id)

    return c2s, exs

# Writes resampled identifiers to disk
# @param resamples - list of resampled ids
# @param out_path - disk path to write data to
def write_data(resamples, out_path):
    print "Writing " + out_path
    with open(out_path, "w") as f:
        for i in resamples:
            f.write(str(i) + "\n")

# @param c2s - list of lists: clusters to samples mapping
# @param exs - list of exemplar identifiers
# @param min_cov - integer min coverage for each cluster
# @param min_count - minimum total resample count
# @return tuple with first item containing flattened list
#   of selected sample identifiers, second item is a list of lists
#   where each sublist contains the identifiers selected for each cluster
def resample(c2s, exs, min_cov, min_count):
    print "Resampling from " + str(len(exs)) + " clusters"
    if (min_cov < 1) :
        print "Min coverage must be at least 1."
        return [], []

    if (min_count < 1):
        print "Min count must be at least 1."
        return [], []

    # Flattened list of resampled identifiers
    resamples = []

    # List of list - each sublist stores identifiers that have already been selected from each cluster
    resamples_by_clust = []

    # Note - want to "deep" copy most of these lists
    # Sample from min coverage
    print "Sampling to obtain minimum coverage"
    for cid, clust in enumerate(c2s):
        resamples_by_clust.append([])
        # clust is a list of sample ids belonging to this cluster
        samples = []
        # If cluster is too small for min_cov, then just sample entire cluster
        if (min_cov >= len(clust)):
            samples.extend(clust)
            # Now keep sampling to make sure we have balanced representation
            while (len(samples) < min_cov):
                # Note: random.sample returns a list, we just want the first element (as this is a singleton list)
                s = (random.sample(clust, 1))[0]
                samples.append(s)
        else:
            # Else, select a random subsample
            samples.extend(random.sample(clust, min_cov))
            # And make sure exemplar is always included
            if (exs[cid] not in samples):
                samples[0] = exs[cid]
        # Add samples to output
        resamples.extend(samples)
        # Track which samples have been selected from cluster
        resamples_by_clust[cid].extend(samples)

    # Sample to obtain min count
    print "Sampling to obtain minimum count"
    n_clust = len(c2s)
    while (len(resamples) < min_count):
        # Select a cluster
        cid = random.randint(0, n_clust-1)
        # Now, select sample from cluster
        s = -1
        # Determine if unselected samples exist
        n_prev_selected = len(resamples_by_clust[cid])
        n_clust_size = len(c2s[cid])
        if (n_prev_selected >= n_clust_size):
            # All samples previously selected, just select any at random
            s = (random.sample(c2s[cid], 1))[0]
        else:
            # Some exist, only select from the unselected subpopulation
            unselected = list((set(c2s[cid]) - set(resamples_by_clust[cid])))
            s = (random.sample(unselected, 1))[0]
        resamples.append(s)
        resamples_by_clust[cid].append(s)

    # Finally, shuffle to locally decorrelate samples
    print "Performing final shuffle"
    random.shuffle(resamples)

    return resamples, resamples_by_clust

# Visual representation of how clusters were sampled - similar to silhouette plots
# @param resamples_by_clust - list of lists, each sublist contains sample identifiers
#   selected for that cluster
def render_sample_allocation(resamples_by_clust):
    print "Resampled allocations for each cluster:"
    for subpop in resamples_by_clust:
        alloc = ""
        for i in xrange(len(subpop)):
            alloc += "+"
        print alloc

###########################################
# Main
###########################################

# Main script entry point
def __main__():
    print "======================= Cluster Resampler ======================="

    # Seed random
    random.seed()

    # Initialize command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-c2s', '--path_to_c2s_file',
        default=DEFAULT_C2S_FILE,
        help='Path to cluster to sample mapping file')
    parser.add_argument('-ex', '--path_to_exemplar_file',
        default=DEFAULT_EX_FILE,
        help='Path to exemplar file')
    parser.add_argument('-rs', '--path_to_out_rs_file',
        default=DEFAULT_RS_FILE,
        help='Path to write out resampling')
    parser.add_argument('-mcov', '--min_coverage',
        default=DEFAULT_MIN_COVERAGE,
        help='Minimum number of samples to keep from each cluster')
    parser.add_argument('-mcout', '--min_output_count',
        default=DEFAULT_MIN_OUTPUT_COUNT,
        help='Minimum number of total samples to output')

    # Parse command line
    args = parser.parse_args()

    # Print command line
    print '\t-c2s = ' + args.path_to_c2s_file
    print '\t-ex = ' + args.path_to_exemplar_file
    print '\t-rs = ' + args.path_to_out_rs_file
    print '\t-mcov = ' + str(args.min_coverage)
    print '\t-mcout = ' + str(args.min_output_count)

    # Load data
    c2s, exs = load_data(args.path_to_c2s_file, args.path_to_exemplar_file)

    # Resample
    resamples, resamples_by_clust = resample(c2s, exs, int(args.min_coverage), int(args.min_output_count))

    # Visualize sampling allocation
    render_sample_allocation(resamples_by_clust)

    # Write re-sampled data
    write_data(resamples, args.path_to_out_rs_file)

    print "Resampling finished."

# Run main if we are the active script
if __name__ == '__main__':
    __main__()
