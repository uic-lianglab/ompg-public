#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Remove duplicate files from a folder."""

#############################################################################
# Uses sha1 and md5 hashes to filter duplicates within a folder.
#
# Usage:
#   python <path_to_this_script> <path_to_folder> <path_to_output_map> <0|1>
#   - <path_to_this_script>: path to this .py script to run with python
#       interpreter
#   - <path_to_folder>: path to folder to remove duplicates for
#   - <path_to_output_map>: path to write csv entries of form:
#           |file name|,|identical file name which was retained|
#       In other words, the first column lists the file name and the second
#       column lists the file that can serve as a template for that file
#       because they are identical. If set to delete, then only the template
#       files listed in second column are retained.
#   - <0|1>: 0 -> do not delete duplicates, 1 -> delete duplicates
#
# Some relevant discussions:
# http://superuser.com/questions/556923/whats-the-percentage-of-having-two-files-with-the-same-byte-size-length-giving
# http://security.stackexchange.com/questions/5229/mathematically-theoretically-what-is-the-chance-that-2-different-inputs-would
#############################################################################

#############################################################################
# Imports
#############################################################################

import csv
import hashlib
import os
import sys

#############################################################################
# Globals
#############################################################################

# Path to script directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Enumerated arguments
ARG_DUP_FOLDER_PATH = 1  # Path to folder to remove duplicates for
ARG_OUT_MAP_PATH = 2  # Path for writing file to unique template mapping
ARG_B_DELETE = 3  # 0 -> leave duplicates untouched, 1 -> delete them!
CMD_LINE_SIZE = 4  # Expected command line length

# Defaults
DEF_BLOCK_SIZE = 65536
DEF_BASE_PATH = os.path.join(SCRIPT_DIR, "..", "ompg", "output",
                             "loop_6_r228d_2iww", "pH7", "prot_raw")
DEF_DUP_FOLDER_PATH = os.path.join(DEF_BASE_PATH, "psf")
DEF_OUT_MAP_PATH = os.path.join(DEF_BASE_PATH, "psf_templates.csv")
DEF_B_DELETE = 1

#############################################################################
# Methods
#############################################################################


# @param fpath - path of file to get md5, sha1, and file size of
# @param block_size - blocks of this size bytes are read in at a time
# @return tuple (md5 hash, sha1 hash, file size)
def get_file_id(fpath, block_size=DEF_BLOCK_SIZE):
    """Get sha1, md5, and file size attributes."""
    f = open(fpath, 'rb')
    h_md5 = hashlib.md5()
    h_sha1 = hashlib.sha1()
    buf = f.read(block_size)
    while len(buf) > 0:
        h_md5.update(buf)
        h_sha1.update(buf)
        buf = f.read(block_size)
    f.close()
    return (h_md5.hexdigest(), h_sha1.hexdigest(), os.path.getsize(fpath))


# Based on code from:
# http://pythoncentral.io/finding-duplicate-files-with-python/
# @param dup_folder_path - path to folder to search for duplicates
# @param out_map_path - path to write mapping from duplicates to unique files
# @param b_delete - 0 -> do not modify duplicates, 1 -> delete duplicates
def remove_dups(dup_folder_path, out_map_path, b_delete):
    """Search for duplicates within dup_folder_path."""
    # Current number of files processed
    num_files = 0
    # Current number of duplicates
    num_dups = 0
    # List of (file name, template file name) tuples
    file_to_template = []
    # Mapping from file identifier to first encountered file index with that id
    fids = {}
    # Process folder
    print "Processing folder: " + dup_folder_path
    for fname in os.listdir(dup_folder_path):
        print str(num_files) + ": Reading " + fname
        template_fname = fname
        fpath = os.path.join(dup_folder_path, fname)
        fid = get_file_id(fpath)
        if fid in fids:
            num_dups = num_dups + 1
            template_ix = fids[fid]
            template_fname = file_to_template[template_ix][0]
            print "\tDuplicate of " + str(template_ix) + ": " + template_fname
            if (b_delete):
                os.remove(fpath)
        else:
            fids[fid] = num_files
        file_to_template.append((fname, template_fname))
        num_files = num_files + 1
    # Print summary
    perc_dups = 100.0 * float(num_dups) / float(max(num_files, 1))
    print "Found " + str(num_dups) + " of " + str(num_files) + " = " + \
        str(perc_dups) + "% duplicate files."
    # Write file to template map
    if (num_files > 0):
        print "Writing duplicate mapping to: " + out_map_path
        with open(out_map_path, "w") as f:
            writer = csv.writer(f)
            writer.writerows(file_to_template)

#############################################################################
# Main
#############################################################################


def __main__():
    # Initialize with defaults
    dup_folder_path = DEF_DUP_FOLDER_PATH
    out_map_path = DEF_OUT_MAP_PATH
    b_delete = DEF_B_DELETE

    # Check command line
    if (len(sys.argv) == CMD_LINE_SIZE):
        dup_folder_path = sys.argv[ARG_DUP_FOLDER_PATH]
        out_map_path = sys.argv[ARG_OUT_MAP_PATH]
        b_delete = int(sys.argv[ARG_B_DELETE])
    else:
        print "WARNING, expected " + str(CMD_LINE_SIZE) + \
            " arguments but only " + str(len(sys.argv)) + " read."
        print "Ignoring command line and using default values."

    # Validate arguments
    if (not os.path.exists(dup_folder_path)):
        print "ERROR, duplicate folder path:\n\t" + dup_folder_path + \
            "\ndoes not exist. Exiting."
        sys.exit(0)
    if ((b_delete != 0) and (b_delete != 1)):
        print "ERROR, b_delete argument " + str(b_delete)
        + " is not equal to 0 or 1. Exiting."
        sys.exit(0)

    # Print parameters
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "REMOVING DUPLICATES:"
    print "\tfolder to check: " + dup_folder_path
    print "\toutput map: " + out_map_path
    print "\tb_delete: " + str(b_delete)

    # Remove duplicates!
    remove_dups(dup_folder_path,
                out_map_path,
                b_delete)

if __name__ == '__main__':
    __main__()
