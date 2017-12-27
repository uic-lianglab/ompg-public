#!/usr/bin/python
# -*- coding: utf-8 -*-

# Script processes a folder of barrel renderings and outputs a csv file with
# columns:
#   NAME - name of PNG file
#   L1 - % of lumen pixels mapped to loop 1
#   L2 - % of lumen pixels mapped to loop 2
#   L3 - % of lumen pixels mapped to loop 3
#   L5 - % of lumen pixels mapped to loop 5
#   L6 - % of lumen pixels mapped to loop 6
#   L7 - % of lumen pixels mapped to loop 7
#   BG - % of lumen pixels napped to background (void)
#   WINL - Identifier of max lumen pixel coverage among loops
#   WINALL - Identifier of max lumen pixel coverage among loops + background

# Arguments:
#   <path_to_folder>: Folder contains PNG images to process
#   <out_csv_path>: Location to write summary csv file

###########################################
# Imports
###########################################

import os
import sys
import numpy as np
from scipy import misc
from scipy.spatial import distance

###########################################
#Globals
###########################################

# Path containing this script
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

# Command line arguments are in this order within sys.argv list
ARG_IN_PNG_DIR = 1
ARG_OUT_CSV_PATH = 2
CMD_LINE_SIZE = 3

# Project root directory
ROOT_DIR =  os.path.join(SCRIPT_DIR, "..", "..")
ROOT_DIR = os.path.abspath(ROOT_DIR)

# Assumed base directories
OUTPUT_DIR = os.path.join(ROOT_DIR, "ompg", "output")
RELAX_DIR = os.path.join(OUTPUT_DIR, "relax")
CAPT_BASE_DIR = os.path.join(RELAX_DIR, "capt")

# Defaults for testing
DEF_IN_PNG_DIR = os.path.join(CAPT_BASE_DIR, "pH5", "wt", "regsc",
                              "sids.merge.close.PNGs")
DEF_OUT_CSV_PATH = os.path.join(CAPT_BASE_DIR, "pH5", "wt", "regsc",
                                "sids.merge.close.occlusion.csv")

# Region labels
ID_L1 = 0
ID_L2 = 1
ID_L3 = 2
ID_L5 = 3
ID_L6 = 4
ID_L7 = 5
ID_BG = 6
ID_BARREL = 7

# RGB tags for each loop and barrel
RGB_L1 = np.array([191, 0, 191, 255], dtype=int)
RGB_L2 = np.array([0, 0, 255, 255], dtype=int)
RGB_L3 = np.array([0, 255, 0, 255], dtype=int)
RGB_L5 = np.array([50, 153, 50, 255], dtype=int)
RGB_L6 = np.array([255, 255, 0, 255], dtype=int)
RGB_L7 = np.array([255, 0, 0, 255], dtype=int)
RGB_BG = np.array([255, 255, 255, 255], dtype=int)
RGB_BARREL = np.array([0, 255, 255, 255], dtype=int)

# These arrays must be in order defined by ID_L1 to ID_BARREL
RGB_LST = [RGB_L1, RGB_L2, RGB_L3, RGB_L5, RGB_L6, RGB_L7, RGB_BG, RGB_BARREL]
RGB_LAB = [ID_L1, ID_L2, ID_L3, ID_L5, ID_L6, ID_L7, ID_BG, ID_BARREL]
RGB_STR = ["L1", "L2", "L3", "L5", "L6", "L7", "BG", "BARREL"]

#############################################################################
# Methods
#############################################################################


# @return 2-D matrix with region label for each pixel
def label_im(im):
    """Assigns each pixel to a region label."""
    nrow = im.shape[0]
    ncol = im.shape[1]
    lab = np.empty([nrow, ncol], dtype=int)
    
    for irow in range(0, nrow):
        for icol in range(0, ncol):
            elem = im[irow, icol, :]
            min_dist2 = distance.sqeuclidean(elem, RGB_LST[0])
            min_ix = 0
            for imin in range(1, len(RGB_LST)):
                dist2 = distance.sqeuclidean(elem, RGB_LST[imin])
                if (dist2 < min_dist2):
                    min_dist2 = dist2
                    min_ix = imin
            lab[irow, icol] = RGB_LAB[min_ix]
    
    return lab


def get_row_spans(row, label, min_len):
    """@return [(first pix, last+1 pix), ...] for regions with target label."""
    STATE_SEARCH=0
    STATE_FOUND=1    
    state = STATE_SEARCH
    out_spans = []
    span = [0,0]
    e = -1
    for i in range(0, len(row)):
        e = row[i]
        if (e == label) and (state == STATE_SEARCH):
            span[0] = i
            state = STATE_FOUND
        elif (e != label) and (state == STATE_FOUND):
            span[1] = i
            if (span[1] - span[0] + 1) >= min_len:
                out_spans.append((span[0], span[1]))
            state = STATE_SEARCH
                 
    if (e == label) and (state == STATE_FOUND):
        span[1] = len(row)
        if (span[1] - span[0] + 1) >= min_len:
            out_spans.append((span[0], span[1]))

    return out_spans


def get_row_occ_count(lab_row):
    """Determines occlusion counts for labeled row."""
    # Define minimum number of pixels span lengths
    # in order to filter noise:
    # Minimum number of pixels defining a barrel span
    MIN_BARREL_SPAN = 15
    # Minimum number of pixels defining a lumen span
    MIN_LUMEN_SPAN = 10
    occ = np.zeros(len(RGB_LAB))
    spans = get_row_spans(lab_row, ID_BARREL, MIN_BARREL_SPAN)
    if (len(spans) == 2):
        lum_start = spans[0][1]
        lum_end = spans[1][0]
        if (lum_end - lum_start) >= MIN_LUMEN_SPAN:
            for i in range(lum_start, lum_end):
                occ[lab_row[i]] = occ[lab_row[i]] + 1
    return occ
        

def get_occ(lab):
    """Determines occlusion profile for labeled image."""
    occ = np.zeros(len(RGB_LAB)) 
    for row in lab:
        occ = np.add(occ, get_row_occ_count(row))
    # Crop barrel pixel counts - they can creep in due to
    # minor spans being ignored and counted as lumen
    occ = occ[0:(len(occ)-1)]
    occ = occ / np.sum(occ)
    return occ


def occlusion_analysis(in_png_dir, out_csv_path):
    """Determines occlusion % for each loop."""
    # Process each image
    with open(out_csv_path, 'w') as f:
        f.write("NAME,L1,L2,L3,L5,L6,L7,NONE,WINL,WINALL\n")
        for png_name in os.listdir(in_png_dir):
            png_path = os.path.join(in_png_dir, png_name)
            # Assume each file is a PNG
            if os.path.isfile(png_path):
                print "Processing: " + png_name
                png = misc.imread(png_path)
                lab = label_im(png)
                occ = get_occ(lab)
                id_winl = np.argmax(occ[0:len(occ)-1])
                id_winall = np.argmax(occ)
                line = png_name + "," + str(occ[ID_L1]) + "," + \
                    str(occ[ID_L2]) + "," + str(occ[ID_L3]) + "," + \
                    str(occ[ID_L5]) + "," + str(occ[ID_L6]) + "," + \
                    str(occ[ID_L7]) + "," + str(occ[ID_BG]) + "," + \
                    RGB_STR[id_winl] + "," + RGB_STR[id_winall] + "\n"
                f.write(line)


#############################################################################
# Main
#############################################################################


def __main__():
    
    in_png_dir = DEF_IN_PNG_DIR
    out_csv_path = DEF_OUT_CSV_PATH

    # Check if user arguments exist
    arg_count = len(sys.argv)
    if (arg_count == CMD_LINE_SIZE):
        in_png_dir = sys.argv[ARG_IN_PNG_DIR]
        out_csv_path = sys.argv[ARG_OUT_CSV_PATH]
        
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "OCCLUSION ANALYSIS"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "-in_png_dir: " + in_png_dir
    print "-out_csv_path: " + out_csv_path
    
    occlusion_analysis(in_png_dir, out_csv_path)

    print "Occlusion analysis script finished."


if __name__ == '__main__':
    __main__()
