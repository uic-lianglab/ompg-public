#!/usr/bin/python
# -*- coding: utf-8 -*-

# THIS VERSION HAS HARDCODED MEMBRANE Z-BILAYER COORDS FOR 2iww AND 2iwv.
# THIS MAKES IT EASIER TO PROCESS DIRECTORIES CONTAINING MIXED LOOPS FROM BOTH
# OF THESE TEMPLATES.

# Scripts main purpose is to conclude if an oriented,
# barrel membrane protein is in a closed or open
# conformation based on castp output

# Classfier is based on if there's a pocket
# that has at least 2 subpockets classified as mouths.
# Furthermore, 1 mouth/subpocket must be primarily
# above the outer bilayer and the other mouth/subpocket
# must be primarily below the inner bilayer

import argparse
import os
from collections import defaultdict
import sys


# OMPs are assumed oriented along the z-axis
# In other words, z is normal to the plane of the bilayer

# Bilayer z values
IX_OMP_BILAYER_IN_Z = 0
IX_OMP_BILAYER_OUT_Z = 1
OMP_BILAYER_Z = {"2iwv": (-12.650, 12.650), "2iww": (-12.35, 12.35)}

# Obtain path to script directory
SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))

DEFAULT_OMP_BASE_DIR = os.path.join(SCRIPT_DIR, '..', 'ompg', 'output', 'merge', 'castp', 'pH7', 'wt')

# Default folder for outer membrane (omp) pdb data
DEFAULT_OMP_IN_PDB_DATA_DIR = os.path.join(DEFAULT_OMP_BASE_DIR, 'pdb')

# Default folder for outer membrane (omp) castp data
DEFAULT_OMP_IN_CASTP_DATA_DIR = os.path.join(DEFAULT_OMP_BASE_DIR, '2.75')

# Default path to output results
DEFAULT_OMP_OUT_RESULTS_PATH = os.path.join(DEFAULT_OMP_BASE_DIR, "oc.results.csv")

# Default minimum number of mouth atoms above the bilayer
# for a pocket to classified as spanning the outer bilayer
DEFAULT_MIN_MOUTH_ATOMS_ABOVE_BILAYER = 2

# Default minimum number of mouth atoms below the bilayer
# for a pocket to be classified as spanning the inner bilayer
DEFAULT_MIN_MOUTH_ATOMS_BELOW_BILAYER = 2

# For a pocket with mouths on both sides of the bilayer,
# default minimum solvent accessible surface area
# necessary for the protein to be considered 'open'
DEFAULT_MIN_SA_AREA = 50.0 # small value for now

# Return value for a closed protein
RV_OMP_CLOSED = 0

# Return value for an open protein
RV_OMP_OPEN = 1


# A 3-D atom position
class PDBAtom:
    
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


# Utility class for obtaining paths to all necessary
# files (PDB, CASTp) necessary for classification
class CASTpPathManager:

    def __init__(self, in_castp_data_dir, in_pdb_data_dir):
        self.in_castp_data_dir = in_castp_data_dir
        self.in_pdb_data_dir = in_pdb_data_dir

    # Based on identifier, get the pdb file name
    @staticmethod
    def to_pdb_fname(omp_id):
        return omp_id + '.pdb'

    # Based on the identifier, get the mouth info file name
    # containing the meta/summary information
    # which maps pockets to mouths to areas
    @staticmethod
    def to_mouth_info_fname(omp_id):
        return omp_id + '.mouthInfo'

    # Based on the identifier, get the mouth file name
    # containing positions for each atom belonging to a mouth
    @staticmethod
    def to_mouth_coords_fname(omp_id):
        return omp_id + '.mouth'

    # Based on the identifier, get the subpocket info file name
    # containing meta data specifying the
    # type of subpocket (1=mouth, 2=neck, 3=body)
    @staticmethod
    def to_subpocket_info_fname(omp_id):
        return omp_id + '.subpocket'

    # Based on the identifier, get the subpocket file name
    # containing positions for each atom of the subpocket
    @staticmethod
    def to_subpocket_coords_fname(omp_id):
        return omp_id + '.subpoc'

    # Based on the identifier, get the pocket info file name
    # containing meta data specifying the
    # areas, lengths, and number of mouths
    @staticmethod
    def to_pocket_info_fname(omp_id):
        return omp_id + '.pocInfo'

    # Based on the identifier, get the pocket file name
    # containing positions for each atom of the pocket
    @staticmethod
    def to_pocket_coords_fname(omp_id):
        return omp_id + '.poc'

    # Internal method to convert a file name to a full path
    # to castp data
    def to_castp_path_(self, fname):
        return os.path.join(self.in_castp_data_dir, fname)

    # Internal method to convert a file name to a full path
    # to pdb data
    def to_pdb_path_(self, fname):
        return os.path.join(self.in_pdb_data_dir, fname)

    def get_pdb_path(self, omp_id):
        return self.to_pdb_path_(CASTpPathManager.to_pdb_fname(omp_id))

    def get_mouth_info_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_mouth_info_fname(omp_id))

    def get_mouth_coords_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_mouth_coords_fname(omp_id))

    def get_subpocket_info_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_subpocket_info_fname(omp_id))

    def get_subpocket_coords_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_subpocket_coords_fname(omp_id))

    def get_pocket_info_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_pocket_info_fname(omp_id))

    def get_pocket_coords_path(self, omp_id):
        return self.to_castp_path_(CASTpPathManager.to_pocket_coords_fname(omp_id))

    def pocket_and_mouth_info_paths_exist(self, omp_id):
        return os.path.exists(self.get_pocket_info_path(omp_id)) and \
            os.path.exists(self.get_mouth_info_path(omp_id))


# Struct for storing mouth info fields
class CASTpMouthInfo:

    def __init__(self):
        self.molecule = ''
        self.poc_id = '' # Pocket containing mouths
        self.num_mouths = 0 # Number of mouths in pocket
        self.area_sa = 0.0 # Solvent accessible area
        self.area_ms = 0.0 # Molecular surface area
        self.len_sa = 0.0 # Length of solvent accessible surface
        self.len_ms = 0.0 # Length of molecular surface
        self.ntri = 0 # Number of triangles


# Struct for storing pocket info fields
class CASTpPocketInfo:

    def __init__(self):
        self.molecule = ''
        self.poc_id = '' # Pocket containing mouths
        self.num_mouths = 0 # Number of mouths in pocket
        self.area_sa = 0.0 # Solvent accessible area
        self.area_ms = 0.0 # Molecular surface area
        self.vol_sa = 0.0 # Volume of solvent accessible surface
        self.vol_ms = 0.0 # Volume of molecular surface
        self.length = 0.0 # Length of pocket
        self.cnr = 0


# Utility class for parsing a line from a mouth|subpocket coords file
class CASTpCoordsAtomLineParser:

    def __init__(self):
        self.rec = ''  # Record name (ATOM)
        self.serial = ''  # Atom serial number (Integer)
        self.atom = ''  # Atom name (Atom)
        self.alt_loc = ''  # alternate location indicator (Character)
        self.res = ''  # Residue name
        self.chain = ''  # Chain identifier (Character)
        self.res_num = ''  # Residue sequence number (Integer)
        self.ins = ''  # Code for insertion of residues (AChar)
        self.x = ''  # Orthogonal coords for X in Angstroms (Real-8.3)
        self.y = ''  # Orthogonal coords for Y in Angstroms (Real-8.3)
        self.z = ''  # Orthogonal coords for Z in Angstroms (Real-8.3)
        self.occ = ''  # Occupancy (Real-6.2)
        self.temp = ''  # Temperature factor (Default = 0.0) (Real-6.2)
        self.poc_id = '' # Format: <PocketId>|<PocketId.SubpocketId>
        self.castp_rec = '' # Should always be 'POC'

    def parse(self, line):
        self.rec = line[0:6].strip()
        self.serial = line[6:11].strip()
        self.atom = line[12:16].strip()
        self.alt_loc = line[16]
        self.res = line[17:20].strip()
        self.chain = line[21].strip()
        self.res_num = line[22:26].strip()
        self.ins = line[26]
        self.x = line[30:38].strip()
        self.y = line[38:46].strip()
        self.z = line[46:54].strip()
        self.occ = line[54:60].strip()
        self.temp = line[60:66].strip()
        # Note, seg, elem, and charge are missing from
        # reference pdbs 2iwv and 2iww, so if seg,
        # elem, and charge exist in the pdbs, then may
        # need to include them as well

        # HACK: tokenizing by whitespace and then
        # extracting last two tokens.  This is so
        # this class can be reused for both mouth and
        # subpocket position files
        tmp = line.split()
        assert(len(tmp) >= 2)
        self.poc_id = tmp[-2]
        self.castp_rec = tmp[-1]
        # sanity check
        assert (self.castp_rec == 'POC') or (self.castp_rec == 'M4P')


# Stores the positions for atoms it encounters
class CASTpStoreCoordsVisitor:
    def __init__(self):
        self.coords = defaultdict(list)

    def visit(self, parser):
        atom = PDBAtom()
        atom.x = float(parser.x)
        atom.y = float(parser.y)
        atom.z = float(parser.z)
        self.coords[parser.poc_id].append(atom)


# Used for sanity checking the subpocket
# atoms with the mouth atoms to make sure
# they are the same (subpocket contains all)
# mouth atoms
class CASTpStoreAtomSerialVisitor:
    def __init__(self):
        self.atom_serials = defaultdict(list)

    def visit(self, parser):
        self.atom_serials[parser.poc_id].append(parser.serial)


# Function parses a subpockets coordinates file and
# returns a dictionary where the key is a subpocket id
# and the value is a list of atom positions.  Note that
# subpocket ids are in the form <pocket_id.subpocket_id>
def visit_castp_coords(path_to_coords, visitor):
    fo_in = open(path_to_coords, 'r')
    parser = CASTpCoordsAtomLineParser()
    for line in fo_in:
        parser.parse(line)
        visitor.visit(parser)
    fo_in.close()


# @param in_data_dir string of name of directory
#   containing OMP pbds and CASTp data
def get_omp_ids_from_dir(in_pdb_data_dir):
    out_omp_ids = []
    if os.path.exists(in_pdb_data_dir):
        for file in os.listdir(in_pdb_data_dir):
            if file.endswith(".pdb") and not file.endswith(".poc.pdb"):
                out_omp_ids.append( os.path.splitext(file)[0] )
    return out_omp_ids


# Function determines which subpockets have mouths
# and returns a dictionary where the key is a pocket id
# and the value is a list of subpocket ids.  Note that
# subpocket ids are in the form <pocket_id.subpocket_id>
def extract_subpockets_with_mouths(path_to_subpocket_info):

    # Enumerated column indices (there are more but we only need these)
    IX_SUBPOC = 0 # sub-pocket id
    IX_TYPE = 1 # 1-mouth, 2-neck and 3-body

    # The identifier for a mouth
    ID_MOUTH = "1"

    # The number of columns we expect to encounter
    # for a line that we are interested in
    EXPECTED_COLUMNS = 6

    out_subpockets = defaultdict(list)
    fo_in = open(path_to_subpocket_info, 'r')
    for line in fo_in:
        # Remove leading and trailing whitespace
        line.strip()

        # Ignore comment lines
        if line.startswith("#"):
            continue

        # Split based on ';' and ignore any line that
        # doesn't have expected number of entries
        tokens = line.split(";")
        if len(tokens) != EXPECTED_COLUMNS:
            continue

        # See if we are parsing a subpocket labeled as a mouth
        if tokens[IX_TYPE] == ID_MOUTH:
            # The full subpocket id is prefixed with the pocket id (e.g 28.2)
            full_subpoc_id = tokens[IX_SUBPOC]
            # Determine pocket id by splitting on '.'
            result = full_subpoc_id.split(".")
            # Make sure id is of form <poc_id>.<subpoc_id>
            if (len(result) != 2):
                print "WARNING: unexpected token for line:\n\t" + \
                    line + "\n\tin file: " + path_to_subpocket_info
                continue
            # Append full subpocket id
            poc_id = result[0]
            out_subpockets[poc_id].append(full_subpoc_id)
    fo_in.close()
    return out_subpockets


# Function parses a mouth info file and
# returns a dictionary mapping the pocket
# identifier to the mouth statistics
def extract_mouth_info(path_to_info):
    # Enumerated column indices
    IX_MI_PREFIX = 0 # Line starts with 'MTH:'
    IX_MI_MOLECULE = 1 # Name of pdb
    IX_MI_POC_ID = 2 # Pocket identifier
    IX_MI_N_MTH = 3 # Number of mouths
    IX_MI_AREA_SA = 4 # Solvent accessible area
    IX_MI_AREA_MS = 5 # Molecular surface area
    IX_MI_LEN_SA = 6
    IX_MI_LEN_MS = 7 # Length of molecular surface
    IX_MI_NTRI = 8 # Number of triangles

    # The number of columns we expect to encounter
    # for a line that we are interested in
    EXPECTED_COLUMNS = 9

    out_infos = {}
    fo_in = open(path_to_info, 'r')
    # Skip past first line
    next(fo_in)
    for line in fo_in:
        # Remove leading and trailing whitespace
        line.strip()
        tokens = line.split()
        assert(len(tokens) == EXPECTED_COLUMNS)
        assert(tokens[IX_MI_PREFIX] == 'MTH:')

        info = CASTpMouthInfo()
        info.molecule = tokens[IX_MI_MOLECULE]
        info.poc_id = tokens[IX_MI_POC_ID]
        info.num_mouths = int(tokens[IX_MI_N_MTH])
        info.area_sa = float(tokens[IX_MI_AREA_SA])
        info.area_ms = float(tokens[IX_MI_AREA_MS])
        info.len_sa = float(tokens[IX_MI_LEN_SA])
        info.len_ms = float(tokens[IX_MI_LEN_MS])
        info.ntri = int(tokens[IX_MI_NTRI])

        assert info.poc_id not in out_infos
        out_infos[info.poc_id] = info

    fo_in.close()
    return out_infos


# Function parses a pocket info file and
# returns a dictionary mapping the pocket
# identifier to the pocket statistics
def extract_pocket_info(path_to_info):
    # Enumerated column indices
    IX_PI_PREFIX = 0 # Line starts with 'POC:'
    IX_PI_MOLECULE = 1 # Name of pdb
    IX_PI_POC_ID = 2 # Pocket identifier
    IX_PI_N_MTH = 3 # Number of mouths
    IX_PI_AREA_SA = 4 # Solvent accessible area
    IX_PI_AREA_MS = 5 # Molecular surface area
    IX_PI_VOL_SA = 6 # Volume of solvent accessible surface
    IX_PI_VOL_MS = 7 # Volume of molecular surface
    IX_PI_LENGTH = 8 # Length of pocket
    IX_PI_CNR = 9

    # The number of columns we expect to encounter
    # for a line that we are interested in
    EXPECTED_COLUMNS = 10

    out_infos = {}
    fo_in = open(path_to_info, 'r')
    # Skip past first line
    next(fo_in)
    for line in fo_in:
        # Remove leading and trailing whitespace
        line.strip()
        tokens = line.split()
        assert(len(tokens) == EXPECTED_COLUMNS)
        assert(tokens[IX_PI_PREFIX] == 'POC:')

        info = CASTpMouthInfo()
        info.molecule = tokens[IX_PI_MOLECULE]
        info.poc_id = tokens[IX_PI_POC_ID]
        info.num_mouths = int(tokens[IX_PI_N_MTH])
        info.area_sa = float(tokens[IX_PI_AREA_SA])
        info.area_ms = float(tokens[IX_PI_AREA_MS])
        info.vol_sa = float(tokens[IX_PI_VOL_SA])
        info.vol_ms = float(tokens[IX_PI_VOL_MS])
        info.length = float(tokens[IX_PI_LENGTH])
        info.cnr = int(tokens[IX_PI_CNR])

        assert info.poc_id not in out_infos
        out_infos[info.poc_id] = info

    fo_in.close()
    return out_infos


# Checks to make sure the poc_id keys match
def validate_pocket_and_mouth_info(poc_info, mouth_info):
    pk = set(poc_info.keys()) #set(poc_info.viewkeys())
    mk = set(mouth_info.keys()) #set(mouth_info.viewkeys())
    i = mk & pk
    assert i == pk
    assert i == mk
    for poc_id in mouth_info:
        assert poc_id in poc_info
        assert mouth_info[poc_id].num_mouths == \
            poc_info[poc_id].num_mouths


# Utility function to verify that subpockets
# contain all the mouth atoms defined in the
# mouth position/coords file
def validate_subpocket_mouths(omp_id, subpockets_with_mouths, path_manager):
    # Get subpocket atom serials
    subpocket_visitor = CASTpStoreAtomSerialVisitor()
    visit_castp_coords(path_manager.get_subpocket_coords_path(omp_id),
                       subpocket_visitor)
    subpocket_to_atoms_map = subpocket_visitor.atom_serials

    # Get mouth positions
    mouth_visitor = CASTpStoreAtomSerialVisitor()
    visit_castp_coords(path_manager.get_mouth_coords_path(omp_id),
                       mouth_visitor)
    mouth_pocket_to_atoms_map = mouth_visitor.atom_serials

    for mouth_poc_id, mouth_serials_list in mouth_pocket_to_atoms_map.items():
        # See if pocket id has been identified as containing a mouth
        if not (mouth_poc_id in subpockets_with_mouths):
            print "Warning: mouth pocket " + mouth_poc_id + " does not have associated subpockets"
            continue
        # Obtain subpockets ids that contain mouths
        subpoc_ids = subpockets_with_mouths[mouth_poc_id]
        # Iterate over each mouth atom and verify it is
        # found in a subpocket
        for mouth_serial in mouth_serials_list:
            atom_found = False
            # Check the subpockets for the atom
            for subpoc_id in subpoc_ids:
                subpocket_serials_list = subpocket_to_atoms_map[subpoc_id]
                atom_found = mouth_serial in subpocket_serials_list
                if (atom_found):
                    break
            if not atom_found:
                print "Warning: mouth atom " + mouth_serial + \
                    " not found for pocket " + mouth_poc_id


class OMPClosedClassifierApp:

    @staticmethod
    def append_cmd_args(parser):
        parser.add_argument('-omc_ic', '--omc_in_castp_data_dir',
                            default=DEFAULT_OMP_IN_CASTP_DATA_DIR,
                            help='Path to directory containing OMP CASTp data')

        parser.add_argument('-omc_ip', '--omc_in_pdb_data_dir',
                            default=DEFAULT_OMP_IN_PDB_DATA_DIR,
                            help='Path to directory containing OMP PDB data')

        parser.add_argument('-omc_op', '--omc_out_results_path',
                            default=DEFAULT_OMP_OUT_RESULTS_PATH,
                            help='Path for output results')

        parser.add_argument('-omc_min_o', '--omc_min_atoms_above',
                            default=DEFAULT_MIN_MOUTH_ATOMS_ABOVE_BILAYER,
                            help='min number of mouth atoms above bilayer necessary to classify a pocket as spanning outer bilayer')

        parser.add_argument('-omc_min_i', '--omc_min_atoms_below',
                            default=DEFAULT_MIN_MOUTH_ATOMS_BELOW_BILAYER,
                            help='min number of mouth atoms below bilayer necessary to classify a pocket as spanning inner bilayer')

        parser.add_argument('-omc_min_sa', '--omc_min_sa',
                            default=DEFAULT_MIN_SA_AREA,
                            help='min solvent accessible surface area to classify a pocket spannign the bilayer as open')


    @staticmethod
    def print_cmd_args(args):
        print '======================= OMP closed classifier ======================='
        print '\t-omc_ic = ' + args.omc_in_castp_data_dir
        print '\t-omc_ip = ' + args.omc_in_pdb_data_dir
        print '\t-omc_op = ' + args.omc_out_results_path
        print '\t-omc_min_o = ' + str(args.omc_min_atoms_above)
        print '\t-omc_min_i = ' + str(args.omc_min_atoms_below)
        print '\t-omc_min_sa = ' + str(args.omc_min_sa)


    @staticmethod
    def classify_omp(omp_id, out_memb_z_coord, in_memb_z_coord, min_atoms_above, min_atoms_below, min_sa, path_manager):
        # Find a pocket that
        # - spans the bilayer
        # - contains 1 mouth above bilayer
        # - contains 1 mouth below bilayer

        # If necessary pocket/mouth file is missing,
        # It's assumed that CASTp found no pockets
        # Defaults to returning closed but outputs warning
        if not path_manager.pocket_and_mouth_info_paths_exist(omp_id):
            print "Warning, unable to find one or more of the following paths:\n" + \
                "\tPocket information at: " + path_manager.get_pocket_info_path(omp_id) + "\n" + \
                "\tMouth information at: " + path_manager.get_mouth_info_path(omp_id) + "\n" + \
                "\nAssuming no pocket exists and therefore structure is CLOSED."
            return RV_OMP_CLOSED

        # Obtain information about pocket statistics
        poc_infos = \
            extract_pocket_info(path_manager.get_pocket_info_path(omp_id))

        # Obtain information about mouth statistics
        mouth_infos = \
            extract_mouth_info(path_manager.get_mouth_info_path(omp_id))

        # Sanity check pocket and mouth info
        validate_pocket_and_mouth_info(poc_infos, mouth_infos)

        # Filter out any pockets that have less than two mouths
        for poc_id, info in poc_infos.items():
            if info.num_mouths < 2:
                del poc_infos[poc_id]

        # If we have no viable pockets, return that we are closed
        if len(poc_infos.keys()) == 0:
            return RV_OMP_CLOSED

        # Get mouth positions
        mouth_visitor = CASTpStoreCoordsVisitor()
        visit_castp_coords(path_manager.get_mouth_coords_path(omp_id),
                           mouth_visitor)
        mouth_atoms = mouth_visitor.coords

        # Filter pockets that don't have mouth atoms
        # both above and below the bilayer
        for poc_id, poc_info in poc_infos.items():
            assert poc_info.num_mouths >= 2
            num_mouth_atoms_above_outer_bilayer = 0
            num_mouth_atoms_below_inner_bilayer = 0
            assert poc_id in mouth_atoms
            pocket_mouth_atoms = mouth_atoms[poc_id]
            for atom in pocket_mouth_atoms:
                num_mouth_atoms_above_outer_bilayer += int(atom.z > out_memb_z_coord)
                num_mouth_atoms_below_inner_bilayer += int(atom.z < in_memb_z_coord)

            spans_bilayer = num_mouth_atoms_above_outer_bilayer > min_atoms_above and \
                num_mouth_atoms_below_inner_bilayer > min_atoms_below

            if not spans_bilayer:
                del poc_infos[poc_id]

        # Check if the mouths spanning the bilayer meet the minimum
        # amount of solvent accessible surface area
        if (len(poc_infos.keys()) > 0):
            for poc_id in poc_infos:
                if mouth_infos[poc_id].area_sa >= min_sa:
                    return RV_OMP_OPEN

        return RV_OMP_CLOSED


    @staticmethod
    def classify_omps(in_castp_data_dir, in_pdb_data_dir, out_results_path, min_atoms_above, min_atoms_below, min_sa):
        omp_ids = get_omp_ids_from_dir(in_pdb_data_dir)

        if (omp_ids):
            path_manager = CASTpPathManager(in_castp_data_dir, in_pdb_data_dir)

            f = open(out_results_path, 'w')

            for omp_id in omp_ids:

                # Determine bilayer z values based on template key
                # Default assuming "2iwv" template
                in_memb_z_coord = OMP_BILAYER_Z["2iwv"][IX_OMP_BILAYER_IN_Z]
                out_memb_z_coord = OMP_BILAYER_Z["2iwv"][IX_OMP_BILAYER_OUT_Z]
                # Check if we should switch templates
                if "2iww" in omp_id:
                    in_memb_z_coord = OMP_BILAYER_Z["2iww"][IX_OMP_BILAYER_IN_Z]
                    out_memb_z_coord = OMP_BILAYER_Z["2iww"][IX_OMP_BILAYER_OUT_Z]
                elif "2iwv" not in omp_id:
                    print "Warning: unable to match template 2iwv or 2iww for " + \
                        omp_id + "\n\tDefaulting to 2iwv."

                rv = OMPClosedClassifierApp.classify_omp(
                    omp_id,
                    out_memb_z_coord,
                    in_memb_z_coord,
                    min_atoms_above,
                    min_atoms_below,
                    min_sa,
                    path_manager)

                print "(" + omp_id + ", " + str(rv) + ")"
                f.write(omp_id+","+str(rv)+"\n")

            f.close();
        print "Finished classifying " + str(len(omp_ids)) + " pdbs."


    @staticmethod
    def main(args):
        OMPClosedClassifierApp.classify_omps(
            args.omc_in_castp_data_dir,
            args.omc_in_pdb_data_dir,
            args.omc_out_results_path,
            int(args.omc_min_atoms_above),
            int(args.omc_min_atoms_below),
            float(args.omc_min_sa))


def __main__():
    app = OMPClosedClassifierApp()
    parser = argparse.ArgumentParser()
    app.append_cmd_args(parser)
    args = parser.parse_args()
    app.print_cmd_args(args)
    app.main(args)


if __name__ == '__main__':
    __main__()
