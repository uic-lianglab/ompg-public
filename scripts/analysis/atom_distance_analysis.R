# Purpose is to analyze atom distances for specific backbone atom
# ranges.

###################################################################
# Libraries
###################################################################

library(Rpdb)

###################################################################
# Common path utilities
###################################################################

# Cached directory to script
# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
SCRIPT_DIR = getSrcDirectory(function(x) {
  x
})

# @return full path to directory containing this script
get_script_dir <- function() {
  return(SCRIPT_DIR)
}

# @return root directory of project
get_root_dir <- function() {
  file.path(get_script_dir(), '..', '..')
}

###################################################################
# OMPG data paths
###################################################################

get_template_pdb_2iwv_fpath <- function() {
  file.path(get_root_dir(), 'ompg', 'templates', '2iwv.pdb')
}

###################################################################
# OMPG data
###################################################################

BB_ATOM_NAMES = c("N", "CA", "C", "O", "CB")

PDB_2IWV = read.pdb(file=get_template_pdb_2iwv_fpath(), HETATM=FALSE)

build_bb_atom_map <- function(pdb) {
  num_res = max(pdb$atoms$resid)
  bb_atom_map = matrix(nrow=num_res, ncol=length(BB_ATOM_NAMES))
  colnames(bb_atom_map) = BB_ATOM_NAMES
  num_rec = length(pdb$atoms$resid)
  for (i in 1:num_rec) {
    atom_name = pdb$atoms$elename[i]
    if (atom_name %in% BB_ATOM_NAMES) {
      res_ix = pdb$atoms$resid[i]
      bb_atom_map[res_ix, atom_name] = i 
    }
  }
  return(bb_atom_map)
}

BB_ATOM_MAP_2IWV = build_bb_atom_map(PDB_2IWV)

###################################################################
# Summary
###################################################################

print_summary_bb_dist <- function(bb_atom_name_1, bb_atom_name_2, offset, pdb, bb_atom_map) {
  stopifnot(bb_atom_name_1 %in% BB_ATOM_NAMES)
  stopifnot(bb_atom_name_2 %in% BB_ATOM_NAMES)
  num_rec = length(pdb$atoms$resid)
  dists = c()
  atom_name_2 = bb_atom_name_2
  for (atom_ix_1 in 1:num_rec) {
    atom_name_1 = pdb$atoms$elename[atom_ix_1]
    if (atom_name_1 == bb_atom_name_1) {
      res_ix_1 = pdb$atoms$resid[atom_ix_1]
      res_ix_2 = res_ix_1 + offset
      # Skip if no residue at offset
      if (res_ix_2 < 1) {
        next
      }
      if (res_ix_2 > nrow(bb_atom_map)) {
        next
      }
      atom_ix_2 = bb_atom_map[res_ix_2, atom_name_2]
      # Skip if atom does not exist (e.g. CB for GLY)
      if (is.na(atom_ix_2)) {
        next
      }
      # Calculate distance between the atoms
      atom_pos_1 = c(pdb$atoms$x1[atom_ix_1],
                     pdb$atoms$x2[atom_ix_1],
                     pdb$atoms$x3[atom_ix_1])
      atom_pos_2 = c(pdb$atoms$x1[atom_ix_2],
                     pdb$atoms$x2[atom_ix_2],
                     pdb$atoms$x3[atom_ix_2])
      dist2 = atom_pos_1 - atom_pos_2
      dist2 = dist2 * dist2
      dist2 = sum(dist2)
      dist = sqrt(dist2)
      dists = c(dists, dist)
    }
  }
  print(paste0("SUMMARY FOR ",
               bb_atom_name_1,
               ", ",
               bb_atom_name_2,
               " + ",
               offset,
               ":"))
  print(summary(dists))
}

print_summary_bb_dist('N', 'N', 1, PDB_2IWV, BB_ATOM_MAP_2IWV)
print_summary_bb_dist('C', 'C', 1, PDB_2IWV, BB_ATOM_MAP_2IWV)
print_summary_bb_dist('O', 'CA', 1, PDB_2IWV, BB_ATOM_MAP_2IWV)
