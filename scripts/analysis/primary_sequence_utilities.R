# Utilities for the extracting and displaying regions
# of interest within the OmpG primary sequence

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
# Paths
###################################################################

# @return path to template PDBs directory
get_templates_pdb_dir <- function() {
  file.path(get_root_dir(), "ompg", "templates")
}

# @return path to template PDBs with membranes directory
get_membranes_pdb_dir <- function() {
  file.path(get_templates_pdb_dir(), "membranes")
}

# @return PDB identifier for open conformation
open_id <- function() {
  return("2iwv")
}

# PDB identifer for closed conformation
close_id <- function() {
  return("2iww")
}

# @return path to membrane PDB file
get_membrane_pdb_path <- function(pdb_id = open_id()) {
  file.path(get_membranes_pdb_dir(), paste0(pdb_id, ".pdb"))
}

###################################################################
# Amino Acid Code Mappings
###################################################################

# Amino Acid 1 letter to 3 letter code mapping
# Example usage: AA1to3["A"], will output "ALA"
AA1to3 = c(
  "ALA",
  "ARG",
  "ASN",
  "ASP",
  "CYS",
  "GLU",
  "GLN",
  "GLY",
  "HIS",
  "ILE",
  "LEU",
  "LYS",
  "MET",
  "PHE",
  "PRO",
  "SER",
  "THR",
  "TRP",
  "TYR",
  "VAL"
)

# Amino Acid 3 letter to 1 letter code mapping
# Example usage: AA3to1["ALA"], will output "A"
AA3to1 = c(
  "A",
  "R",
  "N",
  "D",
  "C",
  "E",
  "Q",
  "G",
  "H",
  "I",
  "L",
  "K",
  "M",
  "F",
  "P",
  "S",
  "T",
  "W",
  "Y",
  "V"
)

# Maps 1-letter code to 3-letter code
names(AA1to3) = AA3to1

# Maps 3-letter code to 1-letter code
names(AA3to1) = AA1to3

###################################################################
# Amino Acid Code Mappings
###################################################################

# Obtain primary sequence of TM region
# @return data.frame with TM residues with columns
#   $resid - residue number
#   $resname1 - 1-letter AA name
#   $resname3 - 3-letter AA name
get_tm_region <- function(pdb_id = open_id()) {
  p_path = get_membrane_pdb_path(pdb_id)
  p = read.pdb(p_path)
  p_atoms = p$atoms[p$atoms$recname == "ATOM", ]
  p_hetatms = p$atoms[p$atoms$recname == "HETATM", ]
  memb_min_z = min(p_hetatms$x3)
  memb_max_z = max(p_hetatms$x3)
  
  select_tm_atoms = ((p_atoms$x3 >= memb_min_z) &
                       (p_atoms$x3 <= memb_max_z))
  tm_atoms = p_atoms[select_tm_atoms, ]
  
  # Obtain unique residue names (3-letter)
  select_tm_res = !duplicated(tm_atoms$resid)
  tm_res = tm_atoms[select_tm_res, c("resid", "resname")]
  colnames(tm_res)[2] = "resname3"
  resname1 = AA3to1[tm_res$resname3]
  tm_res = data.frame(resid = tm_res$resid,
                      resname1,
                      resname3 = tm_res$resname3)
  return(tm_res)
}

# Obtain primary sequence of non-tm region
# @return data.frame with with non-TM residues with columns
#   $resid - residue number
#   $resname1 - 1-letter AA name
#   $resname3 - 3-letter AA name
get_non_tm_region <- function(pdb_id = open_id()) {
  tm_res = get_tm_region(pdb_id)
  
  p_path = get_membrane_pdb_path(pdb_id)
  p = read.pdb(p_path)
  p_atoms = p$atoms[p$atoms$recname == "ATOM", ]
  
  select_res = !duplicated(p_atoms$resid)
  res = p_atoms[select_res, c("resid", "resname")]
  
  select_non_tm_res = !(res$resid %in% tm_res$resid)
  non_tm_res = res[select_non_tm_res,]
  non_tm_res = data.frame(
    resid = non_tm_res$resid,
    resname1 = AA3to1[non_tm_res$resname],
    resname3 = non_tm_res$resname
  )
  return(non_tm_res)
}

# Obtain TM region marked as bit column
# @return data.frame with with non-TM residues with columns
#   $resid - residue number
#   $resname1 - 1-letter AA name
#   $resname3 - 3-letter AA name
#   $tm - 1 if TM, 0 o/w
get_marked_tm_region <- function(pdb_id = open_id()) {
  tm_res = get_tm_region(pdb_id)
  
  p_path = get_membrane_pdb_path(pdb_id)
  p = read.pdb(p_path)
  p_atoms = p$atoms[p$atoms$recname == "ATOM", ]
  
  select_res = !duplicated(p_atoms$resid)
  res = p_atoms[select_res, c("resid", "resname")]
  
  select_tm_res = (res$resid %in% tm_res$resid)
  res = data.frame(
    resid = res$resid,
    resname1 = AA3to1[res$resname],
    resname3 = res$resname,
    tm = as.numeric(select_tm_res)
  )
  return(res)
}

# Obtain Loop region marked as numeric column
# @return data.frame with with non-TM residues with columns
#   $resid - residue number
#   $resname1 - 1-letter AA name
#   $resname3 - 3-letter AA name
#   $loop - 1, 2, 3, 5, 6, 7
get_marked_loop_regions <- function(pdb_id = open_id()) {
  p_path = get_membrane_pdb_path(pdb_id)
  p = read.pdb(p_path)
  p_atoms = p$atoms[p$atoms$recname == "ATOM", ]
  
  select_res = !duplicated(p_atoms$resid)
  res = p_atoms[select_res, c("resid", "resname")]
  
  LOOP_START = c(18, 54, 97, 177, 217, 259)
  LOOP_END = c(29, 65, 106, 188, 234, 267)
  LOOP_ID = c(1, 2, 3, 5, 6, 7)
  
  
  loop = rep(0, nrow(res))
  
  for (i in 1:length(LOOP_START)) {
    loop_start = LOOP_START[i]
    loop_end = LOOP_END[i]
    loop_id = LOOP_ID[i]
    
    select_loop = ((res$resid >= loop_start) &
                     (res$resid <= loop_end))
    loop[select_loop] = loop_id
  }
  
  res = data.frame(
    resid = res$resid,
    resname1 = AA3to1[res$resname],
    resname3 = res$resname,
    loop = loop
  )
  return(res)
}
