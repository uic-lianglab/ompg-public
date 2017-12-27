###################################################################
# load_data_utils.R
#
# Utilities for loading, parsing, and formatting specific CSV data
# sets from disk. CSV datasets loaded in this manner will have
# cached binary rdata archives written to disk to improve future
# load times.
###################################################################

###################################################################
# Data path utilities
###################################################################

# @return file name without path information or last .csv extension
get_fid_from_csv <- function(csv_path) {
  return(gsub("\\.csv$", "", basename(csv_path)))
}

# @return file name without path information or last .rdata extension
get_fid_from_rdata <- function(rdata_path) {
  return(gsub("\\.rdata$", "", basename(rdata_path)))
}

# @return modified file path of form:
#   <directory>/<basename of file>.<attrib>.rdata
insert_attrib_in_rdata <- function(rdata_path, attrib) {
  d = dirname(rdata_path)
  fid = get_fid_from_rdata(rdata_path)
  return(file.path(d, paste0(fid, ".", attrib, ".rdata")))
}

# Determines binary RDATA path from input CSV path
# @return rdata path
get_rdata_path_from_csv <- function(csv_path) {
  return(gsub("\\.csv$", ".rdata", csv_path))
}

# Determines ASCII CSV path from RDATA path
# @return csv path
get_csv_path_from_rdata <- function(rdata_path) {
  return(gsub("\\.rdata$", ".csv", rdata_path))
}

# Determines derived from rdata path for file with new extension
get_ext_path_from_rdata <- function(rdata_path, new_ext_with_dot) {
  return(gsub("\\.rdata$", new_ext_with_dot, rdata_path))
}

# Determines binary RDATA path from input CSV path but cropped
# based on energy ranks
# @return rdata path
get_cropped_rdata_path_from_csv <- function(csv_path, max_len) {
  return(gsub("\\.csv$",
              paste0(".", max_len, ".rdata"),
              csv_path))
}

# @return path to superset data, must be checked for existence
get_superset_path_from_csv <-
  function(csv_path, target_len, ext = ".rdata") {
    d = dirname(csv_path)
    # Return names of <ext> files in current directory
    l = list.files(path = d,
                   pattern = paste0(".*\\.", strsplit(
                     x = ext, split = ".", fixed = TRUE
                   )[[1]][2]))
    
    fid = get_fid_from_csv(csv_path)
    fid_w_dot = paste0(fid, ".")
    
    nchar_fid_w_dot = nchar(fid_w_dot)
    nchar_ext = nchar(ext)
    
    min_count = .Machine$integer.max
    out_fpath = ""
    
    for (f in l) {
      # Length accounts for . and at least an integer
      # Find file names of form <fid>.<length>.rdata
      # Skip files if too small
      if (nchar(f) < (nchar_fid_w_dot + nchar_ext + 1)) {
        next
      }
      
      # Skip files that don't start with fid
      if (substr(f, 1, nchar_fid_w_dot) != fid_w_dot) {
        next
      }
      
      # Extract rank count:
      # Strip .rdata extension
      f2 = strsplit(x = f,
                    split = ext,
                    fixed = TRUE)[[1]][1]
      # Strip fid
      f3 = strsplit(x = f2,
                    split = fid_w_dot,
                    fixed = TRUE)[[1]][2]
      
      # Skip file names which were not able to be parsed
      count = as.integer(f3)
      fpath = file.path(d, f)
      if (is.na(count)) {
        print(paste0("Warning, unable to parse for possible superset: ", fpath))
        next
      }
      
      if ((count >= target_len) && (count < min_count)) {
        min_count = count
        out_fpath = fpath
      }
    }
    
    return(out_fpath)
  }

###################################################################
# Convert csv to rdata utilities
###################################################################

# Load an existing rdata file
# @rdata_path - path on disk to rdata file, data is assumed to
# have been saved by 'save_rdata'
# @return loaded data
load_rdata <- function(rdata_path) {
  # Loads into variable named 'data'
  print(paste0("Loading ", rdata_path))
  load(file = rdata_path)
  return(data)
}

# Save to an rdata file - data saved by this method can be
# loaded using load_rdata(<path>)
# @param data - the data to save to rdata archive on disk
# @param rdata_path - path to write rdata archive
# @param b_create_subdirs - TRUE if subdirectories should be created
save_rdata <- function(data, rdata_path, b_create_subdirs = TRUE) {
  print(paste0("Saving ", rdata_path))
  d = dirname(rdata_path)
  if (b_create_subdirs && !dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  save(data, file = rdata_path)
}

# Ensures consistent line endings for plain-text
# @param open should be set to binary mode:
#   "wb" -> open for writing in binary mode
#   "ab" -> open for appending in binary mode
save_csv <- function(x,
                     csv_path,
                     row.names = FALSE,
                     col.names = FALSE,
                     open = "wb",
                     create_subdirs = TRUE) {
  print(paste0("Saving ", csv_path))
  d = dirname(csv_path)
  if (create_subdirs && !dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  # http://r.789695.n4.nabble.com/Write-table-eol-argument-td3159733.html
  # Ensure consistent line endings (no carriage return) by opening file in
  # binary mode
  f = file(csv_path, open = open)
  write.table(
    x = x,
    file = f,
    sep = ",",
    eol = "\n",
    row.names = row.names,
    col.names = col.names
  )
  close(f)
}

# Load CSV data set. Will check for existing cached rdata first
# and load that if found. Else, will load csv and create rdata.
# @param csv_path - path to csv data file to be loaded
# @param header - TRUE if data has header row
# @param to_matrix - TRUE if data is homogeneous and should be
#   converted to matrix
# @param overwrite - Force load from csv path and overwrite cache
# @param ... - additional arguments to feed to read.csv()
load_csv_cached <- function(csv_path,
                            header = TRUE,
                            to_matrix = FALSE,
                            overwrite = FALSE,
                            ...) {
  rdata_path = get_rdata_path_from_csv(csv_path = csv_path)
  # Skip file if rdata exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0("Skipping csv load as rdata already exists: ",
                 rdata_path))
    return(load_rdata(rdata_path))
  }
  
  print(paste0("Loading ", csv_path))
  res = read.csv(file = csv_path, header = header, ...)
  
  if (to_matrix) {
    res = as.matrix(res)
  }
  
  save_rdata(data = res, rdata_path = rdata_path)
  return(res)
}

# Load PDB mappings from disk
# @param pdb_path - path to a PDB file
# @return data.frame with columns:
#   $res_name_map: mapping from residue sequence number to 3 letter AA code
#   $loop_map: mapping from residue sequence number to LOOP_# string
load_pdb_info <- function(pdb_path) {
  library(Rpdb)
  
  prot = read.pdb(pdb_path)
  prot = prot$atoms
  # Strip heteroatoms
  prot = prot[prot$recname == "ATOM",]
  
  # Mapping from residue sequence number to residue name
  res_name_map = rep("NONE", length(unique(prot$resid)))
  res_name_map[unique(prot$resid)] = prot$resname[!duplicated(prot$resid)]
  res_name_map = as.character(res_name_map)
  
  # Mapping from residue sequence number to loop region
  loop_map = rep("NONE", length(res_name_map))
  loop_map[18:29] = "LOOP_1"
  loop_map[54:65] = "LOOP_2"
  loop_map[97:106] = "LOOP_3"
  loop_map[177:188] = "LOOP_5"
  loop_map[217:234] = "LOOP_6"
  loop_map[259:267] = "LOOP_7"
  
  return(data.frame(res_name_map, loop_map, stringsAsFactors = FALSE))
}

# Utility for loading sorted energy values from disk. If no k-value specified,
# will crop outliers for least stable extreme (+1.5 * 1QR).
#
# Writes binary version of CSV data for faster load times in future.
#
# @param csv_path - path to raw ASCII data
# @param k - energy rank cutoff, only top k energy profiles will be kept
#   if <= 0, then no cutoff is applied
# @param overwrite - force overwrite of binary cached rdata
# @return data.frame with columns:
#   $NAME - PDB identifier - unique name of protein sample
#   $PH - pH at which protonation states were set
#   $POTENTIAL - total potential energy (summation of computed energies)
#   $ELECT - electrostatic energy
#   $VDW - Van der Waals or LJ energy
#   $BOND - bond length energy
#   $ANGLE - bond angle energy
#   $DIHED - dihedral angle energy
#   $IMPRP - improper energy
#   $BOUNDARY - boundary energy
#   $MISC - miscellaneous energy
ener_csv_2_rdata <- function(csv_path,
                             k = -1,
                             overwrite = FALSE) {
  # Determine path to binary cache
  rdata_path = c()
  if (k <= 0) {
    # No cutoff specified, get path to full data set
    rdata_path = get_rdata_path_from_csv(csv_path)
  }
  else {
    # Get path to cropped data set
    rdata_path = get_cropped_rdata_path_from_csv(csv_path, k)
  }
  
  # Skip file if rdata exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping csv conversion as rdata already exists: ",
      rdata_path
    ))
    return(load_rdata(rdata_path))
  }
  
  # Check if superset exists
  if (!overwrite && (k > 0)) {
    rdata_superset_path = get_superset_path_from_csv(csv_path, k)
    if (file.exists(rdata_superset_path)) {
      superset = load_rdata(rdata_superset_path)
      data = superset[1:k, ]
      # Write binary, compressed data
      save_rdata(data, rdata_path)
      return(data)
    }
  }
  
  # Load ASCII .csv
  print(paste0("Reading: ", csv_path))
  data = read.csv(file = csv_path,
                  header = TRUE,
                  stringsAsFactors = FALSE)
  
  if (k > 0) {
    # Crop energy profiles
    data = data[1:k, ]
  }
  else {
    # Crop outliers:
    
    # IQR Rule for Outliers
    # 1. Arrange data in order.
    # 2. Calculate first quartile (Q1), third quartile (Q3) and the interquartile
    # range (IQR = Q3 - Q1).
    # 3. Compute Q1 - 1.5 x IQR, Compute Q3 + 1.5 x IQR. Anything outside this range is an outlier.
    
    # Trim  Q3 + 1.5 x IQR as outliers
    s = summary(data$POTENTIAL)
    q1 = s["1st Qu."]
    q3 = s["3rd Qu."]
    iqr = q3 - q1
    max_ener = q3 + (1.5 * iqr)
    filt = data$POTENTIAL <= max_ener
    n_trim = length(filt) - sum(filt)
    perc_trim = 100.0 * n_trim / length(filt)
    print(paste0("Trimming ", n_trim, " outliers (", perc_trim, "%)"))
    print(paste0("Outliers have energy > ", max_ener))
    data = data[filt, ]
    
    # Write superset filename if no k value specified
    max_k_len = nrow(data)
    rdata_max_k_path = get_cropped_rdata_path_from_csv(csv_path, max_k_len)
    # Write binary, compressed data
    save_rdata(data, rdata_max_k_path)
  }
  
  # Write binary, compressed data
  save_rdata(data, rdata_path)
  return(data)
}

# Utility for loading open/closed status from disk.
#
# Writes binary version of CSV data for faster load times in future.
#
# @param csv_path - path to CSV file containing open (1) vs close (0) status
# @param ener_ids - identifiers used for sorting open-close status by energy rank
# @param overwrite - if TRUE, will overwrite any stored rdata on disk
# @return  logical vector with TRUE = open, FALSE = close. The names
#   attribute of the returned logical is the same as ener_ids parameter
oc_csv_2_rdata <- function(csv_path, ener_ids, overwrite = FALSE) {
  target_len = length(ener_ids)
  rdata_path = get_cropped_rdata_path_from_csv(csv_path, target_len)
  # Skip file if rdata exists
  if (!overwrite && file.exists(rdata_path)) {
    print(paste0(
      "Skipping csv conversion as rdata already exists: ",
      rdata_path
    ))
    return(load_rdata(rdata_path))
  }
  
  data = c()
  oc = c()
  
  # Check if a superset file exists
  rdata_superset_path = get_superset_path_from_csv(csv_path, target_len)
  if (!overwrite && file.exists(rdata_superset_path)) {
    print(paste0("Loading superset rdata: ", rdata_superset_path))
    oc = load_rdata(rdata_superset_path)
  }
  else {
    # Load ASCII .csv
    print(paste0("Reading: ", csv_path))
    data = read.csv(file = csv_path,
                    header = FALSE,
                    stringsAsFactors = FALSE)
    
    # Convert to logical values = TRUE if open, FALSE if closed
    oc = as.logical(data[, 2])
    # Sort by energy ranking
    names(oc) = data[, 1]
  }
  
  # Filter data set
  ener_filt = ener_ids %in% names(oc)
  ener_ids = ener_ids[ener_filt]
  data = oc[ener_ids]
  
  # Write binary, compressed data
  save_rdata(data, rdata_path)
  return(data)
}

# Utility for loading electrostatic interaction data from disk.
#
# Writes binary version of CSV data for faster load times in future.
#
# @param csv_path - path to CSV file containing raw electrostatic interactions
#   CSV file should have columns: RESA, RESB, <sample_id1>, ... <sample_idN>
#   WARNING: csv_path needs to be of form <dir>/<fid>.csv
#   The file name cannot have a length - e.g> <fid>.<length>.csv else CSV
#   supersets will not be found
# @param ener_ids - identifiers used for sorting and cropping electrostatic data
# @param oc - logical vector with open-close status for each energy identifier
#   in ener_ids with TRUE -> open, FALSE -> close
# @param min_obs_count - Interaction pairs with fewer than this number of
#   observations are removed from final data set
# @param overwrite - if TRUE, will overwrite any stored rdata caches on disk
# @return list with members:
#   $design - a data.frame with columns:
#     <interaction_id1>, ..., <interaction_idN> where each interaction
#       identifier is created by merging the residue sequence numbers of the
#       interaction pair into a single character string.
#     The rownames are the same as ener_ids.
#     Therefore, each column is the electrostatic interaction energy for that
#     interaction pair at each of the protein samples
#   $meta - a data.frame with meta-information about each interaction pair.
#     Columns include:
#     $RESA.NO - residue sequence number for first residue in interaction pair
#     $RESB.NO - residue sequence number for second residue in interaction pair
#     $RESA.NAME - three letter AA code for first residue in interaction pair
#     $RESB.NAME - three letter AA code for second residue in interaction pair
#     $RESA.LOOP - name of loop region or 'NONE' for first residue in interaction pair
#     $RESB.LOOP - name of loop region or 'NONE' for second residue in interaction pair
#     $N - total number of observations for each interaction pair
#     $N.OPEN - total number of observations from open conformations for each interaction pair
#     $N.CLOSe - total number of observations from closed conformations for each interaction pair
#     $ION.ION - If TRUE, both residues in interaction pair are ionic residues, FALSE o/w
#     $POL.POL - If TRUE, both residues in interaction pair are polar residues, FALSE o/w
#     $LOOP.LOOP - IF TRUE, both residues in interaction pair residue in named loop regions, FALSE o/w
#     $LOOP.ION - logical equivalent to LOOP.LOOP & ION.ION
#     $LOOP.POL - logical equivalent to LOOP.LOOP & POL.POL
#     $RESA.ION - If TRUE, then first residue in interaction pair is ionic residue, FALSE o/w
#     $RESB.ION - If TRUE, then second residue in interaction pair is ionic residue, FALSE o/w
#     $RESA.POL - If TRUE, then first residue in interaction pair is polar residue, FALSE o/w
#     $RESB.POL - If TRUE, then second residue in interaction pair is polar residue, FALSE o/w
#     $MED.OPEN - Median among observed open samples at each interaction pair, NA if no observed samples
#     $MED.CLOSE - Median among observed closed samples at each interaction pair, NA if no observed samples
#     $MIN.OPEN - Minimum value among observed open samples at each interaction pair, NA if no observed samples
#     $Q1.OPEN - 25% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MEA.OPEN - Mean value among observed open samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.OPEN - 75% quantile value among observed open samples at each interaction pair, NA if no observed samples
#     $MAX.OPEN - Maximum value among observed open samples at each interaction pair, NA if no observed samples
#     $MIN.CLOSE - Minimum value among observed closed samples at each interaction pair, NA if no observed samples
#     $Q1.CLOSE - 25% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MEA.CLOSE - Mean value among observed closed samples at each interaction pair, *NaN* (not NA) if no observed samples
#     $Q3.CLOSE - 75% quantile value among observed closed samples at each interaction pair, NA if no observed samples
#     $MAX.CLOSE - Maximum value among observed closed samples at each interaction pair, NA if no observed samples
#     $FAVORS - At each interaction pair, values are: (note: NA -> not observed)
#       "OPEN" if (MED.OPEN < MED.CLOSE) or (CLOSE is NA and MED.OPEN < 0) or (OPEN is NA and MED.CLOSE > 0)
#       "CLOSE" if (MED.CLOSE < MED.OPEN) or (OPEN is NA and MED.CLOSE < 0) or (CLOSE is NA and MED.OPEN > 0)
#       "NONE" o/w
elec_csv_2_rdata <-
  function(csv_path,
           template_pdb_info,
           ener_ids,
           oc,
           min_obs_count = 100,
           overwrite = FALSE) {
    if (min_obs_count <= 0) {
      print("Warning, min observation count <= 0, setting to 1 instead")
      min_obs_count = 1
    }
    
    # Check if data already exists
    dr = dirname(csv_path)
    fid = get_fid_from_csv(csv_path)
    # Add path modifiers here
    # HACKY - this csv path should never actually exist - just needed for
    # transforming to rdata path
    mod_csv_path = file.path(dr,
                             paste0(fid,
                                    ".moc", min_obs_count,
                                    ".csv"))
    
    target_len = length(ener_ids)
    rdata_path = get_cropped_rdata_path_from_csv(mod_csv_path, target_len)
    # Skip file if rdata exists
    if (!overwrite && file.exists(rdata_path)) {
      print(paste0(
        "Skipping csv conversion as rdata already exists: ",
        rdata_path
      ))
      return(load_rdata(rdata_path))
    }
    
    # Check if binary version of pre-transformed ASCII data exists
    # HACKY - this csv path should never actually exist - just needed for
    # transforming to rdata path
    raw_csv_path = file.path(dr, paste0(fid, ".raw.csv"))
    raw_rdata_path = get_cropped_rdata_path_from_csv(raw_csv_path, target_len)
    raw_data = data.frame()
    if (!overwrite && file.exists(raw_rdata_path)) {
      # Binary exists
      raw_data = load_rdata(raw_rdata_path)
    }
    else {
      # Check if superset data exists
      raw_rdata_superset_path = get_superset_path_from_csv(raw_csv_path, target_len)
      if (!overwrite && file.exists(raw_rdata_superset_path)) {
        raw_data = load_rdata(raw_rdata_superset_path)
      }
      else {
        # Load ASCII .csv
        # Check if superset CSV file exists
        final_csv_path = csv_path
        raw_csv_superset_path = get_superset_path_from_csv(csv_path, target_len, ext =
                                                             ".csv")
        if (file.exists(raw_csv_superset_path)) {
          final_csv_path = raw_csv_superset_path
        }
        
        print(paste0("Reading: ", final_csv_path))
        raw_data = read.csv(file = final_csv_path,
                            header = TRUE,
                            stringsAsFactors = FALSE)
        unfilt_len = ncol(raw_data) - 2
        unfilt_rdata_path = get_cropped_rdata_path_from_csv(raw_csv_path, unfilt_len)
        save_rdata(raw_data, unfilt_rdata_path)
      }
      
      # Keep only the raw columns with matching energy ids
      ener_id_filt = ener_ids %in% colnames(raw_data)
      if (!all(ener_id_filt)) {
        print("Warning, energy identifiers missing in electrostatic interactions data set.")
      }
      ener_ids = ener_ids[ener_id_filt]
      raw_data = data.frame(RESA = raw_data$RESA, RESB = raw_data$RESB, raw_data[, ener_ids])
      
      # Save to disk so we don't have to parse ASCII in future
      save_rdata(raw_data, raw_rdata_path)
    }
    
    # The input columns are: ResA.No, ResB.No, SampleId1, SampleId2, ...
    # For working with random forests, etc, we need to transform
    # the data to be a tuple where 1st element is design matrix for
    # random forest with columns (Id1:ResA1,ResB1), (Id2:ResA1,ResB1), ...
    # with rownames being the sample identifiers (so transpose of input CSV)
    # and 2nd element is a list of meta information such as:
    # Ids, ResA.No, ResB.No, ResA.Name, ResB.Name, ResA.Loop, ResB.Loop, N, N.Open, N.Close, Favors, Favors.Strong, Med.Open, Med.Close
    # also: Min.Close, Q1.Close, Q3.Close, Mean.Close, Max.Close, ...
    # also: IsLoopLoop, IsIonIon, IsLILI
    
    oc = oc[ener_ids]
    
    # Create design frame (must transpose and create a unique key for each res pair)
    res_keys = paste0(raw_data$RESA, ".", raw_data$RESB)
    rownames(raw_data) = res_keys
    design = t(raw_data[, ener_ids])
    
    # Logical data.frame for observed interactions
    observed = !is.na(design)
    # Number of observations at each interaction pair
    N = colSums(observed)
    
    # Filter interactions with too few total observations
    keep = (N >= min_obs_count)
    res_keys = res_keys[keep]
    raw_data = raw_data[res_keys,]
    design = design[, res_keys]
    observed = observed[, res_keys]
    N = N[keep]
    
    observed_open = observed[oc,]
    # Number of observations at each interaction pair that are open
    N.OPEN = colSums(observed_open)
    # Number of observations at each interaction pair that are closed
    N.CLOSE = N - N.OPEN
    
    # Create meta info frame
    # Sequence number for first residue in pair
    RESA.NO = raw_data$RESA
    # Sequence number for second residue in pair
    RESB.NO = raw_data$RESB
    # 3-letter AA code for first residue in pair
    RESA.NAME = template_pdb_info$res_name_map[raw_data$RESA]
    # 3-letter AA code for second residue in pair
    RESB.NAME = template_pdb_info$res_name_map[raw_data$RESB]
    # Name of loop region (or NONE) for first residue in pair
    RESA.LOOP = template_pdb_info$loop_map[raw_data$RESA]
    # Name of loop region (or NONE) for second residue in pair
    RESB.LOOP = template_pdb_info$loop_map[raw_data$RESB]
    
    # Names of ionic (charged) residues
    ion_names = c("ARG", "LYS", "ASP", "GLU", "HIS")
    # Names of polar (partial charged) residues
    polar_names = c("GLN", "ASN", "SER", "THR", "TYR", "CYS", "MET", "TRP")
    
    # TRUE if first residue in pair is ionic
    RESA.ION = RESA.NAME %in% ion_names
    # TRUE if second residue in pair is ionic
    RESB.ION = RESB.NAME %in% ion_names
    # TRUE if first residue in pair is polar
    RESA.POL = RESA.NAME %in% polar_names
    # TRUE if second residue in pair is polar
    RESB.POL = RESB.NAME %in% polar_names
    
    # TRUE if both residues in pair are ionic (ion-ion interaction)
    ION.ION = RESA.ION & RESB.ION
    # TRUE if both residues in pair are polar
    POL.POL = RESA.POL & RESB.POL
    # TRUE if both residues in pair are found in loop regions
    LOOP.LOOP = (RESA.LOOP != "NONE") & (RESB.LOOP != "NONE")
    # TRUE if both residues in pair are within loop regions and ionic
    LOOP.ION = LOOP.LOOP & ION.ION
    # TRUE if both residues in pair are within loop regions and polar
    LOOP.POL = LOOP.LOOP & POL.POL
    
    # Compute summaries
    
    get_summary_feat <- function(s, ix) {
      out = strsplit(x = s[ix, ],
                     split = ":",
                     fixed = TRUE)
      out = as.data.frame(out, stringsAsFactors = FALSE)
      # Suppress "Warning: NAs introduced by coercion"
      out = suppressWarnings(as.numeric(out[2, ]))
      return(out)
    }
    
    get_summary <- function(d) {
      IX_MIN = 1
      IX_Q1 = 2
      IX_MED = 3
      IX_MEA = 4
      IX_Q3 = 5
      IX_MAX = 6
      s = summary(d)
      MIN = get_summary_feat(s, IX_MIN)
      Q1 = get_summary_feat(s, IX_Q1)
      MED = get_summary_feat(s, IX_MED)
      MEA = get_summary_feat(s, IX_MEA)
      Q3 = get_summary_feat(s, IX_Q3)
      MAX = get_summary_feat(s, IX_MAX)
      return(data.frame(MIN, Q1, MED, MEA, Q3, MAX, stringsAsFactors = FALSE))
    }
    
    s_df.open = get_summary(design[oc,])
    s_df.close = get_summary(design[!oc,])
    
    # Determine favorability
    # Medians
    open_med = s_df.open$MED
    close_med = s_df.close$MED
    # Median exists
    open_exists = !is.na(open_med)
    close_exists = !is.na(close_med)
    all_exists = open_exists & close_exists
    # Exclusive median existance
    open_only = open_exists & (!close_exists)
    close_only = close_exists & (!open_exists)
    
    # Relative favorability
    favors_open = rep(FALSE, length(open_med))
    # If interaction is observed for only one class, then favor if negative potential
    favors_open[open_only] = open_med[open_only] < 0
    # Both medians exists, check to see which is relatively more stable
    favors_open[all_exists] = open_med[all_exists] < close_med[all_exists]
    
    # Perform similar operation at closed samples
    favors_close = rep(FALSE, length(close_med))
    favors_close[close_only] = close_med[close_only] < 0
    favors_close[all_exists] = close_med[all_exists] < open_med[all_exists]
    
    # Tag positive interactions that are missing from other class as
    # favoring the other class
    favors_open[close_only] = close_med[close_only] > 0
    favors_close[open_only] = open_med[open_only] > 0
    
    # These sets should be disjoint
    stopifnot(!any(favors_open & favors_close))
    
    # Combine logical vectors to form single FAVOR column
    FAVORS = rep("NONE", length(open_med))
    FAVORS[favors_open] = "OPEN"
    FAVORS[favors_close] = "CLOSE"
    
    # Create meta data.frame
    meta = data.frame(
      RESA.NO,
      RESB.NO,
      RESA.NAME,
      RESB.NAME,
      RESA.LOOP,
      RESB.LOOP,
      N,
      N.OPEN,
      N.CLOSE,
      ION.ION,
      POL.POL,
      LOOP.LOOP,
      LOOP.ION,
      LOOP.POL,
      RESA.ION,
      RESB.ION,
      RESA.POL,
      RESB.POL,
      MED.OPEN = s_df.open$MED,
      MED.CLOSE = s_df.close$MED,
      MIN.OPEN = s_df.open$MIN,
      Q1.OPEN = s_df.open$Q1,
      MEA.OPEN = s_df.open$MEA,
      Q3.OPEN = s_df.open$Q3,
      MAX.OPEN = s_df.open$MAX,
      MIN.CLOSE = s_df.close$MIN,
      Q1.CLOSE = s_df.close$Q1,
      MEA.CLOSE = s_df.close$MEA,
      Q3.CLOSE = s_df.close$Q3,
      MAX.CLOSE = s_df.close$MAX,
      FAVORS,
      stringsAsFactors = FALSE
    )
    
    rownames(meta) = res_keys
    
    # Combine into list structure
    data = list()
    data[["design"]] = design
    data[["meta"]] = meta
    save_rdata(data, rdata_path)
    return(data)
  }
