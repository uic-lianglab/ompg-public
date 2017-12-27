# -----------------------------------------------
# Script for generating psf structural file(s) for
# a pdb structure at parameter pH, given a propka
# pKa file.
#
# Usage:
# psfgen <path_to_this_script.tcl> <path_to_pdb> <path_to_propka_pKas> <pH>
# - where 1st argument is path to x-ray conformation pdb
# - and 2nd argument is path to output propka pKas for input pdb
# - and 3rd argument is pH of solution
# -----------------------------------------------

# -----------------------------------------------
# User arguments
# -----------------------------------------------

# Verify at least 3 arguments read
if {$argc != 5} {
    puts "Expected 3 arguments but $argc read from command line."
    puts "Usage: psfgen <in_path_to_this_script.tcl> <in_pdb_path> <in_propka_pKas_path> <pH> <out_coords_path> <out_topology_path>"
    puts "\t<in_pdb_path> - path to input protein coordinates .pdb file"
    puts "\t<in_propka_pKas_path> - path to pKas file as computed by PROPKA tool"
    puts "\t<pH> - pH of solution"
    puts "\t<out_coords_path> - path to write output coordinates .coor/.pdb file"
    puts "\t<out_topology_path> - path to write output topology .psf file"
    puts "Exiting."
    exit 0
}

# input coordinates
set pdb_path [lindex $argv 0]

# PROPKA path
set propka_path [lindex $argv 1]

# pH setting
set pH [lindex $argv 2]

# Output NAMD coordinates
set out_coords_path [lindex $argv 3]

# Output NAMD topology
set out_topo_path [lindex $argv 4]

# Output commmand line args
puts "Arguments:"
puts "\tpdb_path: $pdb_path"
puts "\tpropka_path: $propka_path"
puts "\tpH: $pH"
puts "\tout_coords_path: $out_coords_path"
puts "\tout_topo_path: $out_topo_path"

# -----------------------------------------------
# Configure forcefield
# -----------------------------------------------

# Set all paths relative to script location
set script_dir [file dirname [file normalize $argv0]]

# Base directory for NAMD tool
set base_dir [file join $script_dir ..]

# forcefield
set ff "charmm36"

# Input forcefield topology
set ff_path [file join $base_dir forcefield $ff topo.rtf]

# -----------------------------------------------
# Determine protonation state
# -----------------------------------------------

# Parse ionizable residue pKas. If
#   pH < pKa  -> protonated
#   pH >= pKa -> de-protonated

# Slurp up the propka data file
set propka_fp [open $propka_path r]
set propka_file_data [read $propka_fp]
close $propka_fp
 
# Split propka file into list of lines
set propka_lines [split $propka_file_data "\n"]
# Determine pKa summary interval
set ix_line 0
set ix_pKa_start -1
set ix_pKa_end -1
foreach line $propka_lines {
    # Skip to pKa prediction summary
    if {[string equal "SUMMARY OF THIS PREDICTION" $line]} {
        # ix_line is at section header, following line is header row
        # skip past both of these lines
        set ix_pKa_start [expr {$ix_line + 2}]
    } elseif { [expr {$ix_pKa_start >= 0}] } {
        # Determine terminal length of pKa summary
        if {[string equal -length 7 "-------" $line] } {
            set ix_pKa_end [expr {$ix_line - 1}]
            break;
        }
    }
    incr ix_line
}

# Process each pKa summary entry row
set res_name_lst {}
set res_id_lst {}
set res_pKa_lst {}
for {set ix_pKa_entry $ix_pKa_start} {$ix_pKa_entry <= $ix_pKa_end} {incr ix_pKa_entry} {
    # Tokenize by white space
    # http://stackoverflow.com/questions/13380914/how-to-split-a-string-into-a-list-of-words-in-tcl-ignoring-multiple-spaces
    # Example line:
    #    ASP   5 A    3.38      3.80
    set nfo [regexp -all -inline {\S+} [lindex $propka_lines $ix_pKa_entry]]
    # Column 0 -> residue 3-letter code
    # Column 1 -> residue pdb identifier/index
    # Column 2 -> residue chain identifier
    # Column 3 -> residue predicated pKa
    # Column 4 -> residue model (standard) pKa
    lappend res_name_lst [lindex $nfo 0]
    lappend res_id_lst [lindex $nfo 1]
    lappend res_pKa_lst [lindex $nfo 3]
}

set num_ions [llength $res_name_lst]

# -----------------------------------------------
# Read forcefield topology file
# -----------------------------------------------

topology $ff_path

# -----------------------------------------------
# Construct segment
# -----------------------------------------------

# Set default histidine protonation
# - HSD (neutral) - proton on ND nitrogen
# - HSE (neutral) - proton on NE nitrogen
# - HSP (+ charge) - proton on both nitrogens
pdbalias residue HIS HSD

# Build protein segment
segment OMPG {
    pdb $pdb_path
    # Mutate any charged histidines to HSP
    for {set ix_ion 0} {$ix_ion < $num_ions} {incr ix_ion} {
        set res_name [lindex $res_name_lst $ix_ion]
        if {[string equal "HIS" $res_name]} {
            set res_pKa [lindex $res_pKa_lst $ix_ion]
            if { [expr {$pH < $res_pKa}] } {
                set res_id [lindex $res_id_lst $ix_ion]
                mutate $res_id HSP
            }
        }
    }
}

# To (de)protonate non-histidine residues, we need to apply patches
# Available patches:
# - ASPP - protonated aspartic acid (neutral)
# - GLUP - protonated glutamic acid (neutral)
# - HS2 - switch from HSD to HSE (neutral)
# - LSN - deprotonated lysine (neutral)
for {set ix_ion 0} {$ix_ion < $num_ions} {incr ix_ion} {
    set res_name [lindex $res_name_lst $ix_ion]
    set res_pKa [lindex $res_pKa_lst $ix_ion]
    set res_id [lindex $res_id_lst $ix_ion]
    set is_prot [expr {$pH < $res_pKa}]
    if {$is_prot} {
        if {[string equal "ASP" $res_name]} {
            patch ASPP OMPG:$res_id
        } elseif {[string equal "GLU" $res_name]} {
            patch GLUP OMPG:$res_id
        }
    } else {
        if {[string equal "LYS" $res_name]} {
            patch LSN OMPG:$res_id
        }
    }  
}

# -----------------------------------------------
# Initial hydrogen placement
# -----------------------------------------------

# Aliasing for certain residues may be necessary as OMPG segment is using
# CHARMM naming conventions for residues which may differ from read in pdb
pdbalias atom ILE CD1 CD
coordpdb $pdb_path OMPG
# Guess cordinates for hydrogens
# Note: structure will need to be minimized as
# hydrogens will clash with other atoms and inflate
# the LJ potential
guesscoord

# Apparently when applying a patch, you must regenerate angles and dihedrals?
# http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/2659.html
# http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2005-2006/2904.html
# http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/4137.html
regenerate angles dihedrals 

# -----------------------------------------------
# Write NAMD topology and coordinates
# -----------------------------------------------

# Write structure and coordinate file
writepsf $out_topo_path
writepdb $out_coords_path
