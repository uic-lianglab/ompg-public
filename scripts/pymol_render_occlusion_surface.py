# Usage:
# pymol -qrc <path_to_this_script> -- <path_to_in_pdb> <path_to_out_png>
# Script wil take an input PDB, render it as a surface, cancel effects
# of shading by cranking up ambient lighting, orient camera looking towards
# +z axis, and then export PNG file. Note: on windows, pymol should be in
# your PATH variable and likely will need to call pymolwin rather than pymol.
# http://www.pymolwiki.org/index.php/Command_Line_Options
# -q Quiet launch. Suppress splash screen & other chatter.
# -r Run a Python program (in __main__) on startup
# -c Command line mode, no GUI. For batch operations.
# http://www.pymolwiki.org/index.php/Scripting_FAQs
# The '--' is to tell pymol to stop parsing the command line
# and feed the remaining arguments to the python script.
from pymol import cmd

import sys

# Expected command line arguments
ARG_IN_PDB_PATH = 1
ARG_OUT_PNG_PATH = 2

# Obtain command line PDB path
pdb_path = sys.argv[ARG_IN_PDB_PATH]
# Obtain command line PNG path
png_path = sys.argv[ARG_OUT_PNG_PATH]

print "Loading " + pdb_path
# Load PDB
cmd.load(pdb_path)
# Select loop regions
cmd.do("hide everything")
cmd.do("bg_color white")
cmd.do("color cyan")
cmd.do("sele loop_1, resi 18-29")
cmd.do("color purple, loop_1")
cmd.do("sele loop_2, resi 54-65")
cmd.do("color blue, loop_2")
cmd.do("sele loop_3, resi 97-106")
cmd.do("color green, loop_3")
cmd.do("sele loop_5, resi 177-188")
cmd.do("color forest, loop_5")
cmd.do("sele loop_6, resi 217-234")
cmd.do("color yellow, loop_6")
cmd.do("sele loop_7, resi 259-267")
cmd.do("color red, loop_7")
cmd.deselect()
# Show as Connolly surface (probe radius = 1.4 A = H20)
# Note: probe radius can be changed by:
# cmd.set("solvent_radius", <new radius in Angstrom>)
# http://www.pymolwiki.org/index.php/Surface
cmd.show("surface")
# Crank up the ambient lighting to get rid of shading effects
# http://pymol.org/dokuwiki/doku.php?id=setting:light
cmd.set("ambient", 3.0)
# Turn off specular highlights
cmd.set("specular", "off")
# Set direct light intensity to 0
cmd.set("direct", 0.0)
# Set number of camera lights to 0
cmd.set("light_count", 1)
# Turn off depth cue "fog" effect
cmd.set("depth_cue", 0)
# Set to orthographic rendering mode
cmd.set("orthoscopic", "on")
# Assume beta-barrel is oriented along z-axis
# with top of barrel at +z. We want to orient the
# camera so that is below the barrel looking toward
# +z axis (looking from bottom to top). This was
# done using the following:
#   cmd.turn("y", 180.0) # rotates camera
#   get_view # dumps view matrix
# The matrix can then be fed to set_view in order
# to recreate the scene
# http://www.pymolwiki.org/index.php/Get_View
   
cmd.set_view (\
"""
    -1.000000000,    0.000000000,    0.000000087,\
     0.000000000,    1.000000000,    0.000000000,\
    -0.000000087,    0.000000000,   -1.000000000,\
     0.000000000,    0.000000000, -160.547195435,\
     0.440752983,    0.527666092,    7.564197540,\
   126.576637268,  194.517761230,   20.000000000
"""
)

cmd.refresh()
# Export to PNG
print "Exporting " + png_path
cmd.png(png_path)
