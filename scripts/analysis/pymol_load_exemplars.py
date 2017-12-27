# Usage - open first extracted exemplar in pymol then
# copy and paste below code to load remaining

import os

# Assume working directory contains all exemplars!
pdbs = os.listdir("./")

for f in pdbs: cmd.load(f)

hide everything
show cartoon
bg_color white
color grey90
sele loop_1, resi 18-29
color purple, loop_1
sele loop_2, resi 54-65
color blue, loop_2
sele loop_3, resi 97-106
color green, loop_3
sele loop_5, resi 177-188
color forest, loop_5
sele loop_6, resi 217-234
color yellow, loop_6
sele loop_7, resi 259-267
color red, loop_7

# Assumes first loaded pdb is named '000.pdb'
# Hide all barrel atoms except for first PDB, this greatly improves render
# times, particularly when using the 'ray' command.
sele other_barrels, !(000 | loop_1 | loop_2 | loop_3 | loop_5 | loop_6 | loop_7)
hide (other_barrels)

########################################################################
# WARNING: DO NOT SET MORE THAN ONE VIEW AT A TIME BEFORE EXPORTING PNG!

# Different views: Copy and paste then save before copy and pasting next viewitems
# Captured using "get_view" entered into PyMol console

# Top
cmd.set_view("1.000000000, 0.000000000, 0.000000000, 0.000000000, 1.000000000, 0.000000000, 0.000000000, 0.000000000, 1.000000000, 0.000000000, -0.000000030, -147.997177124, -0.684683383, 2.913268328, 7.910408974, 114.226310730, 181.768066406, -20.000000000")
cmd.refresh()

# Save PNG before setting next view!

# Top with barrel (side)

#cmd.set_view("-0.431708395, -0.606993198, 0.667222500, 0.881666064, -0.127730027, 0.454258054, -0.190507159, 0.784374535, 0.590307653, 0.000000086, -0.000006557, -182.345901489, -1.191768885, -2.780379295, 12.131720543, 154.707092285, 212.884704590, -20.000000000")
### cut below here and paste into script ###
cmd.set_view ("-0.286463380, -0.453160763, 0.844147086, 0.958045423, -0.144100681, 0.247758001, 0.009367788, 0.879705310, 0.475428283, -0.000001786, 0.000012144, -178.465377808, -1.137110710, -0.282293200, 10.770358086, 136.047012329, 220.883712769, -20.000000000")
### cut above here and paste into script ###
cmd.refresh()

# Save PNG before setting next view!
 
# Bottom
cmd.set_view ("-0.930513263, 0.363615155, 0.043900818, 0.358967513, 0.929214954, -0.087765515, -0.072705910, -0.065908164, -0.995173037, 0.000008056, 0.000000402, -145.299743652, -0.754320681, 3.150283337, 23.035011292, 116.211463928, 174.389083862, -20.000000000")
cmd.refresh()

# Save PNG before setting next view!
