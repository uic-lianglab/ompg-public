#!/bin/bash

# Script calls all scripts for relaxed C beta data generation

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

pushd ./

cd $SCRIPT_DIR

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling NAMD minimize"

./a_NAMD_minimize.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling NAMD capture"

./b_NAMD_capture.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Cat NAMD energy"

./c_1_cat_NAMD_energy.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Munge NAMD interactions"

./c_2_munge_NAMD_interactions.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling coor to CASTp PDB converter"

./d_1_coor_to_castp_pdb.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling CASTp on converted PDBs"

./e_call_castp.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Classifiying PDBs as open or close state"

./f_open_closed_classify.sh

popd

echo "Finished pipeline."
