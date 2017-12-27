#!/bin/bash

# Script calls all scripts for multiloop clustering

# Path to script directory
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

pushd ./

cd $SCRIPT_DIR

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling PDB Dirs 2 Cluster CSV"

./a_pdb_dirs_2_cluster_csv.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling Cluster CSV 2 Cluster Bin"

./b_cluster_csv_2_cluster_bin.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling PCA Reduce"

./c_pca_reduce.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling SHWAP Cluster"

./d_cluster_shwap.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling Resample IDs"

./e_resample_ids.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling Resample PDBs"

./f_resample_pdbs.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling Generate NAMD Topology"

./h_gen_NAMD_topo.sh

echo ""
echo "---------------------------------------------------"
echo "---------------------------------------------------"
echo "Calling PSF Clean"

./i_psf_clean.sh

popd

echo "Finished pipeline."

