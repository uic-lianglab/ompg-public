multi-loop Distance-guided chain-Growth Monte Carlo method (m-DiSGro)

m-DiSGro is an efficient protein loop sampling and structure prediction tool which simultaneously 
models of two or more interacting loops based on distance-guided sequential chain-growth Monte Carlo method.


If you use m-DiSGro, cite the following reference:

Ke Tang, et al, (2015) Conformational Sampling and Structure Prediction of Multiple Interacting Loops
in Soluble and β-Barrel Membrane Proteins Using Multi-Loop Distance-Guided
Chain Growth Monte Carlo Method. Bioinformatics. 


1. Data Preparation

   Input files must be in the PDB format. The atoms of absent loop
   regions need to be labeled as "H" atoms, whose coordinates are assigned to
   0.000 for each axis.  


2. Running m-DiSGro: 

./disgro -mode [Mode] -f [Input PDB] -n [Input # of trials] -nds [# of sampling states] -multiloops [The file recording the start and end information of loops] -confkeep [# of retained conformations] -nscc [# of side chain states] -pdbout [# of output models]


Description:

-mode: mode of the computation, currently modes: smc, torsion
       smc: loop modeling using sequential chain-groth Monte Carlo.
-f: input coordinate file of the protein
-n: number of conformations
-nds: number of sampling states
-start: The starting position for growth mode. Default is -1, starting from the beginning.
-end: The end position for growth mode. Default is -1, ending at the last residue.
-multiloops: The file recording the start and end residues of multiple loops.
-pdbout: The number of the output conformations.
-nscc: Number of side chain states

Example:

./mdisgro -mode smc -f xxx.pdb -n 10000 -nds 32 -start xx -end xx -eval -confkeep 1000 -pdbout 100

./mdisgro -mode smc -f xxx.pdb -n 10000 -nds 32 -multiloops xxx_mulist -confkeep 500 -nscc 20 -multilooplib pdb_output -pdbout 20

Contact
jliang@uic.edu
