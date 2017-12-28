Pretzel uses Distance-guided Sequential chain-Growth Monte Carlo(DiSGro)

Files required for running Pretzel. 

atomProp2.txt:                atom and residue information e.g. atom type.
BBT_phi_psi_pair_NEW.txt:     joint probability of phi and psi torsion angles.
frag.N.C_pdf_32_19.txt:       empirical distance-guidence distribution of N and C atoms.
frag.C.CA_pdf_32_19.txt:      empirical distance-guidence distribution of C and CA atoms. 
LOODIS_ed4_8_V3.txt:          distance-dependent empirical potential function.
LoopGeo_37_pdf_21.txt:        useful information for loop geometry.
parameter.txt:                coefficients of energy terms.  
SCT_PF.txt:                   side-chain torsions distribution. For calculating smoothed counts, 
                              the mean and standard deviations of Gaussian kernels are set to 
                              observed torsion values and the standard deviations of torsions 
                              derived from the PDB database.

Contact

perezratATuicDOTedu
jliangATuicDOTedu
