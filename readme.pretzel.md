**Pretzel** is based on the multi-loop Distance-guided chain-Growth Monte Carlo method (m-DiSGro)

m-DiSGro is an efficient protein loop sampling and structure prediction tool which simultaneously models two or more interacting loops based on distance-guided sequential chain-growth Monte Carlo method.

If you use Pretzel/m-DiSGro, cite the following reference:

* Perez-Rathke, Alan, Monifa Akilah Verna Fahie, Christina Chisholm, Jie Liang, and Min Chen. "Mechanism of OmpG pH-dependent gating from loop ensemble and single channel studies." Journal of the American Chemical Society (2017).
* Tang, Ke, Samuel WK Wong, Jun S. Liu, Jinfeng Zhang, and Jie Liang. "Conformational sampling and structure prediction of multiple interacting loops in soluble and β-barrel membrane proteins using multi-loop distance-guided chain-growth Monte Carlo method." Bioinformatics 31, no. 16 (2015): 2646-2652.
* Tang, Ke, Jinfeng Zhang, and Jie Liang. "Fast protein loop sampling and structure prediction using distance-guided sequential chain-growth Monte Carlo method." PLoS computational biology 10, no. 4 (2014): e1003539.

### Data Preparation

Input files must be in the PDB format.

### Running Pretzel: 

```
./pretzel -mode [Mode] -f [Input PDB] -n [Input # of trials] -nds [# of sampling states] -multiloops [The file recording the start and end information of loops] -confkeep [# of retained conformations] -nscc [# of side chain states] -pdbout [# of output models]
```

### Description:

* -mode: mode of the computation, currently modes: smc, torsion
       smc: loop modeling using sequential chain-groth Monte Carlo.
* -f: input coordinate file of the protein
* -n: number of conformations
* -nds: number of sampling states
* -start: The starting position for growth mode. Default is -1, starting from the beginning.
* -end: The end position for growth mode. Default is -1, ending at the last residue.
* -multiloops: The file recording the start and end residues of multiple loops.
* -pdbout: The number of the output conformations.
* -nscc: Number of side chain states

### Example:

```
./pretzel -mode smc -f xxx.pdb -n 10000 -nds 32 -start xx -end xx -eval -confkeep 1000 -pdbout 100
```

```
./pretzel -mode smc -f xxx.pdb -n 10000 -nds 32 -multiloops xxx_mulist -confkeep 500 -nscc 20 -multilooplib pdb_output -pdbout 20
```

### Contact

perezratATuicDOTedu
jliangATuicDOTedu
