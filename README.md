# Loop simulations for Outer Membrane Protein G

If using this source code, please cite the following:

* Perez-Rathke, Alan, Monifa Akilah Verna Fahie, Christina Chisholm, Jie Liang, and Min Chen. "Mechanism of OmpG pH-dependent gating from loop ensemble and single channel studies." Journal of the American Chemical Society (2017).
* Tang, Ke, Samuel WK Wong, Jun S. Liu, Jinfeng Zhang, and Jie Liang. "Conformational sampling and structure prediction of multiple interacting loops in soluble and β-barrel membrane proteins using multi-loop distance-guided chain-growth Monte Carlo method." Bioinformatics 31, no. 16 (2015): 2646-2652.
* Tang, Ke, Jinfeng Zhang, and Jie Liang. "Fast protein loop sampling and structure prediction using distance-guided sequential chain-growth Monte Carlo method." PLoS computational biology 10, no. 4 (2014): e1003539.

This code is available from the following git repositories:

* [github](https://github.com/uic-lianglab/ompg-public)
* [bitbucket](https://bitbucket.org/aperezrathke/ompg-public)

## Compiling

### Unix

In Unix, navigate to the *scripts* folder and run

```
./cmake_init.sh
```

followed by

```
./cmake_build.sh.
```

The pretzel executable binary will be located in the CMakeBuild folder.

### Windows

In Windows, navigate to the *vstudio* folder and open the solution (.sln) in Visual Studio 2017 (or later) and then press *ctrl+shift+b* to compile.


## Usage

Proteins are modeled in 4 sequential stages

1. Fragment library generation
2. Multi-loop generation
3. Multi-loop clustering
4. Relaxation and topological gating state assignment

All scripts assume the compiled pretzel executable binary has been copied to the base project directory!

### Fragment library generation

Prior to multi-loop modeling, all single loops must have their fragment libraries generated using pretzel.

For scripts demonstrating fragment library generation as used in JACS paper, please see [scripts/a_frag_lib](scripts/a_frag_lib).

All scripts are prefixed in the order in which they should be run.

### Multi-loop generation

Once fragment libraries are available, all loops structures may be simultaneously modeled using pretzel.

For scripts demonstrating full multi-loop modeling, please see [scripts/b_multi_loop](scripts/b_multi_loop).

### Multi-loop clustering

To approximate uniform coverage over the geometric space of loops, pretzel sampling bias must be removed.

This is currently done by clustering the multi-loop samples using a streaming affinity propagation algorithm.

The source code for the clustering utility (*SCOBY*) is available [here](https://bitbucket.org/aperezrathke/scoby).

For scripts demonstrating clustering, please see [scripts/c_cluster](scripts/c_cluster).

All scripts assume the *SCOBY* binary is located in the base project folder.

### Relaxation and topological gating state assignment

To simulate a physical pH environment, the following steps are performed:

1. PROPKA is used to compute the pKa of each ionizable residue
2. NAMD is used to protonate residues with pKa > pH
3. NAMD is used to relax protons and all side chains
4. NAMD is used to compute energy of each relaxed loop conformation
5. CASTp is used to determine if relaxed loop conformation is in conducting state

Note that we have modified NAMD to facilitate capturing the energy state of each residue.

For scripts demonstrating relaxation and topological gating state assignment, please refer to [scripts/d_relax](scripts/d_relax).

## Additional references

For (streaming) affinity propagation clustering:

* Frey, Brendan J., and Delbert Dueck. "Clustering by passing messages between data points." science 315, no. 5814 (2007): 972-976.
* Zhang, Xiangliang, Cyril Furtlehner, and Michele Sebag. "Data streaming with affinity propagation." Machine learning and knowledge discovery in databases (2008): 628-643.
* Zhang, Xiangliang, Cyril Furtlehner, Cecile Germain-Renaud, and Michele Sebag. "Data stream clustering with affinity propagation." IEEE Transactions on Knowledge and Data Engineering 26, no. 7 (2014): 1644-1656.

For protonation state assignment with PROPKA:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.
* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.

For relaxation and energy scoring with NAMD and CHARMM36 force field:

* Phillips, James C., Rosemary Braun, Wei Wang, James Gumbart, Emad Tajkhorshid, Elizabeth Villa, Christophe Chipot, Robert D. Skeel, Laxmikant Kale, and Klaus Schulten. "Scalable molecular dynamics with NAMD." Journal of computational chemistry 26, no. 16 (2005): 1781-1802.
* Huang, J. and MacKerell, A.D., 2013. CHARMM36 all‐atom additive protein force field: Validation based on comparison to NMR data. Journal of computational chemistry, 34(25), pp.2135-2145.

For topological gating state assignment:

* Dundas, Joe, Zheng Ouyang, Jeffery Tseng, Andrew Binkowski, Yaron Turpaz, and Jie Liang. "CASTp: computed atlas of surface topography of proteins with structural and topographical mapping of functionally annotated residues." Nucleic acids research 34, no. suppl_2 (2006): W116-W118.
