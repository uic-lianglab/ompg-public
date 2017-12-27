# PROPKA 3.1

PROPKA predicts the pKa values of ionizable groups in proteins
(version 3.0) and protein-ligand complexes (version 3.1)
based on the 3D structure.

For proteins without ligands both version should produce the same result.

The method is described in the following papers, which you should cite
in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.

See [propka.ki.ku.dk](http://propka.ki.ku.dk/) for the PROPKA web server,
using the [tutorial](http://propka.ki.ku.dk/~luca/wiki/index.php/PROPKA_3.1_Tutorial).

## Modifications 

This release of PROPKA 3.1 was modified by Oliver Beckstein
<oliver.beckstein@asu.edu> from the released version.

* Included patches from
  https://github.com/schrodinger/propka-3.1/tree/python27-compat to
  make it compatible with Python 2.7

* Packaged for installation with setuptools.

## Examples

Calculate using pdb file

    propka31 1hpx.pdb


## Testing (for developers)

Please run `Tests/runtest.py/` after changes before pushing commits.

## References / Citations

Please cite these references in publications:

* Sondergaard, Chresten R., Mats HM Olsson, Michal Rostkowski, and Jan H. Jensen. "Improved Treatment of Ligands and Coupling Effects in Empirical Calculation and Rationalization of pKa Values." Journal of Chemical Theory and Computation 7, no. 7 (2011): 2284-2295.

* Olsson, Mats HM, Chresten R. Sondergaard, Michal Rostkowski, and Jan H. Jensen. "PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions." Journal of Chemical Theory and Computation 7, no. 2 (2011): 525-537.



