/**
 * Utilities for converting from internal DiSGro atom types to simple C, H, O, N, S
 */

#include "build.h"
#include "vdw_utils.h"

/**
 * Select van der Waals radii scheme
 */

/**
 * Radii based on values found in atomProp2.txt
 * This is most conservative scheme as values are fairly large overall.
 * Perhaps values for carbon are based off of r_min value in Lennard-
 * Jones potential; in which case, r_vdw ~= r_min / 1.122
 * Tang, Ke et al. "Fast protein loop sampling and structure prediction using distance-guided sequential chain-growth Monte Carlo method." PLoS Comput Biol 10.4 (2014): e1003539.
 */
#define DISGRO_VDW_SCHEME_KE 0

/**
 * Radii based on http ://en.wikipedia.org/wiki/Van_der_Waals_radius
 * which cites : http ://pubs.acs.org/doi/pdf/10.1021/j100785a001
 * Bondi, A. (1964). "Van der Waals Volumes and Radii".J.Phys.Chem. 68 (3) : 441 - 51. doi : 10.1021 / j100785a001.
 */
#define DISGRO_VDW_SCHEME_BONDI 1

/**
 * Radii based on United Atom Radii (ProtOr) which attempt to pad for implicit hydrogens
 * https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwtables.html
 * Tsai, Jerry, et al. "The packing density in proteins: standard radii and volumes." Journal of molecular biology 290.1 (1999): 253-266.
 */
#define DISGRO_VDW_SCHEME_UAR 2

/**
 * Radii based on Pymol tool
 * can be found by entering in pymol console:
 *  iterate (name <atom_type>), print vdw
 * e.g.
 *  iterate (name C), print vdw
 * will print van der Waals radii for all carbons in loaded protein
 */
#define DISGRO_VDW_SCHEME_PYMOL 3

/**
 * Radii based on CHARMM36 forcefield Lennard-Jones potential r_min
 * https://en.wikipedia.org/wiki/Lennard-Jones_potential
 * http://chemwiki.ucdavis.edu/Core/Physical_Chemistry/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Intermolecular_Forces/Specific_Interactions/Lennard-Jones_Potential
 *
 * We can derive van der Waals radii through relation:
 *  r_vdw = (r_min/2)/(2^(1/6) ~= (r_min/2)/1.122 where r_min is not a radius
 *      but rather the distance at which the LJ potential is at a minimum
 *      (to convert to radius, we divide r_min by 2).
 * References:
 *
 * Robert B. Best, R.B., Xiao Zhu, X., Shim, J., Lopes, P.
 * Mittal, J., Feig, M. and MacKerell, A.D., Jr. Optimization of the
 * additive CHARMM all-atom protein force field targeting improved
 * sampling of the backbone phi, psi and sidechain chi1 and chi2
 * dihedral angles. In preparation
 *
 * MacKerell, A.D., Jr., Feig, M. and Brooks, III, C.L. "Improved
 * treatment of the protein backbone in empirical force fields," Journal
 * of the American Chemical Society, 126: 698-699, 2004
 *
 * MacKerell, Jr., A. D.; Bashford, D.; Bellott, M.; Dunbrack Jr., R.L.;
 * Evanseck, J.D.; Field, M.J.; Fischer, S.; Gao, J.; Guo, H.; Ha, S.;
 * Joseph-McCarthy, D.; Kuchnir, L.; Kuczera, K.; Lau, F.T.K.; Mattos,
 * C.; Michnick, S.; Ngo, T.; Nguyen, D.T.; Prodhom, B.; Reiher, III,
 * W.E.; Roux, B.; Schlenkrich, M.; Smith, J.C.; Stote, R.; Straub, J.;
 * Watanabe, M.; Wiorkiewicz-Kuczera, J.; Yin, D.; Karplus, M.  All-atom
 * empirical potential for molecular modeling and dynamics Studies of
 * proteins.  Journal of Physical Chemistry B, 1998, 102, 3586-3616.
 */
#define DISGRO_VDW_SCHEME_CHARMM36 4

/**
 * Scheme where r_vdw = max(r_uar, r_charmm36)
 */
#define DISGRO_VDW_MAX_UAR_CHARMM36 5

// Set VDW scheme here
#define DISGRO_VDW_SCHEME DISGRO_VDW_MAX_UAR_CHARMM36

/**
 * Mapping from DiSGro Atom::_type to simplified AtomFlavor (C, H, O, N, S)
 */
const VdwUtils::AtomFlavor VdwUtils::disgro_atom_type_to_CHONS[DISGRO_MAX_ATOM_TYPES] = {
    AF_UNDEF, // 0 - mapping doesn't appear to be used according to atomProp2.txt
    // For reference - here are the radii with atomProp2.txt
    // as of 9/27/2015. Here is the format:
    // 1st column: Van der Waals radius (angstrom)
    // 2nd through 8th columns: welldepth lk_volume lk_lambda lk_dgfree acceptor donor hbondH
    // 9th column: atom type (this is stored in short Atom::_type)
    AF_C,     //    2.0000, 0.1200, 14.7, 3.5,  0.00,  0 0 0  1  # carbonyl C in Asn and Gln and guanidyl C in Arg
    AF_C,     //    2.0000, 0.1200,  8.3, 3.5, -1.40,  0 0 0  2  # carboxyl C in Asp and Glu
    AF_CA,    //    2.0000, 0.0486, 23.7, 3.5, -0.25,  0 0 0  3  # aliphatic C with one H (Val, Ile, Thr)
    AF_CA,    //    2.0000, 0.1142, 22.4, 3.5,  0.52,  0 0 0  4  # aliphatic C with two H (other residues)
    AF_CA,    //    2.0000, 0.1811, 30.0, 3.5,  1.50,  0 0 0  5  # aliphatic C with three H (Ala)
    AF_C,     //    2.0000, 0.1200, 18.4, 3.5,  0.08,  0 0 0  6  # aromatic ring C (His, Phe, Tyr, Trp)
    AF_N,     //    1.7500, 0.2384,  4.4, 3.5,  -8.9,  0 1 0  7  # N in Trp side-chain
    AF_N,     //    1.7500, 0.2384,  4.4, 3.5,  -4.0,  1 0 0  8  # N in His side-chain
    AF_N,     //    1.7500, 0.2384, 11.2, 3.5,  -10.0, 0 1 0  9  # N in Asn and Gln side-chain
    AF_N,     //    1.7500, 0.2384, 11.2, 6.0,  -20.0, 0 1 0  10 # N in Lys side-chain, N-terminus?
    AF_N,     //    1.7500, 0.2384, 11.2, 6.0,  -11.0, 0 1 0  11 # N in Arg side-chain
    AF_N,     //    1.7500, 0.2384,  0.0, 3.5,  -1.55, 0 1 0  12 # N in Pro backbone
    AF_O,     //    1.5500, 0.1591, 10.8, 3.5,  -6.77, 1 1 0  13 # hydroxyl O in Ser, Thr and Tyr
    AF_O,     //    1.5500, 0.1591, 10.8, 3.5,  -10.0, 1 0 0  14 # carbonyl O in Asn and Gln
    AF_O,     //    1.5500, 0.2100, 10.8, 6.0,  -10.0, 1 0 0  15 # carboyxl O in Asp and Glu
    AF_S,     //    1.9000, 0.1600, 14.7, 3.5,  -4.1,  0 0 0  16 # sulfur in Cys and Met
    AF_N,     //    1.7500, 0.2384,  4.4, 3.5,  -5.0,  0 1 0  17 # backbone N'
    AF_CA,    //    2.0000, 0.0486, 23.7, 3.5,   1.0,  0 0 0  18 # backbone CA
    AF_C,     //    2.0000, 0.1400, 14.7, 3.5,   1.0,  0 0 0  19 # backbone C'
    AF_O,     //    1.5500, 0.1591, 10.8, 3.5,  -5.0,  1 0 0  20 # backbone O'
    AF_P,     //    1.9000, 0.3182, 14.7, 3.5,  -4.1,  0 0 0  21 # nucleic acid P
    AF_H,     //    0.1000, 0.0500,  0.0, 3.5,   0.0,  0 0 1  22 # polar H
    AF_H,     //    0.7000, 0.0500,  0.0, 3.5,   0.0,  0 0 0  23 # nonpolar H
    AF_H,     //    0.7000, 0.0500,  0.0, 3.5,   0.0,  0 0 0  24 # aromatic H
    AF_UNDEF, //    0.1000, 0.0500,  0.0, 3.5,   0.0,  0 0 1  25 # backbone HN
    AF_UNDEF, //    1.4000, 0.0500,  0.0, 3.5,   0.0,  1 1 0  26 # H2O
    AF_UNDEF, //    0.1,	  0,	   0,	0,	   0,  0 0 0  27 # for atoms that serve as a place holder
};

#if (DISGRO_VDW_SCHEME == DISGRO_VDW_SCHEME_KE)

const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
        2.00, // AF_C = sp2 carbon
        2.00, // AF_CA = sp3 carbon
        0.70, // AF_H = hydrogen
        1.55, // AF_O = oxygen
        1.75, // AF_N = nitrogen
        1.90, // AF_S = sulfur
        1.90, // AF_P = phosphorus
        -1.0  // AF_UNDEF = undefined
};

#elif (DISGRO_VDW_SCHEME == DISGRO_VDW_SCHEME_BONDI)

const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
        1.70, // AF_C = sp2 carbon
        1.70, // AF_CA = sp3 carbon
        1.20, // AF_H = hydrogen
        1.52, // AF_O = oxygen
        1.55, // AF_N = nitrogen
        1.80, // AF_S = sulfur
        1.80, // AF_P = phosphorus
        -1.0  // AF_UNDEF = undefined
};

#elif (DISGRO_VDW_SCHEME == DISGRO_VDW_SCHEME_UAR)

const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
        1.61, // AF_C = sp2 carbon
        1.88, // AF_CA = sp3 carbon
        1.20, // AF_H = hydrogen (not actually part of UAR as model assumes implicit H)
        1.46, // AF_O = oxygen
        1.64, // AF_N = nitrogen
        1.77, // AF_S = sulfur
        1.80, // AF_P = phosphorus (not specified in UAR, using Bondi)
        -1.0  // AF_UNDEF = undefined
};

#elif (DISGRO_VDW_SCHEME == DISGRO_VDW_SCHEME_PYMOL)

const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
    1.70, // AF_C = sp2 carbon (Source: Bondi et al)
    1.88, // AF_CA = sp3 carbon (Source: UAR/ProtOr)
    1.20, // AF_H = hydrogen (Source: Bondi et al)
    1.52, // AF_O = oxygen (Source: Bondi et al)
    1.55, // AF_N = nitrogen (Source: Bondi et al)
    1.80, // AF_S = sulfur (specifically SD of methionine) (Source: UAR/ProtOr)
    1.80, // AF_P = phosphorus (Source: Bondi et al)
    -1.0  // AF_UNDEF = undefined
};

#elif (DISGRO_VDW_SCHEME == DISGRO_VDW_SCHEME_CHARMM36)

// Values are from charmm36 params.prm obtained from Mackerell lab
// http://mackerell.umaryland.edu/charmm_ff.shtml
// R_vdw = (Rmin/2) / (2^(1/6)) rounded up to nearest 0.01 angstrom
const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
    1.79, // AF_C = sp2 carbon
    1.79, // AF_CA = sp3 carbon
    1.20, // AF_H = hydrogen - not modeled, so this value was not computed
    1.52, // AF_O = oxygen
    1.65, // AF_N = nitrogen
    1.79, // AF_S = sulfur (specifically SD of methionine)
    1.92, // AF_P = phosphorus (from MacKerell par_all36_lipid.prm)
    -1.0  // AF_UNDEF = undefined
}; 

#elif (DISGRO_VDW_SCHEME == DISGRO_VDW_MAX_UAR_CHARMM36)

const double VdwUtils::atom_flavor_to_vdw_radius[AF_NUM] = {
    1.79, // AF_C = sp2 carbon
    1.88, // AF_CA = sp3 carbon
    1.20, // AF_H = hydrogen - not modeled
    1.52, // AF_O = oxygen
    1.65, // AF_N = nitrogen
    1.80, // AF_S = sulfur (specifically SD of methionine)
    1.92, // AF_P = phosphorus
    -1.0  // AF_UNDEF = undefined
};

#else

#error Undefined van der Waals scheme!

#endif // DISGRO_VDW_SCHEME
