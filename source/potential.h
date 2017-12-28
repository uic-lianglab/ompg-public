// potential.h
// classes for potential function terms.

#ifndef _POTENTIAL_
#define _POTENTIAL_


#include <string>
#include "reprst.h"
#include "residue.h"

using namespace std;

// Minimum total Boltzmann factor before uniform sampling is triggered
#define MIN_BOLTZF 0.0000000001

// Parameters
#define MAX_NUM_ROT 81
#define MAX_ATOM_TYPE 30

#define LOODIS_DIS_BIN 80      // from 0 - 8 A, 0.1 A per bin
#define PF_DIS_CUT_SQUARE 64
#define H_INLO 0.1             // the interval of loop distance is set to 0.1
#define SC_T_INT 4             // side chain torsion interval in degrees

// this is also used in determine whether two centers are close
#define CUB_SIZE 5.5
#define VDW_CLASH_CUTOFF 0.60
// enum type for energy types
// E_VDWA: atractive part of van der Waals energy
// E_VDWR: repulsive part of van der Waals energy
// E_SOL: solvation
// E_ADP: atom-atom distance potential (statistical, high resolution)
// E_RP: residue pair energy using functional groups FGP.txt
// E_BBT: backbone torsion energy
// E_SCA: side chain atom-atom distance energy, this is used in side
//        chain sampling to store the side chain energy of SIMPL
// E_SCT: side chain torsion energy
// E_HBB: backbone hydrogen bond energy
// E_HBS: sidechain hydrogen bond energy
// E_RAMA: Rama term of Rosetta
// E_SS: secondary structure term
// E_CNT: a fully parameterized VDW
// E_VAA: VDW within the same amino acid, this term affects side chain
// conformations the most.
// E_HAV: h-bond adjusted VDW. For donor and acceptor atom, their VDW
// energy is adjusted
// E_CON:  contact numbers
// E_EP:  loop entropy
// E_SRC:  short-range correlation

enum {
    E_VDWA, E_VDWR, E_VAA, // VDW (0,1,2)
    E_SOL,                 // SOL (3)
    E_HALP,                // HALP (4)
    E_RP,                  // RP (5)
    E_HBB, E_HBS,          // HB (6,7)
    E_BBT,                 // BBT (8)
    E_SCT,                 // SCT (9)
    E_RAMA,                // RAMA H (10)
    E_SIMPL,               // SIMPL (11)
    E_SS,                  // SS (12)
    E_CNT,                 // CNT (13)
    E_ROT,                 // ROT (14)
    E_HAV,                 // HAV (15)
    E_RAMAE,               // RAMA E (16)
    E_RAMAC,               // RAMA C (17)
    E_CON,                 // CON (18)
    E_EP,                  // ENTROPY (19)
    E_SRC,                 // SRC (20)
    E_MODEL,               // MODEL (21)
    E_TEMP,                // Template-based residue-residue distance terms (22)
    E_HB,                  // Hydrogen bond (23)
    E_LOODIS,              // LOODIS (24)
    ENERGY_TYPES
};

// Energy modes
enum {
    EM_VDW, EM_SOL, EM_HALP, EM_RP, EM_HB, EM_BBT, EM_SCT, EM_RAMA, EM_SIMPL,
    EM_SS, EM_CNT, EM_ROT, EM_CON, EM_EP, EM_SRC, EM_MODEL, EM_TEMP, EM_LOODIS,
    ENERGY_MODES
};

class Point;

class Structure;

class Rotamer;

class RotLib;

class SCT;

// all parameters are stored in class PF
class PF {
public:
    // Whether to calculate each energy mode
    static bool cal[ENERGY_MODES];

    static double LOODIS[20][20][LOODIS_DIS_BIN];

    // VDW_AA store for each amino acids the atom pairs that need to be
    // calculated for VDW energy
    static vector<vector<int> > VDW_AA[20];

    // CNT potential terms
    static map<string, IDMAP> Parameter;

    // index of atom pairs for HADIP, the first two dimensions are residues types
    // the second two dimensions are atom types.
    static void InitPar(string);

    // initialize parameters in LOODIS potential
    static void initLOODIS(string filename);
};

#endif
