/* residue .h
Class Residue.
*/

#ifndef _RESIDUE_
#define _RESIDUE_

#include "atom.h"

// Data files
#define FILE_ATOMPROP "data/atomProp2.txt"

#define CC_DIS_CUT 12    // residue center to center distance cutoff,

// Backbone atom offsets: N, CA, C, O, H, CB, Gly has pseudo CB and Pro has a pseudo H
enum {
    ATM_N=0,     // 0
    ATM_CA,      // 1
    ATM_C,       // 2
    ATM_O,       // 3
    ATM_H,       // 4
    ATM_CB,      // 5
    NUM_BB_ATOM  // 6
};

// Residue types
enum { 
    ALA=0,       // 0
    CYS,         // 1
    ASP,         // 2
    GLU,         // 3
    PHE,         // 4
    GLY,         // 5
    HIS,         // 6
    ILE,         // 7
    LYS,         // 8
    LEU,         // 9
    MET,         // 10
    ASN,         // 11
    PRO,         // 12
    GLN,         // 13
    ARG,         // 14
    SER,         // 15
    THR,         // 16
    VAL,         // 17
    TRP,         // 18
    TYR,         // 19
    NUM_CANONICAL_RESIDUES // 20
};

// enum type for secondary structure
enum {
    SS_H, SS_E, SS_C
};

// there is no OXT atom type in regular residues, The OXT atom can be added to the last residue
// with proper type and increment its _numAtom by 1.
// For the first N atom, its _type can also be changed.

class Residue {
public:
    Atom *_atom;      // _atom[0] is ca, _atom[1] is cb, etc.
    short _type;      // residue type represented by int
    short _posn;      // from 1 to NumRes
    short _numAtom;   // number atoms in residue, keep track length of _atom
    double _phi;      // phi angle
    double _psi;      //  psi angle
    double _omega;    //  amide torsion angle, omega
    double _scChi[5]; // side chain Chi angles
    // used to label side-chain state, initilaized to -1, -2 means
    // side-chain is fixed.
    int _scState;
    int _pdbIndex;    // index in pdb file

    class Structure *_parent; // point to the parent structure
    // the following attributes are used for fast detection of collisions
    Atom _bbc;        // backbone center including Cb
    Atom _scc;        // side-chain center
    Atom _center;     // center of the whole residue

    // class variable for names
    static string Name1[25];  // one letter name including bases
    static string Name3[25];  // three letter name
    static SIMAP AIMap;  // Amino acid to integer map
    static SSMAP AAMap;  // amino acid three letter name and one letter name map
    static SIMAP SIMap;  // Secondary structure one letter name to integer map
    // map for atoms in each amino acid, for example, N in ASP is DN and
    // its mapped integer value is the index of the atom in atom array _atom
    static SIMAP AtomMap;
    // map for atom index in the atom level distance potential function.
    static SIMAP AtomIndexMap;
    static ISMAP ResMap;
    // parameter between atoms are stored in static variables, the position of each atom is fixed, so that the properties can be located for any atom of a given resiude.
    static string cType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];   // the string of the atom name
    static int vdwType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];    // the integer type of the atom name
    static int prev_atom[NUM_RES_TYPE][MAX_NUM_ATOM_RES][3];
    // Bond length: length from indexed atom to previous atom
    static double bond_length[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
    // Bond angle: angle formed by indexed atom and previous two atoms
    static double bond_angle[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
    static double torsion[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
    static double size[NUM_RES_TYPE];
    static double sc_size[NUM_RES_TYPE];
    static double bb_size;
    static int numAtom[NUM_RES_TYPE];
    // index of functional atoms for each residue
    static IVEC funcAtom[NUM_RES_TYPE];

    explicit Residue(int n = MAX_NUM_ATOM_RES);

    Residue(const Residue &R);

    Residue *operator=(const Residue &r);

    ~Residue();

    void Destruct();

    void init(int n);

    static void InitMap();

    static void InitPar(char *parFile, string &);

    void out();

    void cal_scc();

    /**
     * @return 1-letter AA code
     */
    inline const std::string& get_name_1() const {
        return Name1[this->_type];
    }

    /**
     * @return 3-letter AA code
     */
    inline const std::string& get_name_3() const {
        return Name3[this->_type];
    }

    /**
     * debug function to print atom x,y,z coordinates
     */
    void debug_print() const;
};

#endif
