// structure.h
/**
 * header file for protein structure. It contains definition for Point, Atom, Residue and Structure.
 */


#ifndef _STRUCTURE_
#define _STRUCTURE_

#include "residue.h"
#include "potential.h"
#include "energy_stats.h"

#define MAX_NUM_SC_ST 64

class Structure {
public:
    int _numRes;
    // residues of the conformation, start from position 1.
    Residue *_res;
    // number of chains in the structure
    int _numChain;
    // store the names of the chain, from the first character of the
    // string, assuming each chain is represented by a single character as
    // defined in PDB format.
    string _chainName;
    // Name of the chain
    string _prot_name;
    // overall energy calculated
    double _energy;
    // Detailed energy information
    EnergyStats _enStats;
    // residues to be sampled
    vector<int> _toBeSampled;
    // true if the structure has no backbone gaps
    bool Closed;
    bool Success;

    // the oxygen atom at the end of the chain. It can copied to the last
    // atom of the last residue
    Atom _ATM_OXT;
    Atom _Confcenter;

    // the sequence of the structure
    static string _sequence;
    vector<string> missSeq;
    vector<int> missSeqPos;
    vector<int> missSeqfrom1;

    // Key functions
    explicit Structure(int n = 1);

    Structure(const Structure &S);

    Structure *operator=(const Structure &S);

    ~Structure();

    void Destruct();

    void init(int n);

    void copyStructure(const Structure &C, int start, int end);

    bool readPdb(const string &, const SSET &SelRes);

    bool readPdb(const vector<string> &, const SSET &SelRes);

    void calCenter(int Start, int End);

    void calCenter(int Start, int End, bool Sidechain);

    void writePdb(string filename, int start, int end, int type) const;

    void writePdb(string filename, int start, int end, int type, int md_num) const;

    void outputAngles(char *filename, int start, int end, int outType);

    void calBBCo(const int resIndex, Residue &res, const double phi, const double psi, const double omega) const;

    void calSCCo(const double *torAngles, Residue &res) const;

    void analyticClosure(const int start, const int Start, const int End, const std::vector<int> &List, const bool Ellipsoid);

    void analyticClosure_h(const int start, const double len_change, const double bon_change, const double t_change,
                           const int Start, const int End, const std::vector<int> &List, const bool Ellipsoid);

    bool IsClosed(const int End);

    void grow_sc(const int start, const int end, const int type, const int numStates, int growType,
                 vector<int> &resToGrow, class SMC &smc, const int LoopChosen);

    void SinglecalCenter(Residue &_res, int type) const;

    void StoreSequence();

    // Marks atoms at start residue as being at origin
    // but avoids marking N and CA atoms as they are given
    void mark_start_residue(const int start);

    // Marks atoms at end residue as being at origin
    // but avoids marking CA and C atoms as they are given
    void mark_end_residue(const int end);

    // Move atoms within a single loop region to origin. This "marks" atoms that have not yet
    //  been grown and avoids their use in potential energy calculations.
    void mark_loop_region(const int start, const int end);

    // Move atoms in loop regions to origin (0,0,0). This "marks" atoms that have not yet
    //  been grown and avoids their use in potential energy calculations.
    void mark_loop_regions(const std::vector<int> &sTart, const std::vector<int> &eNd);

    // Debug function for printing status of backbone atoms in interval [start, end]
    void print_bb_atoms_status(const int start, const int end);
};

#endif
