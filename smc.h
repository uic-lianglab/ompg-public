/* smc.h
header file for smc.cpp, sequential Monte Carlo functions
*/
#ifndef _SMC_
#define _SMC_

#include "structure.h"
#include "rotamer.h"
#include "potential.h"
#include "reprst.h"
#include "util.h"
#include "cal_energy.h"
#include "sample_states.h"
#include "misc_structs.h"
#include "potential_frontend_loodis.h"
#include "collision_frontend.h"

// Empirical distance files
#define FILE_FRAG_C_CA "data/frag.C.CA_pdf_32_19.txt"
#define FILE_FRAG_N_C "data/frag.N.C_pdf_32_19.txt"
#define FILE_LOOPGEO "data/LoopGeo_37_pdf_21.txt"

using namespace std;

class SMC {
public:

    // Totals
    int num_conf; // number of conformations in SMC simulation

    // The template conformation
    Structure Conf;
    // The current conformation being grown
    Structure _conf;
    vector<Structure> LoopStore;
    vector<vector<Structure> > MultiLooplib;
    // Some parameters for SMC sampling:
    // number of distance states sampled for backbone growth
    int num_dist_states;
    // number of states for side chain sampling.
    // if numSCStates=5, then 5 side chain conformations will be sampled for each residue
    int num_sc_states;
    // random seed
    unsigned int rand_seed;
    // name of the protein
    string prot_name;
    // torsion angle type
    int ang_type;
    // Will analytic closure by performed?
    bool should_close;
    // the starting position of the protein sequence for SMC sampling
    // for loop modeling, the C atom of the start residue will be
    // sampled. When sampling start residue, the atoms N and CA of start+1
    // are also sampled.
    int start;
    // The end position of the protein sequence for SMC sampling
    // for loop modeling, the N and CA atoms of end residue are sampled
    // when sampling residue end-1
    int end;

    bool use_multiloop_lib;
    bool use_rot_lib;
    vector<int> Reslist;
    int Joint_Angle[20][TORBIN][TORBIN];
    double etedCon[2][20][32][32];
    double minDistcon[2][20];
    double minDistdel[2][20];
    double DistconBy[2][20];
    double DistdelBy[2][20];

    double MtorsionEEdis[17][37][37]; // The dihedral angle between the plane of middle point of the loop and the plane extended protein body, start from length = 4
    double EEdisBy[17];

    int conf_keep; // the number of kept conformations sorted by energy
    int NumClosedconf;
    int outputconf;
    int max_num_bb_clashes;
    int max_num_sc_clashes;
    int max_muloop_restart;
    int num_muloop_restart;
    int max_frag_checks;

    // parallel vectors for start and end residue positions for simulated loop regions
    vector<int> starts;
    vector<int> ends;

    // Front end interface for potential calculations
    PotentialFrontend pfe;

    // Front end interface for collision checking
    CollisionFrontend cfe;

    SMC(const Structure &nativeConf, const struct Params &params);

    // overall Sequential Monte Carlo algorithm
    void smc();

    void Whole_multiloop(vector<loop_info> &multiloops, vector<Structure> &Topconflist);

    void smc_multiloop(vector<loop_info> &multiloops);

    bool smc_restart_loop(const int LoopChosen, vector<loop_info> &multiloops,
                          const vector<int> &multilooplength, const int min_frag_len, const int tail_len);

    // For glycine and alanine, only backbone atoms are grown (CB is backbone), for all other
    // residue types, side chain atoms are also grown if enabled
    bool grow_one(const int Position, const int tmpEnd, const int LoopChosen);

    // grow residue back bone atoms (no side chains)
    bool grow_one_bb_only(Structure &Conf, const int Position, const int tmpEnd, const int LoopChosen,
                          const Atom &EndPt);

    // grow entire residue including side chains
    bool grow_one_with_sc(Structure &Conf, const int Position, const int tmpEnd, const int LoopChosen,
                          const Atom &EndPt);

    // ellipsoid criterion
    void init_REDCELL(const struct Params &params);

    // for loop modeling
    void fragdis(string fname, int label);

    void geometryinfo(string fname);

    void PreProcess();

    void sample_distance(const Atom &b, const Atom &c, const Atom &B,
                         const double theta, const double lcd,
                         Atom &p, const double lcon, const int label, const int rem);

    void Wholeproc(vector<Structure*> &Topconflist);

    void simpBBT_Init(string fname);

};

#endif
