// main program for folding programs

#include "residue.h"
#include "cal_energy.h"
#include "smc.h"
#include "reprst.h"
#include "potential.h"
#include "util.h"
#include "params.h"
#include "energy_stats_writer.h"
#include "frag_lib_manager.h"
#include "collision_frontend.h"

#include <algorithm>

using namespace std;

/**
 * Utility for getting an output pdb filename
 */
void get_out_pdb_fname(string &out_fname, const Params &params, const int rank) {
    out_fname = params.pdb_outdir + "/" + params.job_prefix + "_" + params.prot_name + "_" + itoa(rank) + ".pdb";
}

/**
 * Utility for getting pdb lib filename
 */
void get_out_pdb_lib_fname(string &out_fname, const Params &params, const Structure &Conf) {
    const string start = itoa(Conf._res[1]._pdbIndex);
    const string end = itoa(Conf._res[Conf._numRes]._pdbIndex);
    out_fname =
        params.pdb_outdir + "/" + params.job_prefix + 
            "_" + params.prot_name + "_" + start + "_" + end + ".lib.pdb";
}

/**
 * Utility parses text file containing loop start, end regions
 */
void parse_loop_intervals(const string &fname,
                          std::vector<int>& starts,
                          std::vector<int>& ends) {
    // Early out if no file is specified
    if (fname.empty()) {
        return;
    }

    // Early out if file not open
    ifstream LOOPLIST;
    LOOPLIST.open(fname.c_str(), ios::in);
    if (!LOOPLIST.good()) {
        std::cout << "Loop interval file not found: " << fname << std::endl;
        return;
    }
    assert(LOOPLIST.is_open());

    // Parse regions
    string tmp;
    vector<string> tmploopend;
    while (std::getline(LOOPLIST, tmp)) {
        split(tmp, ' ', tmploopend);
        if (tmploopend.size() == 1) continue;
        const int start = atoi(tmploopend[0].c_str());
        const int end = atoi(tmploopend[1].c_str());
        tmploopend.clear();
        starts.push_back(start);
        ends.push_back(end);
    }

    // Print parsed regions
    std::cout << "Read file: " << fname << std::endl;
    for (size_t i = 0; i < starts.size(); ++i) {
        std::cout << "START:\t" << starts[i]
                  << "\tEND:\t" << ends[i] << std::endl;
    }
}

/**
 * Parses start and end regions of simulated loop regions from file
 * @TODO - probably should be a member method of params
 */
void init_muloop(Params &params) {
    if (!params.muloop_file.empty()) {
        std::cout << "Reading multi-loop simulation regions.\n";
    }
    parse_loop_intervals(params.muloop_file, params.starts, params.ends);
    if (params.starts.empty() || 
        params.ends.empty() ||
        (params.starts.size() != params.ends.size())) {
        std::cout << "Error reading loop intervals.\n";
        exit(0);
    }
    params.start = params.starts[0];
    params.end = params.ends[0];
}

/**
 * Parses mask files
 * @TODO - probably should be a member method of params
 */
void init_mask(Params &params) {
    if (!params.mask_file.empty()) {
        std::cout << "Reading mask regions.\n";
    }
    parse_loop_intervals(params.mask_file, params.starts_mask, params.ends_mask);
}

/**
 * Map start and end residue pdb identifiers to internal indices
 */
void map_sim_region(const Structure &Conf, std::vector<int> &starts, std::vector<int> &ends, const bool allow_single_res_loop) {

    if (starts.size() != ends.size()) {
        cout << "INVALID LOOP SPECIFICATION - number of start residues does not match number of end residues.\n";
        exit(0);
    }

    if (starts.empty()) {
        return;
    }

    const int INDEX_NOT_FOUND = -1;
    const int n_loops = static_cast<int>(starts.size());
    for (int i_loop = 0; i_loop < n_loops; ++i_loop) {
        const int pdb_start = starts[i_loop];
        const int pdb_end = ends[i_loop];

        // Note pdb_end can equal pdb_start for certain utilities like PostScRefine
        if (pdb_end < pdb_start) {
            cout << "INVALID LOOP SPECIFICATION - end residue must be greater than start residue.\n";
            exit(0);
        }

        if (pdb_start == pdb_end && !allow_single_res_loop) {
            // Currently, this should only be allowed for PostScRefine utility
            cout << "INVALID LOOP SPECIFICATION - end residue must be greater than start residue.\n";
        }

        int mapped_start = INDEX_NOT_FOUND;
        int mapped_end = INDEX_NOT_FOUND;
        for (int i_res = 1; i_res <= Conf._numRes; ++i_res) {
            const Residue &res = Conf._res[i_res];
            if (pdb_start == res._pdbIndex) {
                mapped_start = i_res;
            }

            if (pdb_end == res._pdbIndex) {
                mapped_end = i_res;
                break;
            }
        }

        if (mapped_start == INDEX_NOT_FOUND) {
            cout << "INVALID LOOP SPECIFICATION - start residue " << pdb_start << " not found.\n";
            exit(0);
        }

        if (mapped_end == INDEX_NOT_FOUND) {
            cout << "INVALID LOOP SPECIFICATION - end residue " << pdb_end << " not found.\n";
            exit(0);
        }

        assert(mapped_start <= mapped_end);
        if (mapped_start != pdb_start) {
            cout << "Internally remapping pdb loop start residue " << pdb_start << " to " << mapped_start << ".\n";
            starts[i_loop] = mapped_start;
        }

        if (mapped_end != pdb_end) {
            cout << "Internally remapping pdb loop end residue " << pdb_end << " to " << mapped_end << ".\n";
            ends[i_loop] = mapped_end;
        }
    }
}

/**
 * Patches fragment libraries _posn and possibly other attributes to match
 * protein. Exits program if the fragments do not match underlying conformation
 */
void patch_multiloop_libs(const Structure &Conf, vector<vector<Structure> > &MultiLooplib,
                          const std::vector<int> &sTart,
                          const std::vector<int> &eNd) {

    assert(sTart.size() == eNd.size());
    const int n_loops = static_cast<int>(sTart.size());

    if (MultiLooplib.size() != n_loops) {
        cout << "ERROR - number of fragment libraries does not match number of simulated loops.\n";
        exit(0);
    }

    for (int i_loop = 0; i_loop < n_loops; ++i_loop) {
        const int n_frags = static_cast<int>(MultiLooplib[i_loop].size());
        if (n_frags == 0) {
            // It's possible that some of the loops don't have frag libs, in which case
            // we can ignore them, but unfortunately, downstream code assumes that
            // libs are 1:1 with number of loops so we have to insert dummy libs.
            continue;
        }
        for (int i_conf_res = sTart[i_loop]; i_conf_res <= eNd[i_loop]; ++i_conf_res) {
            const Residue &conf_res = Conf._res[i_conf_res];
            for (int i_frag = 0; i_frag < n_frags; ++i_frag) {
                const int i_frag_res = i_conf_res - sTart[i_loop] + 1;
                Structure &frag = MultiLooplib[i_loop][i_frag];
                if (i_frag_res > frag._numRes) {
                    cout << "ERROR - fragment size is smaller than loop size.\n";
                    exit(0);
                }
                Residue &frag_res = frag._res[i_frag_res];
                if (frag_res._type != conf_res._type) {
                    cout << "ERROR - fragment residue type " << Residue::Name3[frag_res._type]
                    << " does not match conformation type " << Residue::Name3[conf_res._type] << endl;
                    exit(0);
                }
                if (frag_res._numAtom != conf_res._numAtom) {
                    cout << "ERROR - fragment residue atom count does not match conformation residue atom count.\n";
                    exit(0);
                }

                // Patch here:
                frag_res._posn = conf_res._posn;
            }
        }
    }
}

/**
 * Computes clash count of simulated regions
 */
void do_clash_check_only(const Params &params, const Structure &nativeConf) {
    Structure maskedConf = nativeConf;
    // Only mask regions are ignored for clash checks
    maskedConf.mark_loop_regions(params.starts_mask, params.ends_mask);
    // Create a collision interface
    CollisionFrontend cfe;
    cfe.init(maskedConf, params);
    // Check broad phase collisions
    cfe.print_broad_clash_counts(
        nativeConf,
        params.starts,
        params.ends,
        true /*b_ignore_self*/
        );
}

/**
 * Scores parameter PDB
 */
void do_score_only(SMC &smc, const Params &params, Structure &nativeConf) {
    cout << "Scoring native conformation.\n";
    nativeConf._energy = smc.pfe.energy(nativeConf, smc.starts, smc.ends);
    std::vector<Structure*> Topconflist;
    Topconflist.push_back(&nativeConf);
    EnergyStatsWriterScoreOnly().write(Topconflist, params);
}

/**
 * Utility for printing side chain clashes from output of
 * Clash_detection_list() calls. This is *not* same as broad phase clashes
 * detected by do_clash_check_only() utility. In particular, membrane
 * collisions are *not* detected.
 * @return total number of detected clashes
 */
int print_sc_clashes(const Structure &conf, const vector<int> &ResIdx, const vector<int> &ClashNum) {
    assert(ResIdx.size() == ClashNum.size());
    int total_clashes = 0;
    for (size_t i = 0; i < ResIdx.size(); ++i) {
        const int res_idx = ResIdx[i];
        // bounds check
        assert(res_idx >= 1);
        assert(res_idx < conf._numRes);
        const int clash_num = ClashNum[i];
        total_clashes += clash_num;
        const Residue &res = conf._res[res_idx];
        cout << "-Side chain clashes at " << res._pdbIndex << ":" << res.get_name_3() << " = " << clash_num << endl;
    }
    cout << "Total side chain clashes = " << total_clashes << endl;
    return total_clashes;
}

/**
 * Side-chain clash refinement post-processing utility
 */
void do_post_sc_refine(const SMC &smc, const Params &params, Structure &nativeConf) {
    cout << "======================================================\n";
    cout << "Running post processing side chain refinement utility.\n";
    cout << "======================================================\n";

    // Vector stores residues which clash
    vector<int> ResIdx;
    // Vector stores clash counts at each clashing residue
    vector<int> ClashNum;

    cout << "Initial clash check ...\n";
    Clash_detection_list(nativeConf, smc.starts, smc.ends, ResIdx, ClashNum, smc.Reslist);
    int num_clashes = print_sc_clashes(nativeConf, ResIdx, ClashNum);

    for (int i = 0; i < params.post_sc_refine; ++i) {
        // Stop refinement if no clashes detected
        if (num_clashes <= 0) {
            break;
        }

        cout << "Refining side chains (" << i+1 << ") ...\n";
        // Note: vec_start was being passed in for vec_end - assumed this was bug
        SCE_Minimization_list(nativeConf, smc.starts, smc.ends, ResIdx, ClashNum, smc.Reslist, SCE_MIN_DELTA_ROT);
        ResIdx.clear();
        ClashNum.clear();

        cout << "Clash check ...\n";
        Clash_detection_list(nativeConf, smc.starts, smc.ends, ResIdx, ClashNum, smc.Reslist);
        num_clashes = print_sc_clashes(nativeConf, ResIdx, ClashNum);
    }

    // Output refined structure
    if (params.post_sc_refine > 0) {
        const string pdb_path(params.pdb_outdir + "/" + params.prot_name + ".pdb");
        cout << "Writing refined PDB to: " << pdb_path << std::endl;
        nativeConf.writePdb(pdb_path, 1, nativeConf._numRes, 0);
    }
}

/**
 * Perform single loop growth procedure
 */
void do_single_loop_growth(SMC &smc, const Params &params) {
    cout << "Protein Name:\t" << smc.prot_name.substr(0, 4) << endl;
    // Set number of states, populating out to residue length if necessary
    const int FragLength = smc.end - smc.start;
    cout << "Start Residue  " << smc.start << "   :   End Residue  " << smc.end << endl;
    // Launch SMC
    vector<Structure*> Topconflist;
    smc.Wholeproc(Topconflist);
    if (params.output_pdb > 0 && !Topconflist.empty()) {
        const int ConfNum = (params.output_pdb > static_cast<int>(Topconflist.size())) 
            ? static_cast<int>(Topconflist.size()) : params.output_pdb;
        string pdblib_fname = "";
        get_out_pdb_lib_fname(pdblib_fname, params, *(Topconflist[0]));
        cout << "Output generated conformations to " << pdblib_fname << endl;
        for (int i = 0; i < ConfNum; i++) {
            const Structure &top_conf = *(Topconflist[i]);
            top_conf.writePdb(pdblib_fname,
                              1, /* start residue */
                              top_conf._numRes, /* end residue */
                              1, /* don't write side chains */
                              i /* model number */
            );
        }

        EnergyStatsWriterFragLib().write(Topconflist, params);
    }
}

/**
 * Perform multi loop growth procedure
 */
void do_multi_loop_growth(SMC &smc, const Params &params) {
    cout << "Protein Name:\t" << smc.prot_name << "\t" << smc.prot_name.substr(0, 4) << endl;
    vector<Structure> Topconflist;
    vector<loop_info> multiloops;
    loop_info TMP;
    int tmpLength;
    int FragLength = 0;

    for (int i = 0; i < params.starts.size(); i++) {
        TMP.Start = params.starts[i];
        TMP.End = params.ends[i];
        TMP.CurrentPos = params.starts[i];
        multiloops.push_back(TMP);
        tmpLength = params.ends[i] - params.starts[i];
        if (tmpLength > FragLength)
            FragLength = tmpLength;
    }
    
    // Load fragment libraries for each simulated loop region
    if (params.use_multiloop_lib) {
        FragLibManager::load(smc.MultiLooplib, params.mulib_frag_ini_path, params.mulib_frag_dir, params.starts, params.ends);
        patch_multiloop_libs(smc.Conf, smc.MultiLooplib, smc.starts, smc.ends);
    }

    smc.Whole_multiloop(multiloops, Topconflist);
    if (params.output_pdb) {

        // @HACK - Output of whole multiloop should be a vector of pointer to structures sorted by energy
        // However, this a quick patch up for now. I think the better solution would be to have a single
        // Conf structure that would be patched with the current generated loops. Then, this state
        // would have its stats computed and written to disk. Kind of couples things more, but it would
        // be more memory efficient. If memory wasn't an issue, then creating the entire structure is fine.
        // For now, am just converting to a vector of pointers
        vector<Structure*> indirect_Topconflist;

        string pdb_fname = "";
        const int n_loops = (int) Topconflist.size();
        for (int i = 0; i < n_loops; i++) {
            get_out_pdb_fname(pdb_fname, params, i);
            Topconflist[i].writePdb(pdb_fname.c_str(), 1, Topconflist[i]._numRes, 0);
            // @HACK store pointer to full conformation
            indirect_Topconflist.push_back(&Topconflist[i]);
        }

        EnergyStatsWriterMultiLoops().write(indirect_Topconflist, params);
    }
}

/**
 * Main entry point
 */
int main(const int argc, const char *argv[]) {

    // protFile stores protein's coordinates;
    char parFile[50];
    strcpy(parFile, FILE_ATOMPROP);

    // Initialize energy modes and weights: we start off weighting energy types
    // equally on their native scales, and including none of them, and an
    // intercept of zero
    for (int i = 0; i < ENERGY_MODES; i++)
        PF::cal[i] = false;
    PF::cal[EM_LOODIS] = true;

    Params params;
    params.parse_command_line(argc, argv);
    params.print();

    // seed random numbers
    srand(params.rand_seed);

    // initialize class variables
    Atom::InitPar(parFile, params.dir);
    Residue::InitMap();
    Residue::InitPar(parFile, params.dir);
    PF::InitPar(params.param_file);
    if (PF::cal[EM_LOODIS] == true) {
        PF::initLOODIS("data/LOODIS_ed4_8_V3.txt");
    }

    // We're in either fold or grow modes, so we first do a common set of
    // things necessary to initialize the SMC object
    // initialize side chain angles
    SCR::InitSCAng(FILE_SCTORSION2);
    
    // If user defined number of side chain states exceeds MAX_NUM_SC_ST,
    // then we will end up with buffer overflows
    if (MAX_NUM_SC_ST < params.num_sc_states) {
        cout << "Number of side chain sample states "
             << params.num_sc_states << " exceeds max allowed states (" << MAX_NUM_SC_ST << ")" << endl;
        exit(0);
    }

    // Early out if we're just reformatting a fragment library
    if (params.should_export_cl_csv_frag_lib) {
        extern void do_export_cl_csv_frag_lib(const Params &);
        do_export_cl_csv_frag_lib(params);
        exit(0);
    }
    
    // Initialize multi-loop regions
    init_muloop(params);

    // Early out if we're reformatting full PDB sets
    if (params.should_export_cl_csv_pdbs) {
        extern void do_export_cl_csv_pdb(const Params &);
        do_export_cl_csv_pdb(params);
        exit(0);
    }

    // Initialize mask
    init_mask(params);

    // Read the initial structure
    const SSET EMPTY_SET;
    Structure nativeConf(MAX_NUM_RES);
    nativeConf.readPdb(params.prot_file.c_str(), EMPTY_SET);

    // Check mutations residue within sim regions
    params.muts.verify_all_in_sim_region(params.starts, params.ends);

    // Apply any mutations
    params.muts.apply(nativeConf);

    // Map start and end residue PDB identifiers to internal indices
    const bool b_post_sc_refine = params.post_sc_refine >= 0;
    map_sim_region(nativeConf, params.starts, params.ends, b_post_sc_refine);
    map_sim_region(nativeConf, params.starts_mask, params.ends_mask, b_post_sc_refine);

    // Override default protein name (from file) if given explicitly
    if (params.prot_name != "") {
        nativeConf._prot_name = params.prot_name;
    }

    // Early out if we're only checking clash status
    if (params.clash_check_only) {
        do_clash_check_only(params, nativeConf);
        exit(0);
    }

    // Initialize SMC simulation object
    SMC smc(nativeConf, params);
    
    std::cout << "start:\t" << smc.start << "\tend:\t" << smc.end << std::endl;
    if (smc.conf_keep == 1) {
        std::cout << "Warning: only one conformation will be kept." << std::endl;
    }
    std::cout << "Using REDCELL!! Residue number for calculation is reduced from  " 
              << smc.Conf._numRes << " to " << smc.Reslist.size() << endl;

    // @TODO - if incorporating a hybrid statistical potential,
    // then will want to initialize it here or perhaps SMC constructor

    // Only compute energy
    if (params.score_only) {
        do_score_only(smc, params, nativeConf);
    }
    // Side chain refinement
    else if (b_post_sc_refine) {
        do_post_sc_refine(smc, params, nativeConf);
    }
    // Single loop modeling (Fragment Library Generation)
    // @TODO - Address notion that single loop modeling might not just be for libraries
    else if (params.starts.size() == 1) {
        do_single_loop_growth(smc, params);
    }
    else {
         // Multi-loop chain-growth
         do_multi_loop_growth(smc, params);
    }
}
