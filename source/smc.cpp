// SMC code

#include <algorithm>
#include "cal_energy.h"
#include "build.h"
#include "smc.h"
#include "sample_states.h"
#include "util.h"
#include "vla.h"
#include "params.h"

#include <assert.h>
#include <cmath>
#include <iostream>

using namespace std;

// @return TRUE if conformation is within clash threshold, FALSE o/w
bool passes_clash_count(const int max_num_clashes,
                        vector<int> &ResIdx,
                        vector<int> &ClashNum,
                        const bool b_reset) {
    assert(max_num_clashes >= 0);
    bool clash_flag = false;
    const int nResIdx = (int) ResIdx.size();
    for (int t = 0; t < nResIdx; ++t) {
        if (ClashNum[t] > max_num_clashes) {
            clash_flag = true;
            break;
        }
    }
    if (b_reset) {
        ResIdx.clear();
        ClashNum.clear();
    }
    return !clash_flag;
}

// Reports number of closed conformations after sampling has finished
void report_sampling_finished(const time_t start_time, const time_t end_time, const SMC &smc) {
    std::cout << "Conformational sampling done with " << smc.NumClosedconf << " closed samples out of ";
    std::cout << smc.num_conf << " attempts (" << (double) smc.NumClosedconf / (double) smc.num_conf << ")\n";
    std::cout << "\ttime cost:  " << (end_time - start_time) << " s" << endl;
}

SMC::SMC(const Structure &nativeConf, const struct Params &params) {
    this->NumClosedconf = 0;

    this->Conf.init(nativeConf._numRes);
    this->_conf.init(nativeConf._numRes);

    //////////////////////////////////////////////////////////////////////////

    // Set start and ending positions for fragment to be regrown (only for smc)
    this->start = (params.starts[0] == -1) ? 1 : params.starts[0];
    this->end = (params.ends[0] == -1) ? nativeConf._numRes : params.ends[0];
    this->Conf = nativeConf;
    
    // Make sure regions which need to be grown are marked in some fashion
    // to avoid inclusion in potential energy calculations
    this->Conf.mark_loop_regions(params.starts, params.ends);
    // Also mark the mask regions in the same fashion
    this->Conf.mark_loop_regions(params.starts_mask, params.ends_mask);

    // @TODO - store params structure and provide read-only accessors
    this->prot_name = params.prot_name;
    this->rand_seed = params.rand_seed;
    this->num_conf = params.num_conf;
    this->ang_type = params.ang_type;
    this->should_close = params.should_close;
    this->num_dist_states = params.num_dist_states;
    this->num_sc_states = params.num_sc_states;
    this->conf_keep = params.conf_keep;
    this->outputconf = params.output_pdb;
    this->max_num_bb_clashes = params.max_num_bb_clashes;
    this->max_num_sc_clashes = params.max_num_sc_clashes;
    this->max_muloop_restart = params.max_muloop_restart;
    this->max_frag_checks = params.max_frag_checks;
    this->use_multiloop_lib = params.use_multiloop_lib;
    this->use_rot_lib = params.use_rot_lib;
    this->starts = params.starts;
    this->ends = params.ends;

    // Initialize collision checks for static geometry
    this->cfe.init(this->Conf, params);
    
    // Initialize ellipsoid criterion
    this->init_REDCELL(params);

    // Initialize conformation to be grown with template
    this->_conf = this->Conf;

    // bbt torsion angle for joint prob sampling
    this->simpBBT_Init("data/BBT_phi_psi_pair_NEW.txt");
    if (params.num_dist_states > 0) {
        // Load fragment end-to-end distance distributions
        // Initialize conditional distributions for C and N selection
        this->fragdis(FILE_FRAG_N_C, 0);
        this->fragdis(FILE_FRAG_C_CA, 1);
        this->geometryinfo(FILE_LOOPGEO);
    }
}

void SMC::Wholeproc(vector<Structure *> &Topconflist) {
    vector<combiEIndex> Erank;
    combiEIndex temp;
    int i_conf = 0;

    Topconflist.clear();

    cout << "DiSGro in progress ..." << endl;
    PreProcess();
    const time_t start_time = time(NULL);
    /********Do SMC*************/
    for (i_conf = 0; i_conf < num_conf; i_conf++) {
        smc();
        this->pfe.on_sample_finish();
        if ((i_conf % DISGRO_SMC_REPORT_INTERVAL) == 0) {
            cout << "Finished loop " << i_conf << ", num closed = " << this->NumClosedconf << "." << endl;
        }
        if (_conf.Success) {
            temp.energy = _conf._energy;
            temp.index = NumClosedconf;
            Erank.push_back(temp);
        }
    }

    const time_t end_time = time(NULL);
    report_sampling_finished(start_time, end_time, *this);

    if (Erank.size()) {

        for (int j = 0; j < Erank.size(); j++)
            Erank[j].index = j;

        const int fsize = std::min<int>(Erank.size(), outputconf);
        sort(Erank.begin(), Erank.end(), sortfun_E);

        for (int j = 0; j < fsize; j++) {
            const int Eidx = Erank[j].index;

            // Save final conformation
            Topconflist.push_back(&(LoopStore[Eidx]));
        }
    }
    cout << "Finished sorting by conformation energy" << endl;
}

void SMC::Whole_multiloop(vector<loop_info> &multiloops, vector<Structure> &Topconflist) {
    vector<combiEIndex> Erank;
    combiEIndex temp;

    const time_t start_time = time(NULL);
    std::cout << "m-DiSGro in progress ..." << endl;

    for (int t = 0; t < num_conf; t++) {
        smc_multiloop(multiloops);
        this->pfe.on_sample_finish();
        if ((t % DISGRO_SMC_REPORT_INTERVAL) == 0) {
            cout << "Finished sample attempt " << t << ", num closed = " << this->NumClosedconf << "." << endl;
        }
        if (_conf.Success) {
            temp.energy = _conf._energy;
            temp.index = NumClosedconf;
            Erank.push_back(temp);
        }
    }
    const time_t end_time = time(NULL);
    report_sampling_finished(start_time, end_time, *this);

    Structure tmpconf;
    int Eidx, base = 0;

    if (Erank.size()) {
        for (int j = 0; j < Erank.size(); j++)
            Erank[j].index = j;

        const int fsize = std::min<int>(static_cast<int>(Erank.size()), outputconf);
        sort(Erank.begin(), Erank.end(), sortfun_E);
        const double minEnergy = Erank[0].energy;

        for (int j = 0; j < fsize; j++) {
            base = 0;
            tmpconf = Conf;
            Eidx = Erank[j].index;
            cout << "Erank    " << Erank[j].index << "   " << Erank[j].energy << endl;
            for (int k = 0; k < multiloops.size(); k++) {
                for (int p = multiloops[k].Start; p <= multiloops[k].End; p++) {
                    tmpconf._res[p] = LoopStore[Eidx]._res[p - multiloops[k].Start + 1 + base];
                }
                base += multiloops[k].End - multiloops[k].Start + 1;
            }
            tmpconf._energy = LoopStore[Eidx]._energy;
            tmpconf._enStats = LoopStore[Eidx]._enStats;
            Topconflist.push_back(tmpconf);

        }
        cout << "MinE     " << minEnergy << endl;
    }
}

/**
 * @return true if loop was restarted, false o/w
 */
bool SMC::smc_restart_loop(const int LoopChosen, vector<loop_info> &multiloops,
                           const vector<int> &multilooplength, const int min_frag_len, const int tail_len) {
    // Check if we've exceeded our restart budget
    if ((++num_muloop_restart) <= max_muloop_restart) {
        // Defer to collision module
        this->cfe.restart_loop(LoopChosen,
            multiloops,
            multilooplength,
            min_frag_len,
            tail_len,
            *this);
        return true;
    }
    return false;
}

void SMC::smc_multiloop(vector<loop_info> &multiloops) {

    const int multisize = static_cast<int>(multiloops.size());
    assert(multisize > 0);

    bool grow_Success = false;
    int total_len = 0;
    std::vector<int> loop_ids(multisize);
    std::vector<int> multilooplength(multisize);
    vector<int> Closed(multisize, 0);

    for (int i = 0; i < multisize; i++) {
        loop_ids[i] = i;
        multilooplength[i] = (multiloops[i].End - multiloops[i].Start + 1);
        total_len += multiloops[i].End - multiloops[i].Start + 1;
        multiloops[i].CurrentPos = multiloops[i].Start;
    }

    Structure tmploop(total_len);

    _conf = Conf;
    
    // reset loop restart budget
    this->num_muloop_restart = 0;

    while (loop_ids.size() > 0) {
        // Randomly select a loop to sequentially grow:
        const int select_loop_idx = rand() % loop_ids.size();
        assert(select_loop_idx >= 0 && select_loop_idx < loop_ids.size());
        int LoopChosen = loop_ids[select_loop_idx]; // evenly
        start = multiloops[LoopChosen].CurrentPos;
        end = multiloops[LoopChosen].End;

        // Case - only 3 residues remain to be grown -> attempt closure
        if (end <= start + 2) {
            // Attempt analytic closure
            _conf.analyticClosure(end - 2, start, end, Reslist, true /* Ellipsoid */);
            Closed[LoopChosen] = _conf.IsClosed(end - 2);

            // Case: initial analytic closure failed -> try multiple variations to see if they close
            if (!Closed[LoopChosen]) {
                for (int i = 0; i < 300; i++) {
                    _conf.analyticClosure_h(end - 2, 0.06, PI / 32, PI / 32, start, end, Reslist,
                                            true /* Ellipsoid */);
                    Closed[LoopChosen] = _conf.IsClosed(end - 2);
                    if (Closed[LoopChosen])
                        break;
                }
            }

            // Case: analytic closure successful -> remove loop from further processing for this sample
            if (Closed[LoopChosen]) {
                this->pfe.on_analytic_closure(_conf, LoopChosen, multiloops[LoopChosen].End);
                loop_ids.erase(loop_ids.begin() + select_loop_idx);
            } // Case: analytic closure failed -> restart loop from fragment library (or from beginning)
            else {
                if (!this->smc_restart_loop(LoopChosen, multiloops, multilooplength, 6 /*min_frag_len*/, 3 /*tail_len*/)) {
                    break;
                }
            }

        } // Case - attempt to sequentially grow a random residue
        else {
            grow_Success = grow_one(start, end, LoopChosen);

            // Case: grow failed -> restart loop
            if (!grow_Success) {
                if (!this->smc_restart_loop(LoopChosen, multiloops, multilooplength, 7 /*min_frag_len*/, 4 /*tail_len*/)) {
                    break;
                }
            } // Case: grow successful -> increment growth counter
            else {
                ++(multiloops[LoopChosen].CurrentPos);
            }
        }

    }

    // If loops were closed -> attempt clash resolution
    // @TODO - probably want to add Lennard-Jones potential
    _conf.Success = false;
    if (loop_ids.size() == 0) {

        ///////////////////////////////////////////////////////
        // Backbone clash detection
        ///////////////////////////////////////////////////////

        vector<int> ResIdx;
        vector<int> ClashNum;
        BBClash_detection_list(_conf, starts, ends, ResIdx, ClashNum, Reslist);

        _conf.Success = passes_clash_count(this->max_num_bb_clashes, ResIdx, ClashNum, true /*reset*/);

        if (_conf.Success) {

            ///////////////////////////////////////////////////////
            // Add final set of side chains
            ///////////////////////////////////////////////////////
            if (this->num_sc_states > 0) {
                for (int i = 0; i < multisize; i++) {
                    assert(starts[i] < ends[i]);
                    assert(ends[i] - 2 >= starts[i]);
#ifdef DISGRO_BUILD_ENABLE_SC_GROW_ONE
                    // Only place side chains for the analytic closure interval
                    // as grow_one_with_sc has already placed preceding side chains
                    _conf.grow_sc(ends[i] - 2, ends[i], 0, num_sc_states, ang_type, _conf._toBeSampled, *this,
                                  i /*LoopChosen*/);
#else
                    // Must place side chains along entire interval as grow_one has not been
                    // incrementally adding them
                    _conf.grow_sc(starts[i], ends[i], 0, num_sc_states, ang_type, _conf._toBeSampled, *this,
                                  i /*LoopChosen*/);
#endif // DISGRO_BUILD_ENABLE_SC_GROW_ONE
                }

                Clash_detection_list(_conf, starts, ends, ResIdx, ClashNum, Reslist);
                // Note: vec_start was being passed in for vec_end - assumed this was bug
                SCE_Minimization_list(_conf, starts, ends, ResIdx, ClashNum, Reslist, SCE_MIN_DELTA_ROT);
                // @TODO - investigate if repeated calls to clash detection followed by
                //  minimization are necessary
                // e.g while (TotalClashNum > max_clash) { Clash_detection(...)*; SCE_Minimization(...)*; }

                ResIdx.clear();
                ClashNum.clear();
                Clash_detection_list(_conf, starts, ends, ResIdx, ClashNum, Reslist);
                _conf.Success = passes_clash_count(this->max_num_sc_clashes, ResIdx, ClashNum, true /*reset*/);
            }

            ///////////////////////////////////////////////////////
            // Add to loop store
            ///////////////////////////////////////////////////////
            if (_conf.Success) {
                _conf._energy = this->pfe.energy(_conf, starts, ends);
                assert(std::all_of(Closed.cbegin(), Closed.cend(), [](int i) { return i == 1; }));
                int ResCount = 1;
                for (int i = 0; i < multisize; i++) {
                    for (int j = multiloops[i].Start; j <= multiloops[i].End; j++) {
                        tmploop._res[ResCount] = _conf._res[j];
                        ResCount++;
                    }
                }
                tmploop._energy = _conf._energy;
                tmploop._enStats = _conf._enStats;
                LoopStore.push_back(tmploop);
                ++NumClosedconf;
            }

        } // end check for !BBclash_flag
    } // end check for loop_ids.size() == 0
}

// sequential monte carlo
void SMC::smc() {
    Structure tmploop(end - start + 1);
    _conf = Conf;
    _conf._energy = 0;
    bool prev_Close = should_close;
    bool grow_Success = false, grow_preSuccess = false;  // label the residue successfully generated or not

    _conf.Closed = false;
    _conf.Success = false;
    if (end == _conf._numRes) {
        int i;
        for (i = start; i < end; i++) {
            grow_Success = grow_one(i, end, 0 /*LoopChosen*/);

            if (!grow_Success) {
                _conf.Success = false;
                return;
            }
        }

    }
    else {
        for (int i = start; i < end - 2; i++) {
            grow_Success = grow_one(i, end, 0 /*LoopChosen*/);

            if (!grow_Success) {
                _conf.Success = false;
                return;
            }
        }
    }
    /********** Set closure indicators *************/
    Point next_atom[3];
    if (end - start < 4) {
        grow_preSuccess = grow_one_bb_only(_conf, end - 2, end, 0 /*LoopChosen*/, _conf._res[end]._atom[ATM_C]);
        if (grow_preSuccess) {
            grow_Success = grow_one_bb_only(_conf, end - 1, end, 0 /*LoopChosen*/, _conf._res[end]._atom[ATM_C]);

            next_atom[0] = _conf._res[end + 1]._atom[ATM_CA];
            next_atom[1] = _conf._res[end + 1]._atom[ATM_N];
            next_atom[2] = _conf._res[end]._atom[ATM_C];
            calCo(next_atom,
                  Residue::bond_length[_conf._res[end]._type][ATM_C],
                  Residue::bond_angle[_conf._res[end + 1]._type][ATM_N], PI, _conf._res[end]._atom[ATM_CA]);
        }

        if (grow_Success && grow_preSuccess)
            _conf.Closed = _conf.IsClosed(end);
    }
    if (should_close && (_conf.Closed == false)) {
        // The algorithm assumes the final CA is placed in accordance
        // with the following atoms, which grow_one_bb_only does not guarantee
        // currently. Here we place it using torsion angle omega = PI.
        _conf.analyticClosure(end - 2, start, end, Reslist, true /* Ellipsoid */);
    }

    /********** Set closure indicators *************/
    if (should_close && (end != Conf._numRes)) {
        if (grow_Success && grow_preSuccess) {
            _conf.Closed = _conf.IsClosed(end);
        }
        else if (!grow_Success && grow_preSuccess) {
            _conf.Closed = _conf.IsClosed(end - 1);
        }
        else {
            _conf.Closed = _conf.IsClosed(end - 2);

            if (!_conf.Closed) {
                for (int i = 0; i < 300; i++) {
                    _conf.analyticClosure_h(end - 2, 0.06, PI / 32, PI / 32, start, end, Reslist,
                                            true /* Ellipsoid */);
                    _conf.Closed = _conf.IsClosed(end - 2);
                    if (_conf.Closed)
                        break;
                }
            }
        }
    }

    if (end == Conf._numRes) {
        // Replace original values of overwritten members
        should_close = prev_Close;
    }

    /************** If End = final residue, do some special processing  ***********************/
    if (end == Conf._numRes) {
        // Replace original values of overwritten members
        double rPhiAng = frand((-PI), PI);
        double rPsiAng = frand((-PI), PI);
        Point prev_atom[3];
        Residue &ttRes = _conf._res[end];  // the last residue
        // place final ATM_C of chain
        prev_atom[0] = _conf._res[end - 1]._atom[ATM_C];
        prev_atom[1] = _conf._res[end]._atom[ATM_N];
        prev_atom[2] = _conf._res[end]._atom[ATM_CA];
        calCo(prev_atom,
              Residue::bond_length[_conf._res[end]._type][ATM_C],
              Residue::bond_angle[_conf._res[end]._type][ATM_C],
              rPhiAng,
              _conf._res[end]._atom[ATM_C]);
        // place ATM_O atom on ATM_C
        prev_atom[0] = _conf._res[end]._atom[ATM_N];
        prev_atom[1] = _conf._res[end]._atom[ATM_CA];
        prev_atom[2] = _conf._res[end]._atom[ATM_C];
        calCo(prev_atom,
              Residue::bond_length[_conf._res[end]._type][ATM_O],
              Residue::bond_angle[_conf._res[end]._type][ATM_O],
              rPsiAng,
              _conf._res[end]._atom[ATM_O]);

        // place ATM_OXT on ATM_C
        prev_atom[0] = _conf._res[end]._atom[ATM_N];
        prev_atom[1] = _conf._res[end]._atom[ATM_CA];
        prev_atom[2] = _conf._res[end]._atom[ATM_C];
        calCo(prev_atom,
              Residue::bond_length[_conf._res[end]._type][ATM_O],
              Residue::bond_angle[_conf._res[end]._type][ATM_O],
              rPsiAng + PI,
              _conf._ATM_OXT);

        // sample side chains
        if (ttRes._type == GLY) {
            ttRes._numAtom = 6;   // 5 backbone atoms plus OXT
        }
        else {
            // place ATM_CB on ATM_CA
            prev_atom[0] = _conf._res[end]._atom[ATM_N];
            prev_atom[1] = _conf._res[end]._atom[ATM_C];
            prev_atom[2] = _conf._res[end]._atom[ATM_CA];
            calCo(prev_atom,
                  Residue::bond_length[_conf._res[end]._type][ATM_CB],
                  Residue::bond_angle[_conf._res[end]._type][ATM_CB],
                  PI * 122.55 / 180,
                  _conf._res[end]._atom[ATM_CB]);
        }
        // copy ATM_OXT to the last residue. This copy needs to be done here
        ttRes._atom[ttRes._numAtom] = _conf._ATM_OXT;
    }

    _conf.Success = false;
    if (_conf.Closed || (end == Conf._numRes)) {
        this->pfe.on_analytic_closure(_conf, 0 /*LoopChosen*/, end);

        ///////////////////////////////////////////////////////
        // Backbone clash detection
        ///////////////////////////////////////////////////////

        vector<int> ResIdx;
        vector<int> ClashNum;
        BBClash_detection_list(_conf, start, end, ResIdx, ClashNum, Reslist);

        _conf.Success = passes_clash_count(this->max_num_bb_clashes, ResIdx, ClashNum, true /*reset*/);

        if (_conf.Success) {

            ///////////////////////////////////////////////////////
            // Add final set of side chains
            ///////////////////////////////////////////////////////
            if (this->num_sc_states > 0) {
#if DISGRO_BUILD_ENABLE_SC_GROW_ONE
                // Add side chains for analytic closure interval
                _conf.grow_sc(end - 2, end, 0, num_sc_states, ang_type, _conf._toBeSampled, *this, 0 /*LoopChosen*/);
#else
                // Since grow_one has not placed any side chains, add them all here
                _conf.grow_sc(start, end, 0, num_sc_states, ang_type, _conf._toBeSampled, *this, 0 /*LoopChosen*/);
#endif // DISGRO_BUILD_ENABLE_SC_GROW_ONE

                Clash_detection_list(_conf, start, end, ResIdx, ClashNum, Reslist);
                // Note: vec_start was being passed in for vec_end - assumed this was bug
                SCE_Minimization_list(_conf, start, end, ResIdx, ClashNum, Reslist, SCE_MIN_DELTA_ROT);
                // @TODO - investigate if repeated calls to clash detection followed by
                //  minimization are necessary
                // e.g while (TotalClashNum > max_clash) { Clash_detection(...)*; SCE_Minimization(...)*; }

                // For now - discard sample if clash threshold is exceeded
                ResIdx.clear();
                ClashNum.clear();
                Clash_detection_list(_conf, start, end, ResIdx, ClashNum, Reslist);
                _conf.Success = passes_clash_count(this->max_num_sc_clashes, ResIdx, ClashNum, true /*reset*/);
            }

            ///////////////////////////////////////////////////////
            // Add to loop store
            ///////////////////////////////////////////////////////
            if (_conf.Success) {
                _conf._energy = 0.0;
                _conf.calCenter(start, end, true);
                _conf._energy = this->pfe.energy(_conf, this->starts, this->ends);
                for (int i = start; i < end + 1; i++) {
                    tmploop._res[i - start + 1] = _conf._res[i];
                }
                tmploop._energy = _conf._energy;
                tmploop._enStats = _conf._enStats;
                LoopStore.push_back(tmploop);
                ++NumClosedconf;
            }

        } // end check for !BBclash_flag
    } // end check that analytic closure was successful

}

/**
 * Initializes bounding ellipsoid - smc object must have valid template Conf member
 */
void SMC::init_REDCELL(const struct Params &params) {
    // calculate Ellipsoid constant 2a
    vector<double> EllipConst2A;
    for (int j = 0; j < params.starts.size(); j++)
        EllipConst2A.push_back(Ellipsoid_Detect(this->Conf, params.starts[j], params.ends[j], true));

    this->Reslist.clear();
    // test the space after using ellipsoid and surface
    for (int i = 1; i <= this->Conf._numRes; i++) {
        // Assuming this residue is part of simulated loop region, always include it!
        if (this->Conf._res[i]._center.is_at_origin()) {
            this->Reslist.push_back(i);
            continue;
        }

        for (int j = 0; j < params.starts.size(); j++) {
            const double dis = this->Conf._res[i]._center.dis(this->Conf._res[params.starts[j]]._atom[ATM_CA]) +
                               this->Conf._res[i]._center.dis(this->Conf._res[params.ends[j]]._atom[ATM_C]);
            if (dis < EllipConst2A[j] + Residue::size[this->Conf._res[i]._type] + CUB_SIZE) {
                this->Reslist.push_back(i);
                break;
            }
        }
    }
}

void SMC::fragdis(string fname, int label) {
    // Populate etedP collection of distance lengths
    // Conditional distribution
    ifstream input(fname.c_str());
    int Countset = 0, FraglenCount = -1;
    while (!input.eof()) {
        // Read line
        string line;
        getline(input, line);
        if (line[0] == '#') continue;
        // Tokenize by tabs
        vector<string> token;
        split(line, ' ', token);
        if (token.size() == 4) {
            FraglenCount++;
            Countset = 0;
            minDistcon[label][FraglenCount] = atof(token[0].c_str());
            DistconBy[label][FraglenCount] = atof(token[1].c_str());
            minDistdel[label][FraglenCount] = atof(token[2].c_str());
            DistdelBy[label][FraglenCount] = atof(token[3].c_str());
        }
        else if (token.size() == 32) {
            for (int i = 0; i < 32; i++) {
                etedCon[label][FraglenCount][Countset][i] = atof(token[i].c_str());
            }
            Countset++;
        }
    }
}

void SMC::geometryinfo(string fname) {
    // Populate etedP collection of distance lengths
    // Conditional distribution
    ifstream input(fname.c_str());
    int Countset = 0, FraglenCount = -1;
    while (!input.eof()) {
        // Read line
        string line;
        getline(input, line);
        if (line[0] == '#') continue;

        // Tokenize by tabs
        vector<string> token;
        split(line, ' ', token);
        if (token.size() == 1) {
            FraglenCount++;
            Countset = 0;
            EEdisBy[FraglenCount] = atof(token[0].c_str());
        }
        else if (token.size() == 37) {
            for (int i = 0; i < 37; i++) {
                MtorsionEEdis[FraglenCount][Countset][i] = atof(token[i].c_str());
            }
            Countset++;

        }
    }
}

void SMC::simpBBT_Init(string fname) {
    ifstream input(fname.c_str());
    if (!input.is_open()) {
        cout << "Error!! Cannot open BBT sampling file!!   " << fname << endl;
        exit(0);
    }
    int resType, phi_angle, psi_angle;
    while (!input.eof()) {
        // Read line
        string line;
        getline(input, line);
        if (line[0] == '#') continue;
        if (strlen((char *) line.c_str()) == 0) break;
        // Tokenize by tabs
        vector<string> token;
        split(line, '\t', token);
        if (token.size() == 4) {
            resType = atoi(token[0].c_str());
            phi_angle = atoi(token[1].c_str()) + 180;
            psi_angle = atoi(token[2].c_str()) + 180;
            Joint_Angle[resType][phi_angle / BBTbinSize][psi_angle / BBTbinSize] = atoi(token[3].c_str());
        }
    }
}

// Pre decide some of the useful parameters, e.g. Close
void SMC::PreProcess() {
    // Edit end fragment
    if (end == Conf._numRes) {
        cout << "Editing end fragment! " << endl;
        // If we're editing a start or end fragment, we switch distance states
        // to angle states
        // APR - removed functionality 9/16/15s
        cout << "This is currently unsupported." << endl;

        exit(0);
    }
}
