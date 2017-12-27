//// cal_energy.cpp
// calculate the potential energy of a given structure

#include "cal_energy.h"
#include "vdw_utils.h"

#include <algorithm>

bool check_loodis_bounds(const int dim1, const int dim2, const int dim3) {
    bool is_valid = dim1 >= 0 && dim1 < 20;
    assert(is_valid);
    is_valid &= dim2 >= 0 && dim2 < 20;
    assert(is_valid);
    is_valid &= dim3 >= 0 && dim3 < LOODIS_DIS_BIN;
    assert(is_valid);
    return is_valid;
}

// function for calculating energy of one residue with the rest of the protein
// when type==0, it calculate energy of this residue with all atoms of
// all residues between start and end
// when type==1, it is for fragment closure of fress, calculating the
// energy of this residue with residues before the fragment
// when type==2, it is for fragment closure of fress, calculating the
// energy of this residue with residues after the fragment
// this function cannot be used to grow the last two residues of the
// fragment, which is done in two_res_en()
// -------------------------
// Note: not all energy terms are used in growing backbone atoms of fragments.
// It may not be that useful to use all the energy terms. A fragment
// needs to be re-evaluated after backbones have been fully grown and
// side chains placed with all the energy terms.  Here VDW terms are
// important since sterics/collision has to be considered during
// backbone growth.
// -------------------------
double one_res_en_loodis_bb2all(const Structure &conf, const Residue &res, const int start, const int end,
                                const int Start, const int End, const int type) {
    int i, j, k, position, atomNum, disInd;
    double energy = 0, dis, disquare;

    // get the position of res to check for adjacent residues
    position = res._posn;
    for (i = start; i <= end; i++) {
        if (conf._res[i]._center.x == 0 && conf._res[i]._center.y == 0 && conf._res[i]._center.z == 0) continue;
        if ((type == 0 && res._center.dis(conf._res[i]._center) <
                          Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
            (type != 0 &&
             res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {

            if (i >= Start && i <= End)
                atomNum = 6;
            else
                atomNum = conf._res[i]._numAtom;
            for (j = 0; j < NUM_BB_ATOM; j++) {
                // skip undefined and H atoms
                if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                    continue;
                if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                //if(type != 0 && j >= NUM_BB_ATOM ) continue;  // when sampling only backbone atoms, skip side chain atoms if there are any
                for (k = 0; k < atomNum; ++k) {
                    if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                    if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 && conf._res[i]._atom[k].z == 0) {
                        continue;
                    }
                    // when calculating energy for the chain before the fragment
                    // do not calculate the energy of C and O of this residue with C and O of residue at end position

                    if (type == 1) {
                        // adjacent residues, do not calculate VDW for atoms
                        // separated by three bonds or less
                        if (position - i == 1) {
                            if (k == ATM_C && (j == ATM_C || j == ATM_CB))
                                continue;
                        }
                    }

                    if (position > Start) if (type == 1 && i == end && (k == ATM_C || k == ATM_O) &&
                                              (j == ATM_C || j == ATM_O))
                        continue;
                    disquare = (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x) +
                               (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y) +
                               (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);

                    if (disquare <= PF_DIS_CUT_SQUARE) {
                        dis = sqrt(disquare);
                        disInd = (int) (dis / H_INLO);
                        assert(check_loodis_bounds(res._atom[j]._type - 1, conf._res[i]._atom[k]._type - 1, disInd));
                        energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type - 1][disInd];
                    }
                }
            }
        }
    }

    return energy;
}

// type 0 is used for loop closure calculation , type 1 is used for normal chain-growth
double one_res_en_loodis_bb2all_list(const Structure &conf, const Residue &res,
                                     const int Start, const int End, const int type,
                                     const std::vector<int> &List) {
    int p, i, j, k, position;
    int disInd;
    double energy = 0.0, dis, disquare;
    // atom-atom distance potential
    // get the position of res to check for adjacent residues
    position = res._posn;
    const int List_size = static_cast<int>(List.size());
    for (p = 0; p < List_size; p++) {
        i = List[p];
        if (conf._res[i]._center.x == 0) {
            if (i != position)
                continue;
        }
        if (position != End) {
            if (i == position) {
                if (type) {
                    if (res._atom[1]._type == UNDEF || conf._res[i]._atom[0]._type == UNDEF ||
                        res._atom[1].x == 0.000 || conf._res[i]._atom[0].x == 0.000)
                        continue;
                    else {
                        disquare = (res._atom[1].x - conf._res[i]._atom[0].x) *
                                   (res._atom[1].x - conf._res[i]._atom[0].x) +
                                   (res._atom[1].y - conf._res[i]._atom[0].y) *
                                   (res._atom[1].y - conf._res[i]._atom[0].y) +
                                   (res._atom[1].z - conf._res[i]._atom[0].z) *
                                   (res._atom[1].z - conf._res[i]._atom[0].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            disInd = (int) (dis / H_INLO);
                            energy += PF::LOODIS[res._atom[1]._type - 1][conf._res[i]._atom[0]._type -
                                                                         1][disInd];
                        }
                    }
                }
            }
            else if (i > position && i < End) continue;
            else {
                if (res._bbc.dis(conf._res[i]._center) <
                    Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE) {
                    for (j = 0; j < NUM_BB_ATOM; j++) {
                        // skip undefined and H atoms
                        if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                            continue;
                        if (res._atom[j].x == 0.000) continue;
                        for (k = 0; k < conf._res[i]._numAtom; ++k) {

                            if (i == End) if (type) if (k != ATM_C && k != ATM_CA)
                                continue;
                            if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                            if (conf._res[i]._atom[k].x == 0.000) continue;
                            if (position - i == 1) {
                                if (k == ATM_C && (j == ATM_C || j == ATM_CB))
                                    continue;
                                if (type) {
                                    if (position > Start) if ((k == ATM_C || k == ATM_O) &&
                                                              (j == ATM_C || j == ATM_O || j == ATM_CB))
                                        continue;
                                }
                                else {
                                    if ((k <= 5 && j == ATM_N) ||
                                        ((k == ATM_CA || k == ATM_C || k == ATM_O) && j == ATM_CA))
                                        continue;
                                }
                            }
                            disquare = (res._atom[j].x - conf._res[i]._atom[k].x) *
                                       (res._atom[j].x - conf._res[i]._atom[k].x) +
                                       (res._atom[j].y - conf._res[i]._atom[k].y) *
                                       (res._atom[j].y - conf._res[i]._atom[k].y) +
                                       (res._atom[j].z - conf._res[i]._atom[k].z) *
                                       (res._atom[j].z - conf._res[i]._atom[k].z);
                            if (disquare <= PF_DIS_CUT_SQUARE) {
                                dis = sqrt(disquare);
                                disInd = (int) (dis / H_INLO);
                                energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type -
                                                                             1][disInd];
                            }
                        }
                    }
                }
            }
        }
        else {
            if (i <= position && i >= position - 2) continue;
            if (res._bbc.dis(conf._res[i]._center) < Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE) {
                for (j = 0; j < NUM_BB_ATOM; j++) {
                    // skip undefined and H atoms
                    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                        continue;
                    if (i >= Start && i < End)    //CA and C have been calculated previously
                    {
                        if (j == 1 || j == 2)
                            continue;
                    }
                    if (res._atom[j].x == 0.000) continue;
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0.000) continue;
                        if (position - i == 1) {
                            if ((k == ATM_C && j == ATM_CB) || (k <= 5 && j == ATM_N))
                                continue;
                        }
                        if (i - position == 1) {
                            if ((k == ATM_N && j <= 5) || (k == ATM_CA && (j == ATM_CA || j == ATM_C || j == ATM_O)) ||
                                (j == ATM_C && (k == ATM_C || k == ATM_CB)))
                                continue;
                        }
                        disquare = (res._atom[j].x - conf._res[i]._atom[k].x) *
                                   (res._atom[j].x - conf._res[i]._atom[k].x) +
                                   (res._atom[j].y - conf._res[i]._atom[k].y) *
                                   (res._atom[j].y - conf._res[i]._atom[k].y) +
                                   (res._atom[j].z - conf._res[i]._atom[k].z) *
                                   (res._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            disInd = (int) (dis / H_INLO);
                            energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type -
                                                                         1][disInd];

                        }
                    }
                }
            }
        }
    }

    return energy;
}

// side chain interaction energy for side chain modeling.
// only side chain atoms are considered
// Note: side chain center for the residue _scc should be calculated using Residue::cal_scc()
// before calling this function
double one_res_en_loodis_sc(const Structure &conf, const Residue &res, const int position, const int Start,
                            const int End) {
    int i, j, k, disInd;
    double dis, disquare;
    if (res._type == ALA || res._type == GLY) {
        return 0.0;
    }

    double energy = 0.0;

    // atom-atom distance potential
    for (i = 1; i < position; ++i) {
        if ((res._scc.dis(conf._res[i]._center) <
             (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) || i >= Start) {
            // Here we can use another distance cutoff instead of CC_DIS_CUT which is 12 for whole residues
            // the distance cutoff here for sidechains can be considerablly smaller and residue dependent
            for (k = 0; k < conf._res[i]._numAtom; ++k) {
                if (VdwUtils::should_ignore_atom(conf._res[i]._atom[k])) { continue; }

                for (j = NUM_BB_ATOM; j < res._numAtom; ++j) {
                    if (VdwUtils::should_ignore_atom(res._atom[j])) { continue; }
                    disquare = res._atom[j].disquare(conf._res[i]._atom[k]);

                    if (disquare <= PF_DIS_CUT_SQUARE) {
                        dis = sqrt(disquare);
                        disInd = (int) (dis / H_INLO);
                        assert(check_loodis_bounds(res._atom[j]._type - 1, conf._res[i]._atom[k]._type - 1, disInd));
                        energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type -
                                                                     1][disInd];
                    }
                }
            }
        }
    }

    for (i = End + 1; i <= conf._numRes; ++i) {
        if (res._scc.dis(conf._res[i]._center) <
            (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
            // Here we can use another distance cutoff instead of CC_DIS_CUT which is 12 for whole residues
            // the distance cutoff here for sidechains can be considerablly smaller and residue dependent
            for (k = 0; k < conf._res[i]._numAtom; ++k) {
                if (VdwUtils::should_ignore_atom(conf._res[i]._atom[k])) { continue; }
                for (j = NUM_BB_ATOM; j < res._numAtom; ++j) {
                    if (VdwUtils::should_ignore_atom(res._atom[j])) { continue; }

                    disquare = res._atom[j].disquare(conf._res[i]._atom[k]);

                    if (disquare <= PF_DIS_CUT_SQUARE) {
                        dis = sqrt(disquare);
                        disInd = (int) (dis / H_INLO);
                        assert(check_loodis_bounds(res._atom[j]._type - 1, conf._res[i]._atom[k]._type - 1, disInd));
                        energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type -
                                                                     1][disInd];
                    }
                }
            }
        }
    }

    // C atom in End
    if (position < End) {
        for (i = position + 1; i <= End; ++i) {
            if (res._scc.dis(conf._res[i]._center) <
                (Residue::sc_size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                const int num_bb_atom = std::min(conf._res[i]._numAtom, (short) NUM_BB_ATOM);
                for (k = 0; k < num_bb_atom; ++k) {
                    if (VdwUtils::should_ignore_atom(conf._res[i]._atom[k])) { continue; }
                    for (j = NUM_BB_ATOM; j < res._numAtom; ++j) {
                        if (VdwUtils::should_ignore_atom(res._atom[j])) { continue; }
                        disquare = res._atom[j].disquare(conf._res[i]._atom[k]);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            disInd = (int) (dis / H_INLO);
                            assert(check_loodis_bounds(res._atom[j]._type - 1, conf._res[i]._atom[k]._type - 1,
                                                       disInd));
                            energy += PF::LOODIS[res._atom[j]._type - 1][conf._res[i]._atom[k]._type -
                                                                         1][disInd];
                        }
                    }
                }
            }
        }
    }

    return energy;
}

inline int loodis_get_num_atoms_to_check_against(const Residue &res, const bool is_in_sim_region, const bool use_sc) {
    if (!use_sc && is_in_sim_region) {
        if (res._type != GLY) {
            return NUM_BB_ATOM;
        } else {
            return NUM_BB_ATOM - 1;
        }
    }
    return res._numAtom;
}

// loodis energy calculation, type = false No side-chain, type = true, side chain added
double loodis_e(const Structure &Conf, const vector<int> &starts, const vector<int> &ends, const bool type) {
    int i_loop, ix_res_loop, ix_res_other, j, k, disInd, numatom1, numatom2;
    double dis, r = 0;
    double energy = 0;
    int atomIndex1, atomIndex2;
    const int n_loop = starts.size();

    for (i_loop = 0; i_loop < n_loop; ++i_loop) {
        // Bounds check loop interval
        assert(starts[i_loop] < ends[i_loop]);
        assert(starts[i_loop] >= 1);
        assert(ends[i_loop] < Conf._numRes);
        for (ix_res_loop = starts[i_loop]; ix_res_loop <= ends[i_loop]; ++ix_res_loop) {
            // Check against all other residues
            for (ix_res_other = 1; ix_res_other <= Conf._numRes; ++ix_res_other) {

                // Skip past adjacent and self residues
                const int relative_position = ix_res_loop > ix_res_other ? ix_res_loop - ix_res_other : ix_res_other -
                                                                                                        ix_res_loop;
                assert(relative_position >= 0);
                if (relative_position <= 1) {
                    continue;
                }

                // To avoid double counting: skip past residue indices in simulated regions
                // which are higher than current (outer) index
                const bool is_other_res_in_sim_region = is_in_simulated_region(ix_res_other, starts, ends);
                if ((ix_res_other > ix_res_loop) && is_other_res_in_sim_region) {
                    continue;
                }

                numatom1 = loodis_get_num_atoms_to_check_against(Conf._res[ix_res_loop], true /*is_in_sim_region*/,
                                                                 type > 0);
                numatom2 = loodis_get_num_atoms_to_check_against(Conf._res[ix_res_other], is_other_res_in_sim_region,
                                                                 type > 0);

                if (Conf._res[ix_res_loop]._center.dis(Conf._res[ix_res_other]._center) < CC_DIS_CUT) {
                    for (j = 0; j < numatom1; ++j) {
                        if (VdwUtils::should_ignore_atom(Conf._res[ix_res_loop]._atom[j])) { continue; }
                        for (k = 0; k < numatom2; ++k) {
                            if (VdwUtils::should_ignore_atom(Conf._res[ix_res_other]._atom[k])) { continue; }

                            atomIndex1 = Conf._res[ix_res_loop]._atom[j]._type;
                            atomIndex2 = Conf._res[ix_res_other]._atom[k]._type;
                            dis = Conf._res[ix_res_loop]._atom[j].dis(Conf._res[ix_res_other]._atom[k]);
                            if (dis < (H_INLO * LOODIS_DIS_BIN)) {
                                disInd = (int) (dis / H_INLO);
                                assert(check_loodis_bounds(atomIndex1 - 1, atomIndex2 - 1, disInd));
                                energy += PF::LOODIS[atomIndex1 - 1][atomIndex2 - 1][disInd];
                            }
                        }
                    }
                }
            }
        }
    }
    return energy;
}

// loodis energy calculation, type = 0 No side-chain, type = 1, side chain added, else, loop region itself does not included
double loodis_e_list(const Structure &Conf, const vector<int> &starts, const vector<int> &ends, const vector<int> &List, const int type) {
    int ix_res_loop, ix_res_lst, j, k, i_loop, i_lst, disInd, numatom1, numatom2;
    double dis;
    double energy = 0;
    int atomIndex1, atomIndex2;
    const int n_loop = starts.size();

    const int List_size = static_cast<int>(List.size());

    for (i_loop = 0; i_loop < n_loop; ++i_loop) {
        // Bounds check loop interval
        assert(starts[i_loop] < ends[i_loop]);
        assert(starts[i_loop] >= 1);
        assert(ends[i_loop] < Conf._numRes);
        for (ix_res_loop = starts[i_loop]; ix_res_loop <= ends[i_loop]; ++ix_res_loop) {
            // Some of the inner loop continues only make sense if ix_res_loop is within List
            // as we are assuming that some higher indexed comparison will account for the
            // energy (must be true for all ix_res_loop)
            assert(std::find(List.cbegin(), List.cend(), ix_res_loop) != List.cend());
            for (i_lst = 0; i_lst < List_size; i_lst++) {
                ix_res_lst = List[i_lst];
                // Bounds check list index
                assert(ix_res_lst >= 1 && ix_res_lst <= Conf._numRes);

                // Skip past adjacent and self residues
                const int relative_position = ix_res_loop > ix_res_lst ? ix_res_loop - ix_res_lst : ix_res_lst -
                                                                                                    ix_res_loop;
                assert(relative_position >= 0);
                if (relative_position <= 1) {
                    continue;
                }

                // To avoid double counting: skip past residue indices in simulated regions
                // which are higher than current (outer) index
                const bool is_lst_res_in_sim_region = is_in_simulated_region(ix_res_lst, starts, ends);
                if ((ix_res_lst > ix_res_loop) && is_lst_res_in_sim_region) {
                    continue;
                }

                numatom1 = loodis_get_num_atoms_to_check_against(Conf._res[ix_res_loop], true /*is_in_sim_region*/,
                                                                 type > 0);
                numatom2 = loodis_get_num_atoms_to_check_against(Conf._res[ix_res_lst], is_lst_res_in_sim_region,
                                                                 type > 0);

                if (Conf._res[ix_res_loop]._center.dis(Conf._res[ix_res_lst]._center) < CC_DIS_CUT) {
                    for (j = 0; j < numatom1; ++j) {
                        if (VdwUtils::should_ignore_atom(Conf._res[ix_res_loop]._atom[j])) { continue; }
                        for (k = 0; k < numatom2; ++k) {
                            if (VdwUtils::should_ignore_atom(Conf._res[ix_res_lst]._atom[k])) { continue; }
                            atomIndex1 = Conf._res[ix_res_loop]._atom[j]._type;
                            atomIndex2 = Conf._res[ix_res_lst]._atom[k]._type;
                            dis = Conf._res[ix_res_loop]._atom[j].dis(Conf._res[ix_res_lst]._atom[k]);
                            if (dis < (H_INLO * LOODIS_DIS_BIN)) {
                                disInd = (int) (dis / H_INLO);
                                assert(check_loodis_bounds(atomIndex1 - 1, atomIndex2 - 1, disInd));
                                if (Conf._res[ix_res_loop]._type == 1 && Conf._res[ix_res_lst]._type == 1) {
                                    if ((Conf._res[ix_res_loop]._atom[j]._type == 16) &&
                                        (Conf._res[ix_res_lst]._atom[k]._type == 16) &&
                                        (dis > 1.8 && dis < 2.2)) {
                                        energy += -2.5;
                                    }
                                    else if (((Conf._res[ix_res_loop]._atom[j]._type == 16 &&
                                               Conf._res[ix_res_lst]._atom[k]._type == 4) ||
                                              (Conf._res[ix_res_loop]._atom[j]._type == 4 &&
                                               Conf._res[ix_res_lst]._atom[k]._type == 16)) &&
                                             (dis > 2.7 && dis < 3.3)) {
                                        energy += -2;
                                    }
                                    else
                                        energy += PF::LOODIS[atomIndex1 - 1][atomIndex2 - 1][disInd];
                                }
                                else
                                    energy += PF::LOODIS[atomIndex1 - 1][atomIndex2 - 1][disInd];
                            }
                        }
                    }
                }
            }
        }
    }
    return energy;
}

// Steric Clash detection for single loop + REDCELL
void Clash_detection_list(const Structure &conf, const int Start, const int End, vector<int> &ResIdx,
                          vector<int> &ClashNum, const std::vector<int> &List) {

    int i, p, j, l, k, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    for (l = Start; l <= End; l++) {
        if (conf._res[l]._type == 0 || conf._res[l]._type == 5)
            continue;
        for (p = 0; p < List_size; p++) {
            i = List[p];

            if (i == l) {
                for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                    // skip undefined and H atoms
                    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                        continue;
                    if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                        conf._res[l]._atom[j].z == 0)
                        continue;
                    for (k = 2; k <= 3; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                   (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                   + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                     (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                   + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                     (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }
            }
            else if (i > l && i <= End) continue;
            else if (i >= Start && i < l) {
                for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                    // skip undefined and H atoms
                    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                        continue;
                    if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                        conf._res[l]._atom[j].z == 0)
                        continue;
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                   (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                   + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                     (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                   + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                     (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }

            }
            else {
                if ((conf._res[l]._center.dis(conf._res[i]._center) <
                     Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                    (conf._res[l]._bbc.dis(conf._res[i]._center) <
                     Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                    for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                        // skip undefined and H atoms
                        if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                            continue;
                        if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                            conf._res[l]._atom[j].z == 0)
                            continue;
                        for (k = 0; k < conf._res[i]._numAtom; ++k) {
                            if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                            if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                conf._res[i]._atom[k].z == 0)
                                continue;
                            disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                       (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                       + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                         (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                       + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                         (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                            if (disquare <= PF_DIS_CUT_SQUARE) {
                                dis = sqrt(disquare);
                                r = Atom::radius[conf._res[l]._atom[j]._type] +
                                    Atom::radius[conf._res[i]._atom[k]._type];
                                quot = dis / r;
                                if (quot <= VDW_CLASH_CUTOFF) {
                                    if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                        clash_count++;
                                    else
                                        clash_count += 5;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (clash_count != 0) {
            ResIdx.push_back(l);
            ClashNum.push_back(clash_count);
        }
        clash_count = 0;
    }
}

// Steric Clash Detection for one residue
int Res_clash_detection_list(const Structure &conf, const Residue &res,
                             const int Start, const int End, 
                             const std::vector<int> &List) {
    int p, i, j, k, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    for (p = 0; p < List_size; p++) {
        i = List[p];
        if (i == res._posn) {
            for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                    continue;
                if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                for (k = 2; k <= 3; ++k) {
                    if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                    if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                        conf._res[i]._atom[k].z == 0)
                        continue;
                    disquare = (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                               + (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                               +
                               (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                    if (disquare <= PF_DIS_CUT_SQUARE) {
                        dis = sqrt(disquare);
                        r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                        quot = dis / r;
                        if (quot <= VDW_CLASH_CUTOFF) {
                            if (quot <= VDW_CLASH_CUTOFF && quot > 0.5)
                                clash_count++;
                            else
                                clash_count += 5;
                        }
                    }
                }
            }
        }
        else if (i > res._posn && i <= End) continue;
        else if (i >= Start && i < res._posn) {
            for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                    continue;
                if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                for (k = 0; k < conf._res[i]._numAtom; ++k) {
                    if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                    if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                        conf._res[i]._atom[k].z == 0)
                        continue;
                    disquare = (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                               + (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                               +
                               (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                    if (disquare <= PF_DIS_CUT_SQUARE) {
                        dis = sqrt(disquare);
                        r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                        quot = dis / r;
                        if (quot <= VDW_CLASH_CUTOFF) {
                            if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                clash_count++;
                            else
                                clash_count += 5;
                        }
                    }
                }
            }
        }
        else {
            if ((res._center.dis(conf._res[i]._center) <
                 Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                (res._bbc.dis(conf._res[i]._center) <
                 Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                    // skip undefined and H atoms
                    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                        continue;
                    if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                    //if(type != 0 && j >= NUM_BB_ATOM ) continue;  // when sampling only backbone atoms, skip side chain atoms if there are any
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        disquare =
                                (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                                +
                                (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                                +
                                (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }
            }
        }
    }
    return clash_count;
}

void Clash_detection_list(const Structure &conf, const vector<int> &starts, const vector<int> &ends, vector<int> &ResIdx,
                          vector<int> &ClashNum, const std::vector<int> &List) {
    int i, p, j, l, k, t, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    const int n_starts = static_cast<int>(starts.size());

    for (t = 0; t < n_starts; t++) {
        for (l = starts[t]; l <= ends[t]; l++) {
            if (conf._res[l]._type == 0 || conf._res[l]._type == 5)
                continue;
            for (p = 0; p < List_size; p++) {
                i = List[p];

                if (i == l) {
                    for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                        // skip undefined and H atoms
                        if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                            continue;
                        if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                            conf._res[l]._atom[j].z == 0)
                            continue;

                        for (k = 2; k <= 3; ++k) {
                            if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                            if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                conf._res[i]._atom[k].z == 0)
                                continue;

                            disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                       (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                       + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                         (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                       + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                         (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                            if (disquare <= PF_DIS_CUT_SQUARE) {
                                dis = sqrt(disquare);
                                r = Atom::radius[conf._res[l]._atom[j]._type] +
                                    Atom::radius[conf._res[i]._atom[k]._type];
                                if (conf._res[l]._type == 1 && conf._res[i]._type == 1) {
                                    if ((conf._res[l]._atom[j]._type == 16) && (conf._res[i]._atom[k]._type == 16) &&
                                        (dis > 1.8 && dis < 2.2)) {
                                        continue;
                                    }
                                    else if (((conf._res[l]._atom[j]._type == 16 && conf._res[i]._atom[k]._type == 4) ||
                                              (conf._res[l]._atom[j]._type == 4 &&
                                               conf._res[i]._atom[k]._type == 16)) && (dis > 2.7 && dis < 3.3)) {
                                        continue;
                                    }
                                    else {
                                        quot = dis / r;
                                        if (quot <= VDW_CLASH_CUTOFF) {
                                            if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                                clash_count++;
                                            else {
                                                clash_count += 5;
                                            }
                                        }
                                    }
                                }
                                else {
                                    quot = dis / r;
                                    if (quot <= VDW_CLASH_CUTOFF) {
                                        if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                            clash_count++;
                                        else {
                                            clash_count += 5;
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
                else if (i > l && i <= ends[t]) continue;
                else if (i >= starts[t] && i < l) {
                    for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                        // skip undefined and H atoms
                        if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                            continue;
                        if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                            conf._res[l]._atom[j].z == 0)
                            continue;

                        for (k = 0; k < conf._res[i]._numAtom; ++k) {
                            if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                            if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                conf._res[i]._atom[k].z == 0)
                                continue;

                            disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                       (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                       + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                         (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                       + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                         (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                            if (disquare <= PF_DIS_CUT_SQUARE) {
                                dis = sqrt(disquare);
                                r = Atom::radius[conf._res[l]._atom[j]._type] +
                                    Atom::radius[conf._res[i]._atom[k]._type];
                                if (conf._res[l]._type == 1 && conf._res[i]._type == 1) {
                                    if ((conf._res[l]._atom[j]._type == 16) && (conf._res[i]._atom[k]._type == 16) &&
                                        (dis > 1.8 && dis < 2.2)) {
                                        continue;
                                    }
                                    else if (((conf._res[l]._atom[j]._type == 16 && conf._res[i]._atom[k]._type == 4) ||
                                              (conf._res[l]._atom[j]._type == 4 &&
                                               conf._res[i]._atom[k]._type == 16)) && (dis > 2.7 && dis < 3.3)) {
                                        continue;
                                    }
                                    else {
                                        quot = dis / r;
                                        if (quot <= VDW_CLASH_CUTOFF) {
                                            if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                                clash_count++;
                                            else {
                                                clash_count += 5;
                                            }
                                        }
                                    }
                                }
                                else {
                                    quot = dis / r;
                                    if (quot <= VDW_CLASH_CUTOFF) {
                                        if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                            clash_count++;
                                        else {
                                            clash_count += 5;
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
                else {
                    if ((conf._res[l]._center.dis(conf._res[i]._center) <
                         Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                        (conf._res[l]._bbc.dis(conf._res[i]._center) <
                         Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                        for (j = NUM_BB_ATOM; j < conf._res[l]._numAtom; j++) {
                            // skip undefined and H atoms
                            if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                                continue;
                            if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                                conf._res[l]._atom[j].z == 0)
                                continue;
                            for (k = 0; k < conf._res[i]._numAtom; ++k) {
                                if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                                if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                    conf._res[i]._atom[k].z == 0)
                                    continue;

                                disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                           (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                           + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                             (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                           + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                             (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                                if (disquare <= PF_DIS_CUT_SQUARE) {
                                    dis = sqrt(disquare);
                                    r = Atom::radius[conf._res[l]._atom[j]._type] +
                                        Atom::radius[conf._res[i]._atom[k]._type];
                                    if (conf._res[l]._type == 1 && conf._res[i]._type == 1) {
                                        if ((conf._res[l]._atom[j]._type == 16) &&
                                            (conf._res[i]._atom[k]._type == 16) && (dis > 1.8 && dis < 2.2)) {
                                            continue;
                                        }
                                        else if (((conf._res[l]._atom[j]._type == 16 &&
                                                   conf._res[i]._atom[k]._type == 4) ||
                                                  (conf._res[l]._atom[j]._type == 4 &&
                                                   conf._res[i]._atom[k]._type == 16)) && (dis > 2.7 && dis < 3.3)) {
                                            continue;
                                        }
                                        else {
                                            quot = dis / r;
                                            if (quot <= VDW_CLASH_CUTOFF) {
                                                if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                                    clash_count++;
                                                else {
                                                    clash_count += 5;
                                                }
                                            }
                                        }
                                    }
                                    else {
                                        quot = dis / r;
                                        if (quot <= VDW_CLASH_CUTOFF) {
                                            if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                                clash_count++;
                                            else {
                                                clash_count += 5;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (clash_count != 0) {
                ResIdx.push_back(l);
                ClashNum.push_back(clash_count);
            }
            clash_count = 0;
        }
    }
}

// Backbone Loop
void BBClash_detection_list(const Structure &conf, const int Start, const int End,
                            vector<int> &ResIdx, vector<int> &ClashNum,
                            const std::vector<int> &List) {

    int i, p, j, l, k, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    for (l = Start; l <= End; l++) {
        for (p = 0; p < List_size; p++) {
            i = List[p];
            if (i >= l && i <= End) continue;
            if ((conf._res[l]._center.dis(conf._res[i]._center) <
                 Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                (conf._res[l]._bbc.dis(conf._res[i]._center) <
                 Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                for (j = 0; j < NUM_BB_ATOM; j++) {
                    // skip undefined and H atoms
                    if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                        continue;
                    if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                        conf._res[l]._atom[j].z == 0)
                        continue;
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        if (i == l + 1) {
                            if (j == ATM_CA && k == ATM_CA)
                                continue;
                            if (j == ATM_C && (k == ATM_N || k == ATM_CA))
                                continue;
                            if (j == ATM_O && k == ATM_N)
                                continue;
                        } else if (i == l - 1) {
                            if (k == ATM_CA && j == ATM_CA)
                                continue;
                            if (k == ATM_C && (j == ATM_N || j == ATM_CA))
                                continue;
                            if (k == ATM_O && j == ATM_N) continue;
                        }
                        disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                   (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                   + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                     (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                   + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                     (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[conf._res[l]._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }
            }
        }
        if (clash_count != 0) {
            ResIdx.push_back(l);
            ClashNum.push_back(clash_count);
        }
        clash_count = 0;
    }
}


void BBClash_detection_list(const Structure &conf, const vector<int> &starts, const vector<int> &ends,
                            vector<int> &ResIdx, vector<int> &ClashNum, const std::vector<int> &List) {
    int i, p, j, l, k, t, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    const int n_starts = static_cast<int>(starts.size());

    for (t = 0; t < n_starts; ++t) {
        for (l = starts[t]; l <= ends[t]; ++l) {
            for (p = 0; p < List_size; ++p) {
                i = List[p];

                //if( i >= l) continue;
                if (i >= l && i <= ends[t]) continue;
                else if (i >= starts[t] && i < l) {
                    for (j = 0; j < NUM_BB_ATOM; ++j) {
                        // skip undefined and H atoms
                        if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                            continue;
                        if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                            conf._res[l]._atom[j].z == 0)
                            continue;
                        for (k = 0; k < NUM_BB_ATOM; ++k) {
                            if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                            if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                conf._res[i]._atom[k].z == 0)
                                continue;
                            if (i == l + 1) {
                                if (j == ATM_CA && k == ATM_CA)
                                    continue;
                                if (j == ATM_C && (k == ATM_N || k == ATM_CA))
                                    continue;
                                if (j == ATM_O && k == ATM_N)
                                    continue;
                            } else if (i == l - 1) {
                                if (k == ATM_CA && j == ATM_CA)
                                    continue;
                                if (k == ATM_C && (j == ATM_N || j == ATM_CA))
                                    continue;
                                if (k == ATM_O && j == ATM_N)
                                    continue;
                            }
                            disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                       (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                       + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                         (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                       + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                         (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                            if (disquare <= PF_DIS_CUT_SQUARE) {
                                dis = sqrt(disquare);
                                r = Atom::radius[conf._res[l]._atom[j]._type] +
                                    Atom::radius[conf._res[i]._atom[k]._type];
                                quot = dis / r;
                                if (quot <= VDW_CLASH_CUTOFF) {
                                    if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5) {
                                        clash_count++;
                                    }
                                    else {
                                        clash_count += 5;
                                    }
                                    //  cout <<"Severe Steric Clash same residue case 1!!	   "<<l<<"   "<<conf._res[l]._type<<"   "<<j<<"    "<<i<<"   "<<conf._res[i]._type<<"   "<<k<<"   "<<quot<<endl;
                                }
                            }
                        }
                    }
                }
                else {
                    if ((conf._res[l]._center.dis(conf._res[i]._center) <
                         Residue::size[conf._res[l]._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                        (conf._res[l]._bbc.dis(conf._res[i]._center) <
                         Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                        for (j = 0; j < NUM_BB_ATOM; j++) {
                            // skip undefined and H atoms
                            if (conf._res[l]._atom[j]._type == UNDEF || conf._res[l]._atom[j]._type >= 22)
                                continue;
                            if (conf._res[l]._atom[j].x == 0 && conf._res[l]._atom[j].y == 0 &&
                                conf._res[l]._atom[j].z == 0)
                                continue;
                            for (k = 0; k < conf._res[i]._numAtom; ++k) {
                                if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                                if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                                    conf._res[i]._atom[k].z == 0)
                                    continue;
                                if (i == l + 1) {
                                    if (j == ATM_CA && k == ATM_CA)
                                        continue;
                                    if (j == ATM_C && (k == ATM_N || k == ATM_CA))
                                        continue;
                                    if (j == ATM_O && k == ATM_N)
                                        continue;
                                } else if (i == l - 1) {
                                    if (k == ATM_CA && j == ATM_CA)
                                        continue;
                                    if (k == ATM_C && (j == ATM_N || j == ATM_CA))
                                        continue;
                                    if (k == ATM_O && j == ATM_N) continue;
                                }

                                disquare = (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x) *
                                           (conf._res[l]._atom[j].x - conf._res[i]._atom[k].x)
                                           + (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y) *
                                             (conf._res[l]._atom[j].y - conf._res[i]._atom[k].y)
                                           + (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z) *
                                             (conf._res[l]._atom[j].z - conf._res[i]._atom[k].z);
                                if (disquare <= PF_DIS_CUT_SQUARE) {
                                    dis = sqrt(disquare);
                                    r = Atom::radius[conf._res[l]._atom[j]._type] +
                                        Atom::radius[conf._res[i]._atom[k]._type];
                                    quot = dis / r;
                                    if (quot <= VDW_CLASH_CUTOFF) {
                                        if (quot <= VDW_CLASH_CUTOFF && quot >= 0.5)
                                            clash_count++;
                                        else
                                            clash_count += 5;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (clash_count != 0) {
                ResIdx.push_back(l);
                ClashNum.push_back(clash_count);
            }
            clash_count = 0;
        }
    }
}

// Performs one-to-all clash detection. Currently, this is called within the side chain clash
// minimization routine (SCE_minimization*)
// @param conf - the parent protein conformation
// @param res - the residue to perform clash detection against
// @param Start - vector containing start residue numbers for loops being grown
// @param End - vector containing end residue numbers for loops being grown
// @param loopidx - the index of the loop containing res. Start[loopidx] and End[loopidx]
//  give start and end residue numbers for that loop respectively.
// @param List - The set of residue numbers within a bounding ellipsoid. Only residues within
//  the ellipsoid are checked against for clashes
int Res_clash_detection_list(const Structure &conf, const Residue &res, const vector<int> &starts,
                             const vector<int> &ends, const int loopidx, const std::vector<int> &List) {
    int p, i, j, k, t, clash_count = 0;
    double disquare, dis, r, quot = 0;

    const int List_size = static_cast<int>(List.size());

    // Loop - check for clashes with residues/atoms within ellipsoid
    for (p = 0; p < List_size; ++p) {
        i = List[p];
        // Case - residue to be checked against is within same loop region as param Res
        if (i >= starts[loopidx] && i <= ends[loopidx]) {
            // Case - residue to be checked against is within loop region and is
            // the residue of interest (this is a self-clash check)
            if (i == res._posn) {
                // For the self check - we see if the residue side chain atoms collide
                // against backbone atoms at indices 2 & 3 which is only the carbonyl carbon
                // and oxygen. Also, because we are rotating around CB-CA axis, given that
                // they were not previously clashing with the side chain, then the rotation
                // about that axis should not induce a clash with CB or CA (hence, why they
                // are not checked).

                // @TODO - verify that CA, CB does not need to be checked
                // @TODO - verify that peptide N does not need to be checked
                // @TODO - investigate if all the steric clash detection in general
                // should be replaced instead with some sort of Van der Waals/Lennard-Jones potential
                // as, per discussion with Ke, the clash cost/count is very "ad-hoc"

                // Loop over side chain atoms for param Res
                for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                    // Skip if atom type is undefined or hydrogen
                    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                        continue;
                    // Skip if atom center is at origin (implies uninitialized?)
                    if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                    // Only check against 2 = peptide carbonyl carbon and 3 = peptide carbonyl oxygen
                    for (k = 2; k <= 3; ++k) {
                        // Skip if for some reason atom type is undefined or is hydrogen
                        // @TODO - probably should be an assert
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        // Compute distance between side chain atom center and backbone atom center
                        disquare =
                                (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                                +
                                (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                                +
                                (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                        // Check if squared distance is less than cut off threshold
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            // Compute ratio of distance / (sum of vdw radii)
                            dis = sqrt(disquare);
                            r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot > 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }

            }
                // Case - we are checking residue side chain atoms against another loop residue
                // in same loop as param Res, but this residue is not the same as param Res
            else {
                // @TODO - what about case i > posn and i < end for the current loop?
                //  - These side chains have not undergone SCE_minimization. In single loop
                //  code, they are explicitly ignored.

                // Loop over side chain atoms for param Res
                for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                    // Skip if atom is undefined or is hydrogen
                    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                        continue;
                    // Skip if atom is at origin => uninitialized?
                    if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                    // Check clashes against all atoms
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        disquare =
                                (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                                +
                                (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                                +
                                (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }
            }
        }
            // Case residue to be checked against is not in loop region of param Res
        else {
            // Broad phase check - see if bounding spheres for side chain and backbone are close enough
            // for a narrow phase check
            if ((res._center.dis(conf._res[i]._center) <
                 Residue::size[res._type] + Residue::size[conf._res[i]._type] + CUB_SIZE) ||
                (res._bbc.dis(conf._res[i]._center) <
                 Residue::bb_size + Residue::size[conf._res[i]._type] + CUB_SIZE)) {
                // Avoid clash detection with loops that have yet to be minimized
                int Continue_label = 0;
                for (t = loopidx + 1; t < starts.size(); t++) {
                    if (i >= starts[t] && i <= ends[t]) {
                        Continue_label = 1;
                        break;
                    }
                }
                if (Continue_label == 1) continue;


                // Check against param Res side chain atoms vs all atoms of current residue
                for (j = NUM_BB_ATOM; j < res._numAtom; j++) {
                    // skip undefined and H atoms
                    if (res._atom[j]._type == UNDEF || res._atom[j]._type >= 22)
                        continue;
                    // skip atoms located at origin => uninitialized?
                    if (res._atom[j].x == 0 && res._atom[j].y == 0 && res._atom[j].z == 0) continue;
                    // Check against all atoms (backbone + side chain)
                    for (k = 0; k < conf._res[i]._numAtom; ++k) {
                        if (conf._res[i]._atom[k]._type == UNDEF || conf._res[i]._atom[k]._type >= 22) continue;
                        if (conf._res[i]._atom[k].x == 0 && conf._res[i]._atom[k].y == 0 &&
                            conf._res[i]._atom[k].z == 0)
                            continue;
                        disquare =
                                (res._atom[j].x - conf._res[i]._atom[k].x) * (res._atom[j].x - conf._res[i]._atom[k].x)
                                +
                                (res._atom[j].y - conf._res[i]._atom[k].y) * (res._atom[j].y - conf._res[i]._atom[k].y)
                                +
                                (res._atom[j].z - conf._res[i]._atom[k].z) * (res._atom[j].z - conf._res[i]._atom[k].z);
                        if (disquare <= PF_DIS_CUT_SQUARE) {
                            dis = sqrt(disquare);
                            r = Atom::radius[res._atom[j]._type] + Atom::radius[conf._res[i]._atom[k]._type];
                            quot = dis / r;
                            if (quot <= VDW_CLASH_CUTOFF) {
                                if (quot >= 0.5)
                                    clash_count++;
                                else
                                    clash_count += 5;
                            }
                        }
                    }
                }
            }
        }
    }

    return clash_count;
}
