//****************************************************
// grow_one_util.cpp
//
// Implementation of incremental growth function
//  which grows both backbone and side chain
//  atoms for a single residue.
//****************************************************

#include "build.h"
#include "smc.h"
#include "vla.h"
#include "vdw_utils.h"

#include <assert.h>
#include <limits>
#include <iostream>

namespace {
    // Hack for rand utilities returning values outside target range producing NaNs
    // e.g. FRAND(0,1) due to 0.5 padding can actually return a value > 1. Some utilities
    // expect this value to be strictly less than or equal to 1; however, boxMullerSample*
    // utility will end up taking a negative log if value is greater than 1 which produces
    // a NaN result. The proper fix is to fix the sampling library and replace it with
    // C++11 <random>. However, quick fix for now is to just trap this case and move on.
    // @TODO - fix this for real - use C++11 <random>
    inline bool conditional_mark_dead_label(const double val, int &torProb, int &deadLabel, int &total) {
        if (!std::isfinite(val)) {
            torProb = 0;
            deadLabel = 1;
            ++total;
            return true;
        }
        return false;
    }

    /**
     * Enumerated torsion angles
     */
    enum eTorsion {
        eix_phi = 0,
        eix_psi,
        eix_omega,
        eix_num_torsion
    };

} // end of helper namespace

// loodis energy calculation for a single residue
// @TODO - Verify N atom of next residue to be grown is properly included
double one_res_en_loodis_all_list(const Structure &Conf, const Residue &res, const std::vector<int> &List) {
    int i, j, k, p, disInd, numatom2;
    double dis;
    double energy = 0;
    int atomIndex1, atomIndex2;

    // get the position of res to check for adjacent residues
    // @TODO - verify this is initialized properly
    const int position = res._posn;

    const int numatom1 = res._numAtom;

    const int List_size = static_cast<int>(List.size());

    for (p = 0; p < List_size; p++) {
        i = List[p];

        // exclude neighbor residue
        if (i < position + 2 && i > position - 2)
            continue;

        // exclude residues with x center at 0 -> uninitialized?
        // @TODO - verify this
        if (Conf._res[i]._center.x == 0) {
            if (i != position)
                continue;
        }

        // skip past uninitialized residues and atom positions
        if (res._atom[1]._type == UNDEF || Conf._res[i]._atom[0]._type == UNDEF || res._atom[1].x == 0. ||
            Conf._res[i]._atom[0].x == 0.)
            continue;

        numatom2 = Conf._res[i]._numAtom;

        if (res._center.dis(Conf._res[i]._center) < CC_DIS_CUT) {
            for (j = 0; j < numatom1; ++j) {
                if (VdwUtils::should_ignore_atom(res._atom[j])) { continue; }
                for (k = 0; k < numatom2; ++k) {
                    if (VdwUtils::should_ignore_atom(Conf._res[i]._atom[k])) { continue; }
                    atomIndex1 = res._atom[j]._type;
                    atomIndex2 = Conf._res[i]._atom[k]._type;
                    dis = res._atom[j].dis(Conf._res[i]._atom[k]);
                    if (dis < (H_INLO * LOODIS_DIS_BIN)) {
                        disInd = (int) (dis / H_INLO);
                        if (res._type == 1 && Conf._res[i]._type == 1) {
                            if ((res._atom[j]._type == 16) && (Conf._res[i]._atom[k]._type == 16) &&
                                (dis > 1.8 && dis < 2.2)) {
                                energy += -2.5;
                            }
                            else if (((res._atom[j]._type == 16 && Conf._res[i]._atom[k]._type == 4) ||
                                      (res._atom[j]._type == 4 && Conf._res[i]._atom[k]._type == 16)) &&
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

    return energy;
}

// This version only grows backbone atoms (no side chains)
// A number of candidate conformations for a residue are sampled
// one of them is chosen and its weight is updated.
// This incremental growth method grows atoms in the following fashion:
// CB -> Residue at param Position
// C  -> Residue at param Position
// O  -> Residue at param Position
// N  -> Residue at param Position + 1
// CA -> Residue at param Position + 1
bool SMC::grow_one_bb_only(Structure &CurConf, const int Position, const int tmpEnd, const int LoopChosen,
    const Atom &EndPt) {
    // Initialize positional values
    const int RemLength = tmpEnd - Position;
    const int LargeNumDistanceStates = 160;
    // phi_int, psi_int are the bin indices containing the observed counts
    // for (phi, psi) angle pair. 
    int phi_int, psi_int;

    // Initialize numbers of states
    const int NumStates = this->num_dist_states;

    // The current residue
    const int CurResType = CurConf._res[Position]._type;
    // bb_angles are stored in radians
    DECLARE_VLA_2D(double, bb_angles, NumStates, eix_num_torsion);
    // larbb_angles are stored in
    //  -phi, psi: degrees, omega: radians
    DECLARE_VLA_2D(double, larbb_angles, LargeNumDistanceStates, eix_num_torsion);

    DECLARE_VLA_1D(double, energy, NumStates);
    DECLARE_VLA_1D(int, Phi, LargeNumDistanceStates);
    DECLARE_VLA_1D(int, Psi, LargeNumDistanceStates);
    DECLARE_VLA_1D(int, TorProb, LargeNumDistanceStates);
    DECLARE_VLA_1D(int, deadlabel, LargeNumDistanceStates);
    const int binNum = 180 / BBTbinSize;
    std::vector<Residue> tmpRes(NumStates);        // position n
    vector<double> prob(NumStates, 0.0);
    const double EEdisC = CurConf._res[Position]._atom[ATM_CA].dis(EndPt);
    double EEdisN = 0;
    const double minEEDisC = minDistcon[1][RemLength];
    const double maxEEDisC = minDistcon[1][RemLength] + 31 * DistconBy[1][RemLength];
    const double minEEDisN = minDistcon[0][RemLength - 1];
    const double maxEEDisN = minDistcon[0][RemLength - 1] + 31 * DistconBy[0][RemLength - 1];
    int num_dead_labels = 0;

    SET_VLA_1D(int, deadlabel, LargeNumDistanceStates, 0);

    // Sample end-to-end distances and convert to positions
    {

        // Early out if distance constraints are not satisfied
        // for carbon atom
        if (EEdisC < minEEDisC || EEdisC > maxEEDisC) {
            return false;
        }

        // Loop through states
        for (int i = 0; i < LargeNumDistanceStates; ++i) {
            // Initialize four residues for this sample
            Residue sam;
            // Sample new points for C atom based on distance to final C
            TorProb[i] = -1;

            // Sample carbon atom
            sample_distance(CurConf._res[Position]._atom[ATM_N],
                CurConf._res[Position]._atom[ATM_CA],
                EndPt,
                Residue::bond_angle[CurResType][ATM_C],
                Residue::bond_length[CurResType][ATM_C],
                sam._atom[ATM_C], EEdisC,
                1, RemLength);

            // Skip if sampling failed
            const double atm_c_x = sam._atom[ATM_C].x;
            if (conditional_mark_dead_label(atm_c_x, TorProb[i], deadlabel[i], num_dead_labels)) {
                continue;
            }

            // Sample new points for N atom based on distance to final C
            EEdisN = sam._atom[ATM_C].dis(EndPt);
            if (EEdisN >= minEEDisN && EEdisN <= maxEEDisN) {
                sample_distance(CurConf._res[Position]._atom[ATM_CA],
                    sam._atom[ATM_C],
                    EndPt,
                    Residue::bond_angle[CurResType][ATM_N],
                    Residue::bond_length[CurResType][ATM_N],
                    sam._atom[ATM_N], EEdisN,
                    // -1 since it's NEXT N
                    0, RemLength - 1);
            }
            else {
                TorProb[i] = 0;
                deadlabel[i] = 1;
                ++num_dead_labels;
                continue;
            }

            if (conditional_mark_dead_label(sam._atom[ATM_N].x, TorProb[i], deadlabel[i], num_dead_labels)) {
                continue;
            }

            // Populate remainder of bb_angles with the torsion angles
            // implied by the positional selections for C and N. This is for
            // unity with the existing code base. Omega = pi by default.
            // Note: phi, psi are in degrees (because rotamer selection bins
            // are in degrees) but omega is in radians.
            larbb_angles[i][eix_phi] =
                torsion(CurConf._res[Position - 1]._atom[ATM_C],
                CurConf._res[Position]._atom[ATM_N],
                CurConf._res[Position]._atom[ATM_CA],
                sam._atom[ATM_C],
                torsion_degrees);
            larbb_angles[i][eix_psi] =
                torsion(CurConf._res[Position]._atom[ATM_N],
                CurConf._res[Position]._atom[ATM_CA],
                sam._atom[ATM_C],
                sam._atom[ATM_N],
                torsion_degrees);
            // Sample omega (peptide bond) angle from normal distribution
            box_MullerNsample_single(larbb_angles[i][eix_omega], PI, TO_RAD(4));

            // @TODO - this occurs because FRAND(0,1) can return a value > 1 which
            // causes box_Muller* to take a log of a negative value -> NaN
            // The proper fix is to replace all calls to rand() with C++11 calls
            // from <random> library. But for now, we'll just handle this rare edge case.
            if (conditional_mark_dead_label(larbb_angles[i][eix_omega], TorProb[i], deadlabel[i], num_dead_labels)) {
                continue;
            }
        }
    }

    // need to check if all are dead as could infinite loop forever
    if (num_dead_labels == LargeNumDistanceStates) {
        return false;
    }

    // Typically, we now have LargeNumDistanceStates # of candidates
    // which have phi, psi torsion angles as determined by sampling
    // from empirical distance distribution. The next step is to select
    // a subset of the candidates (specifically a subset of size
    // NumStates) according to the empirical joint distribution of
    // phi, psi angles for the current residue being grown.
    {
        // @TODO - factor this out into a function, it's a bit
        // tricky because of the VLA arrays, but the best candidate
        // refactoring would be Method-to-Class where local variables
        // become the class variables

        // This loop tallies up how many instances of the discretized
        // (phi, psi) pair have been observed. "TorProb" is actually a
        // misnomer (but effectively functions the same way when combined
        // with "Total" variable) as it is actually an observation count
        // and not a normalized probability.
        int Total = 0;
        for (int i = 0; i < LargeNumDistanceStates; ++i) {
            if (TorProb[i] == -1) {
                // Note: 180 / PI converts from radians to degrees
                // Note: -180 degrees maps to bin 0 by adding binNum = 180/BBTbinSize
                //  -> binNum is the number of bins in [0, 180]
                // This gives the bin corresponding to each angle, where each
                //  bin is BBTbinSize degrees.
                phi_int = (int)(larbb_angles[i][eix_phi] / BBTbinSize) + binNum;

                // Decrement bin for negative angles
                // e.g. if phi = -175.5, it will be placed in the bin for
                //  [-175, 170], but really we want it in the bin for [-180, -175].
                //  Note, we shouldn't need to do this for positive angles.
                phi_int -= (larbb_angles[i][eix_phi] / BBTbinSize + binNum - phi_int < 0);

                // Perform same procedure to find bin for psi angle
                psi_int = (int)(larbb_angles[i][eix_psi] / BBTbinSize) + binNum;
                psi_int -= (larbb_angles[i][eix_psi] / BBTbinSize + binNum - psi_int < 0);

                Phi[i] = phi_int;
                Psi[i] = psi_int;

                // Obtain unnormalized observed counts for ([phi], [psi]) intervals
                TorProb[i] = Joint_Angle[CurResType][phi_int][psi_int];
                Total += TorProb[i];
            } // end check for uninitialized torsion probability
        } // end iteration over large num distance states

        // Check if rotamer library usage has been enabled. The phi, psi
        // backbone libraries are histograms with observed phi, psi counts
        // for a training set of several thousand (~10k) loop samples.
        // However, the histograms have not been smoothed and therefore are
        // somewhat sparse and can make it hard to generate clash free loops.
        // In practice, using the rotamer libraries may lead to loop samples
        // with better RMSD to target structures - but with a cost of potentially
        // higher Lennard Jones potential energy.

        // Now actually select a smaller subset of candidates based on empirical
        // (phi, psi) distribution:
        if ((0 == Total) || !use_rot_lib) {
            // This case occurs if the (phi, psi) angles are not empirically
            // observed according to the Ramachandran data (stored in Joint_Angle).
            // In this case, select a random subset and give a random angle
            // in [bin_min, bin_max] for the mapped bin interval.
            for (int i = 0; i < NumStates; ++i) {
                // Select a random bin "uniformly"
                int r = rand() % LargeNumDistanceStates;
                // Skip to next valid sampled state
                while (1 == deadlabel[r])
                    r = ((r + 1) % LargeNumDistanceStates);
                bb_angles[i][eix_phi] = frand((double)((Phi[r] - binNum) * PI / binNum),
                    (double)((Phi[r] - binNum + 1) * PI / binNum));
                bb_angles[i][eix_psi] = frand((double)((Psi[r] - binNum) * PI / binNum),
                    (double)((Psi[r] - binNum + 1) * PI / binNum));
                bb_angles[i][eix_omega] = larbb_angles[r][eix_omega];
            }
        }
        else {
            for (int i = 0; i < NumStates; ++i) {
                // Select a random bin according to empirically observed probability
                int r = rand() % Total;
                int cur = 0;
                for (int j = 0; j < LargeNumDistanceStates; j++) {
                    cur += TorProb[j];
                    if (r < cur) {
                        // Select a random set of torsion angles within selected bin
                        bb_angles[i][eix_phi] = frand((double)((Phi[j] - binNum) * PI / binNum),
                            (double)((Phi[j] - binNum + 1) * PI / binNum));
                        bb_angles[i][eix_psi] = frand((double)((Psi[j] - binNum) * PI / binNum),
                            (double)((Psi[j] - binNum + 1) * PI / binNum));
                        bb_angles[i][eix_omega] = larbb_angles[j][eix_omega];
                        break;
                    }
                }
            }
        }
    } // end select subset of backbone angles

    double dis = 10000, disquare = 10000;

    // Zero out energy arrays
    SET_VLA_1D(double, energy, NumStates, 0);

    // Loop through states, converting sampled angles to positions
    int num_valid_states = 0;
    for (int i_state = 0; i_state < NumStates; ++i_state) {
        tmpRes[num_valid_states] = CurConf._res[Position];
        tmpRes[num_valid_states]._phi = bb_angles[i_state][eix_phi];
        tmpRes[num_valid_states]._psi = bb_angles[i_state][eix_psi];
        tmpRes[num_valid_states]._omega = bb_angles[i_state][eix_omega];
        CurConf.calBBCo(
            Position,
            tmpRes[num_valid_states],
            tmpRes[num_valid_states]._phi,
            tmpRes[num_valid_states]._psi,
            tmpRes[num_valid_states]._omega);
        // Calculate residue center
        CurConf.SinglecalCenter(tmpRes[num_valid_states], 1);

        if (!this->cfe.gobbo_is_bb_clash_free(tmpRes[num_valid_states], CurConf, *this)) {
            continue;
        }

        // calculate energy
        energy[num_valid_states] += one_res_en_loodis_bb2all_list(CurConf, tmpRes[num_valid_states], start, end, 1, Reslist);

        disquare = tmpRes[num_valid_states]._atom[1].disquare(CurConf._res[Position]._atom[0]);

        if (disquare <= PF_DIS_CUT_SQUARE) {
            dis = sqrt(disquare);
            int disInd = (int)(dis / H_INLO);
            energy[num_valid_states] +=
                PF::LOODIS[tmpRes[num_valid_states]._atom[1]._type - 1][CurConf._res[Position]._atom[0]._type - 1][disInd];
        }

        ++num_valid_states;
    }

    // Early out if no states where valid (clash free)
    assert(num_valid_states >= 0);
    assert(num_valid_states <= NumStates);
    if (num_valid_states <= 0) {
        return false;
    }

    double minE = 100000;
    for (int i = 0; i < num_valid_states; ++i) {
        if (energy[i] < minE)
            minE = energy[i];
    }

    // Now find the probabilities of states and select one
    double sum = 0;
    prob.resize(num_valid_states);
    for (int i = 0; i < num_valid_states; ++i) {
        prob[i] = pow(EXPO, double((minE - energy[i]) * 0.5));
        if (std::isinf(prob[i])) {
            cout << "grow_one_bb_only probability error in state " << i << ": "
                << " " << energy[i] << " " << prob[i] << endl;
            exit(0);
        }
        sum += prob[i];
    }

    // in case that all energies are very high, prob[i] are all zero
    // assign all prob[i] to 1
    if (sum <= MIN_BOLTZF) {
        for (int i = 0; i < num_valid_states; ++i) {
            prob[i] = 1;
        }
    }

    const int chosen = SampleOne(prob);

    assert(!CurConf._res[Position]._atom[ATM_N].is_at_origin());
    assert(!CurConf._res[Position]._atom[ATM_CA].is_at_origin());
    assert(CurConf._res[Position]._atom[ATM_C].is_at_origin());
    assert(CurConf._res[Position]._atom[ATM_O].is_at_origin());
    assert(CurConf._res[Position + 1]._atom[ATM_N].is_at_origin());
    assert(CurConf._res[Position + 1]._atom[ATM_CA].is_at_origin());

    // Update potential state if necessary
    this->pfe.on_grow_one_res_chosen(
        CurConf,
        tmpRes[chosen],
        chosen,
        this->starts,
        this->ends,
        Position,
        LoopChosen);

    CurConf._res[Position]._atom[ATM_C].CopyPos(tmpRes[chosen]._atom[ATM_C]);
    CurConf._res[Position]._atom[ATM_O].CopyPos(tmpRes[chosen]._atom[ATM_O]);
    CurConf._res[Position + 1]._atom[ATM_N].CopyPos(tmpRes[chosen]._atom[ATM_N]);
    CurConf._res[Position + 1]._atom[ATM_CA].CopyPos(tmpRes[chosen]._atom[ATM_CA]);

    if (CurConf._res[Position]._type != GLY) {
        assert(CurConf._res[Position]._atom[ATM_CB].is_at_origin());
        CurConf._res[Position]._atom[ATM_CB].CopyPos(tmpRes[chosen]._atom[ATM_CB]);
    }

    CurConf.SinglecalCenter(CurConf._res[Position], 1);
    if (CurConf._res[Position]._type == -1) {
        cout << "grow_one_bb_only atom type error " << Position << endl;
        exit(0);
    }

    return true;
}

// This version grows backbones along with side chains
// A number of candidate conformations for a residue are sampled
// one of them is chosen and its weight is updated.
bool SMC::grow_one_with_sc(Structure &CurConf, const int Position, const int tmpEnd, const int LoopChosen,
                           const Atom &EndPt) {

    // Initialize positional values
    int phi_int, psi_int;

    // Initialize numbers of backbone states
    const int NumBBStates = this->num_dist_states;

    // The current residue
    const int CurResType = CurConf._res[Position]._type;

    // Placeholders for sampling
    DECLARE_VLA_2D(double, bb_angles, NumBBStates, 3);

    const int binNum = 180 / BBTbinSize;
    const double EEdisC = CurConf._res[Position]._atom[ATM_CA].dis(EndPt);
    double EEdisN = 0;

    const int RemLength = tmpEnd - Position;
    const double minEEDisC = minDistcon[1][RemLength];
    const double maxEEDisC = minDistcon[1][RemLength] + 31 * DistconBy[1][RemLength];
    const double minEEDisN = minDistcon[0][RemLength - 1];
    const double maxEEDisN = minDistcon[0][RemLength - 1] + 31 * DistconBy[0][RemLength - 1];
    
    int num_dead_labels = 0;

    /////////////////////////////////////////////////////////////////
    // Sample end-to-end distances and convert to positions.
    // End-to-end distances are sampled from an empirical distribution
    /////////////////////////////////////////////////////////////////

    // Early out if distance constraints are not satisfied
    // for carbon atom
    if (EEdisC < minEEDisC || EEdisC > maxEEDisC) {
        return false;
    }

    // BEGIN LargeNumDistanceStates scope
    {
        const int LargeNumDistanceStates = 160;
        DECLARE_VLA_2D(double, larbb_angles, LargeNumDistanceStates, 3);

        DECLARE_VLA_1D(int, Phi, LargeNumDistanceStates);
        DECLARE_VLA_1D(int, Psi, LargeNumDistanceStates);
        DECLARE_VLA_1D(int, TorProb, LargeNumDistanceStates);
        DECLARE_VLA_1D(int, deadlabel, LargeNumDistanceStates);

        SET_VLA_1D(int, deadlabel, LargeNumDistanceStates, 0);

        // Loop through states
        for (int i = 0; i < LargeNumDistanceStates; i++) {
            // Initialize four residues for this sample
            Residue sam;
            // Sample new points for C atom based on distance to final C
            TorProb[i] = -1;

            // Sample carbon atom
            sample_distance(CurConf._res[Position]._atom[ATM_N],
                            CurConf._res[Position]._atom[ATM_CA],
                            EndPt,
                            Residue::bond_angle[CurResType][ATM_C],
                            Residue::bond_length[CurResType][ATM_C],
                            sam._atom[ATM_C], EEdisC,
                            1, RemLength);

            // Skip if sampling failed (and try again)
            if (conditional_mark_dead_label(sam._atom[ATM_C].x, TorProb[i], deadlabel[i], num_dead_labels)) {
                continue;
            }

            // Sample new points for N atom based on distance to final C
            EEdisN = sam._atom[ATM_C].dis(EndPt);
            if (EEdisN >= minEEDisN && EEdisN <= maxEEDisN) {
                sample_distance(CurConf._res[Position]._atom[ATM_CA],
                                sam._atom[ATM_C],
                                EndPt,
                                Residue::bond_angle[CurResType][ATM_N],
                                Residue::bond_length[CurResType][ATM_N],
                                sam._atom[ATM_N], EEdisN,
                                // -1 since it's NEXT N
                                0, RemLength - 1);
            }
            else {
                TorProb[i] = 0;
                deadlabel[i] = 1;
                ++num_dead_labels;
                continue;
            }

            if (conditional_mark_dead_label(sam._atom[ATM_N].x, TorProb[i], deadlabel[i], num_dead_labels)) {
                continue;
            }

            // Populate remainder of bb_angles with the torsion angles
            // implied by the positional selections for C and N. This is for
            // unity with the existing code base. Omega = pi by default.
            // Note: PI/180 converts from degrees to radians
            larbb_angles[i][eix_phi] =
                    torsion(CurConf._res[Position - 1]._atom[ATM_C],
                            CurConf._res[Position]._atom[ATM_N],
                            CurConf._res[Position]._atom[ATM_CA],
                            sam._atom[ATM_C], torsion_degrees);
            larbb_angles[i][eix_psi] =
                    torsion(CurConf._res[Position]._atom[ATM_N],
                            CurConf._res[Position]._atom[ATM_CA],
                            sam._atom[ATM_C],
                            sam._atom[ATM_N], torsion_degrees);
            // Sample omega (peptide bond) angle from normal distribution
            box_MullerNsample_single(larbb_angles[i][eix_omega], PI, TO_RAD(4));

            // Box muller can return NaN because frand(0,1) can return value > 1 which
            // causes a negative log. This is due to weird FRAND padding with +0.5.
            // The fix is to use the <random> library which handles interval sampling properly.
            conditional_mark_dead_label(larbb_angles[i][eix_omega], TorProb[i], deadlabel[i], num_dead_labels);
        }

        // need to check if all are dead as could infinite loop forever
        if (num_dead_labels == LargeNumDistanceStates) {
            return false;
        }

        /////////////////////////////////////////////////////////////////
        // Of the candidates sampled from end-to-end distance distribution,
        //  select a subset of the generated candidates according to empirical
        //  phi, psi distribution.
        /////////////////////////////////////////////////////////////////

        // Typically, we now have LargeNumDistanceStates # of candidates
        // which have phi, psi torsion angles as determined by sampling
        // from empirical distance distribution. The next step is to select
        // a subset of the candidates (specifically a subset of size
        // NumStates) according to the empirical joint distribution of
        // phi, psi angles for the current residue being grown.
        {
            // @TODO - factor this out into a function, it's a bit
            // tricky because of the VLA arrays, but the best candidate
            // refactoring would be Method-to-Class where local variables
            // become the class variables

            int Total = 0;
            // This loop tallies up how many instances of the discretized
            // (phi, psi) pair have been observed. "TorProb" is actually a
            // misnomer (but effectively functions the same way when combined
            // with "Total" variable) as it is actually an observation count
            // and not a normalized probability.
            for (int i = 0; i < LargeNumDistanceStates; ++i) {

                if (TorProb[i] == -1) {
                    phi_int = (int)(larbb_angles[i][eix_phi] / BBTbinSize) + binNum;
                    psi_int = (int)(larbb_angles[i][eix_psi] / BBTbinSize) + binNum;

                    if (larbb_angles[i][eix_phi] / BBTbinSize + binNum - phi_int < 0)
                        --phi_int;
                    if (larbb_angles[i][eix_psi] / BBTbinSize + binNum - psi_int < 0)
                        --psi_int;
                    Phi[i] = phi_int;
                    Psi[i] = psi_int;

                    TorProb[i] = Joint_Angle[CurResType][phi_int][psi_int];
                    Total += TorProb[i];
                }
            }

            // Check if rotamer library usage has been enabled. The phi, psi
            // backbone libraries are histograms with observed phi, psi counts
            // for a training set of several thousand (~10k) loop samples.
            // However, the histograms have not been smoothed and therefore are
            // somewhat sparse and can make it hard to generate clash free loops.
            // In practice, using the rotamer libraries may lead to loop samples
            // with better RMSD to target structures - but with a cost of potentially
            // higher Lennard Jones potential energy.

            // Now actually select a smaller subset of candidates based on empirical
            // (phi, psi) distribution:
            if ((0 == Total) || !use_rot_lib) {
                for (int i = 0; i < NumBBStates; ++i) {
                    // @TODO - use actual uniform distribution
                    int r = rand() % LargeNumDistanceStates;
                    // Skip to next valid sampled state
                    while (1 == deadlabel[r])
                        r = ((r + 1) % LargeNumDistanceStates);

                    bb_angles[i][eix_phi] = frand((double)((Phi[r] - binNum) * PI / binNum),
                        (double)((Phi[r] - binNum + 1) * PI / binNum));
                    bb_angles[i][eix_psi] = frand((double)((Psi[r] - binNum) * PI / binNum),
                        (double)((Psi[r] - binNum + 1) * PI / binNum));
                    bb_angles[i][eix_omega] = larbb_angles[r][eix_omega];
                    continue;
                }
            }
            else {
                for (int i = 0; i < NumBBStates; ++i) {
                    // Select a random bin
                    int r = rand() % Total;
                    int cur = 0;
                    for (int j = 0; j < LargeNumDistanceStates; j++) {
                        cur += TorProb[j];
                        if (r < cur) {
                            // Select a random set of torsion angles within selected bin
                            bb_angles[i][eix_phi] = frand((double) ((Phi[j] - binNum) * PI / binNum),
                                                    (double) ((Phi[j] - binNum + 1) * PI / binNum));
                            bb_angles[i][eix_psi] = frand((double) ((Psi[j] - binNum) * PI / binNum),
                                                    (double) ((Psi[j] - binNum + 1) * PI / binNum));
                            bb_angles[i][eix_omega] = larbb_angles[j][eix_omega];
                            break;
                        }
                    }
                }
            }
        } // end select subset of backbone angles
    } // END SCOPE for LargeNumDistanceStates (reclaim stack memory for side chain use)

    // Add side chain here:
    // First, sample side chain angles (according to some distribution) to
    // generate numSCStates candidates.
    // Then - for each candidate - perform "clash resolution",
    // finally - score each candidate based on complete one residue
    // potential - and select one candidate stochastically based on score

    // There are SMC::numSCStates candidates that will be generated
    // per backbone candidate
    // Note:
    // - there are 'NumStates' backbone candidates
    // - there are SMC::numSCStates side chains candidates per backbone candidate
    // - there are at most 6 side chain angles per candidate
    // HOPEFULLY THIS WON'T BLOW THE STACK! - ELSE WILL HAVE TO RESORT TO
    // DYNAMIC MEMORY OR INCREASE SYSTEM ALLOCATION TO STACK
    DECLARE_VLA_3D(double, sc_angles, NumBBStates, this->num_sc_states, 6);
    SET_VLA_3D(double, sc_angles, NumBBStates, this->num_sc_states, 6, 0);

    // Sample side chain angles for each backbone candidate
    for (int i = 0; i < NumBBStates; ++i) {
        sample_sc_angles(CurConf._res[Position],
                         sc_angles[i],
                         this->num_sc_states,
                         this->ang_type);
    }

    DECLARE_VLA_2D(Residue, tmpRes, NumBBStates, this->num_sc_states);

    // @TODO - consider side-chain clash resolution prior to potential filtering

    /////////////////////////////////////////////////////////////////
    // Weight each candidate according to potential function
    /////////////////////////////////////////////////////////////////

    // Zero out energy array
    DECLARE_VLA_2D(double, energy, NumBBStates, this->num_sc_states);
    SET_VLA_2D(double, energy, NumBBStates, this->num_sc_states, 0);

    double minE = std::numeric_limits<double>::max();

    // Loop through states, converting sampled angles to positions and
    // compute candidate energy
    {
        int candidate_key = 0;
        for (int i = 0; i < NumBBStates; ++i) {
            for (int j = 0; j < this->num_sc_states; ++j) {
                tmpRes[i][j] = CurConf._res[Position];

                tmpRes[i][j]._phi = bb_angles[i][eix_phi];
                tmpRes[i][j]._psi = bb_angles[i][eix_psi];
                tmpRes[i][j]._omega = bb_angles[i][eix_omega];

                CurConf.calBBCo(Position, tmpRes[i][j],
                                tmpRes[i][j]._phi, tmpRes[i][j]._psi, tmpRes[i][j]._omega);

                // Calculate side chain coordinates based on angles
                CurConf.calSCCo(&(sc_angles[i][j][0]), tmpRes[i][j]);

                // Calculate residue center
                CurConf.SinglecalCenter(tmpRes[i][j], 0 /* type -> include side chains */);

                // Calculate energy
                energy[i][j] =
                        this->pfe.one_res_en_all_list(
                                CurConf,
                                tmpRes[i][j],
                                Reslist,
                                Position,
                                candidate_key,
                                this->starts,
                                this->ends,
                                true, /*b_is_hybrid_res*/
                                LoopChosen
                            );

                // Update minimum energy
                if (energy[i][j] < minE) {
                    minE = energy[i][j];
                }

                ++candidate_key;
            } // end iteration over side chain candidates
        } // end iteration over back bone candidates
    } // end scope for energy ranking of candidates

    /////////////////////////////////////////////////////////////////
    // Compute relative probabilities for each candidate
    /////////////////////////////////////////////////////////////////

    // Now find the probabilities of states and select one
    // Interface expects 1-D array so we have to transform
    // 2-D mapping to 1-D mapping
    std::vector<double> prob(NumBBStates * this->num_sc_states, 0.0);

    double sum = 0;

    for (int i = 0; i < NumBBStates; ++i) {
        for (int j = 0; j < this->num_sc_states; ++j) {
            // Convert 2-D (i,j) indices to 1-D index k
            const int k = i * this->num_sc_states + j;
            // @TODO - consider using exp() instead of pow
            prob[k] = pow(EXPO, double((minE - energy[i][j]) * 0.5));
            assert(std::isfinite(prob[k]));
            if (std::isinf(prob[k])) {
                cout << "grow_one_with_sc probability error in state " << i << ", " << j << ": "
                << " " << energy[i][j] << " " << prob[k] << endl;
                exit(0);
            }
            sum += prob[k];
            assert(std::isfinite(sum));
        }
    }

    // In case that all energies are very high, prob[k] are all zero
    // assign all prob[i] to 1
    if (sum <= MIN_BOLTZF) {
        for (int i = 0; i < NumBBStates; ++i) {
            for (int j = 0; j < this->num_sc_states; ++j) {
                const int k = i * this->num_sc_states + j;
                prob[k] = 1;
            }
        }
    }

    // Select final candidate based on probability weight
    const int chosen_scalar = SampleOne(prob);
    // convert 1-D index back to 2-D index:
    const int chosen_i = chosen_scalar / this->num_sc_states;
    const int chosen_j = chosen_scalar % this->num_sc_states;

    if ((CurConf._res[Position]._type == -1) || (CurConf._res[Position]._type == UNDEF)) {
        cout << "grow_one_bb_only atom type error " << Position << endl;
        exit(0);
    }

    // Copy chosen candidate's coordinates
    const Residue &chosen_res = tmpRes[chosen_i][chosen_j];

    assert(!CurConf._res[Position]._atom[ATM_N].is_at_origin());
    assert(!CurConf._res[Position]._atom[ATM_CA].is_at_origin());
    assert(CurConf._res[Position]._atom[ATM_C].is_at_origin());
    assert(CurConf._res[Position]._atom[ATM_O].is_at_origin());
    assert(CurConf._res[Position + 1]._atom[ATM_N].is_at_origin());
    assert(CurConf._res[Position + 1]._atom[ATM_CA].is_at_origin());

    // Update potential state if necessary
    this->pfe.on_grow_one_res_chosen(CurConf,
                                     chosen_res,
                                     chosen_scalar,
                                     this->starts,
                                     this->ends,
                                     Position,
                                     LoopChosen);

    CurConf._res[Position]._atom[ATM_C].CopyPos(chosen_res._atom[ATM_C]);
    CurConf._res[Position]._atom[ATM_O].CopyPos(chosen_res._atom[ATM_O]);
    CurConf._res[Position + 1]._atom[ATM_N].CopyPos(chosen_res._atom[ATM_N]);
    CurConf._res[Position + 1]._atom[ATM_CA].CopyPos(chosen_res._atom[ATM_CA]);

    if (CurConf._res[Position]._type != GLY) {
        assert(CurConf._res[Position]._atom[ATM_CB].is_at_origin());
        CurConf._res[Position]._atom[ATM_CB].CopyPos(chosen_res._atom[ATM_CB]);

        // Copy side chain positions
        assert(chosen_res._numAtom == CurConf._res[Position]._numAtom);
        for (int i = NUM_BB_ATOM; i < CurConf._res[Position]._numAtom; ++i) {
            CurConf._res[Position]._atom[i].CopyPos(chosen_res._atom[i]);
        }
    }

    // Update center calculations
    CurConf.SinglecalCenter(CurConf._res[Position], 0 /* type -> include side chains */);

    return true;
}

// For glycine and alanine, only backbone atoms are grown (CB is backbone), for all other
// residue types, side chain atoms are also grown.
bool SMC::grow_one(const int Position, const int tmpEnd, const int LoopChosen) {
    bool grow_Success = false;
    // Allow incremental side chain growth for null_with_sc or electrostatic enabled builds
#ifdef DISGRO_BUILD_ENABLE_SC_GROW_ONE
    // grow_one_with_sc adds side chains, however, if the residue is
    // glycine (no side chain) or alanine (only methyl group), then SCRMAP doesn't have proper angles
    // so we avoid side chain sampling in this case.
    if (_conf._res[Position]._type != GLY && _conf._res[Position]._type != ALA) {
        grow_Success = grow_one_with_sc(_conf, Position, tmpEnd, LoopChosen, _conf._res[tmpEnd]._atom[ATM_C]);
    }
    else
#endif // DISGRO_BUILD_ENABLE_SC_GROW_ONE
    {
        grow_Success = grow_one_bb_only(_conf, Position, tmpEnd, LoopChosen, _conf._res[tmpEnd]._atom[ATM_C]);
    }
    return grow_Success;
}

// Marks atoms at start residue as being at origin
// but avoids marking N and CA atoms as they are given
void Structure::mark_start_residue(const int start) {
    // For start residue, N and CA are known
    Residue& start_res = this->_res[start];
    start_res._atom[ATM_C].move_to_origin();
    start_res._atom[ATM_O].move_to_origin();
    if (start_res._type != GLY) {
        start_res._atom[ATM_CB].move_to_origin();
        for (int iAtom=NUM_BB_ATOM; iAtom < start_res._numAtom; ++iAtom) {
            start_res._atom[iAtom].move_to_origin();
        }
    }
}

// Marks atoms at end residue as being at origin
// but avoids marking CA and C atoms as they are given
void Structure::mark_end_residue(const int end) {
    // For end residue, CA and C are known
    Residue& end_res = this->_res[end];
    end_res._atom[ATM_N].move_to_origin();
    end_res._atom[ATM_O].move_to_origin();
    if (end_res._type != GLY) {
        end_res._atom[ATM_CB].move_to_origin();
        for (int iAtom = NUM_BB_ATOM; iAtom < end_res._numAtom; ++iAtom) {
            end_res._atom[iAtom].move_to_origin();
        }
    }
}

void Structure::mark_loop_region(const int start, const int end) {
    // start may be equal to end for certain utilities
    assert(start <= end);
    assert(start >= 1 && start < this->_numRes);
    assert(end >= 1 && end <= this->_numRes);

    // For start residue, N and CA are known
    this->mark_start_residue(start);

    // For end residue, CA and C are known
    this->mark_end_residue(end);

    // For inner loop residues, all atoms are unknown
    for (int iRes=start+1; iRes<=end-1; ++iRes) {
        Residue& res = this->_res[iRes];
        for (int iAtom=0; iAtom < res._numAtom; ++iAtom) {
            res._atom[iAtom].move_to_origin();
        }
    }
}

// Move atoms in loop regions to origin (0,0,0). This "marks" atoms that have not yet
//  been grown and avoids their use in potential energy calculations.
void Structure::mark_loop_regions(const std::vector<int> &sTart, const std::vector<int> &eNd) {
    assert(sTart.size() == eNd.size());
    for (size_t iLoop=0; iLoop<sTart.size(); ++iLoop) {
        this->mark_loop_region(sTart[iLoop], eNd[iLoop]);
    }
    calCenter(1, _numRes);
}

// Debug function for printing status of backbone atoms in interval [start, end]
void Structure::print_bb_atoms_status(const int start, const int end) {
    assert(start <= end);
    assert(start >= 1 && start <= this->_numRes);
    assert(end >= 1 && end <= this->_numRes);

    for(int i_res = start; i_res <= end; ++i_res) {
        const Residue &res = this->_res[i_res];
        cout << "Res " << i_res << ": ";
        cout << "N " << res._atom[ATM_N].is_at_origin() << " ";
        cout << "CA " << res._atom[ATM_CA].is_at_origin() << " ";
        cout << "C " << res._atom[ATM_C].is_at_origin() << " ";
        cout << "Ox " << res._atom[ATM_O].is_at_origin() << std::endl;
    }
}
