/**
 * Utilities for clash detection with broad-phase acceleration against static
 * geometry.
 */

#include "build.h"
#include "structure.h"
#include "smc.h"
#include "collision_frontend.h"
#include "params.h"

/**
 * Utilities common to enabled and disabled clash-free implementations
 */
namespace {

    /**
     * Utility for applying the selected fragment to the currently grown conformation
     */
    void frag_apply(const int select_lib_idx,
        const int LoopChosen,
        vector<loop_info> &multiloops,
        const vector<int> &multilooplength,
        const int tail_len,
        SMC &smc) {
        // Apply fragment
        const int Original_start = multiloops[LoopChosen].Start;
        for (int i = 1; i <= multilooplength[LoopChosen] - tail_len; ++i) {
            const int res_idx = Original_start + i - 1;
            smc._conf._res[res_idx] = smc.MultiLooplib[LoopChosen][select_lib_idx]._res[i];
            assert(smc._conf._res[res_idx]._posn == smc.Conf._res[res_idx]._posn);
        }
        assert((Original_start + multilooplength[LoopChosen] - tail_len - 1) == (smc.end - tail_len));
        multiloops[LoopChosen].CurrentPos = smc.end - tail_len;
        smc._conf.mark_loop_region(multiloops[LoopChosen].CurrentPos, smc.end);
        assert(Original_start == smc.starts[LoopChosen]);
        assert(multiloops[LoopChosen].End == smc.ends[LoopChosen]);
        smc.pfe.on_restart_from_fragment(smc._conf, LoopChosen, Original_start,
            multiloops[LoopChosen].End, 1,
            multilooplength[LoopChosen] - tail_len);
        assert((multiloops[LoopChosen].CurrentPos - 1) > Original_start);
#ifdef DISGRO_BUILD_ENABLE_SC_GROW_ONE
        // Append side chains for library up to second to last residue added
        smc._conf.grow_sc(Original_start, multiloops[LoopChosen].CurrentPos - 1, 0, smc.num_sc_states, smc.ang_type,
            smc._conf._toBeSampled, smc, LoopChosen);
#endif // DISGRO_BUILD_ENABLE_SC_GROW_ONE
    }

    /**
     * Restarts a loop region without using a fragment library
     */
    void restart_loop_no_frag(const int LoopChosen,
        vector<loop_info> &multiloops,
        SMC &smc) {
        const int Original_start = multiloops[LoopChosen].Start;
        multiloops[LoopChosen].CurrentPos = multiloops[LoopChosen].Start;
        smc._conf.mark_loop_region(multiloops[LoopChosen].Start, multiloops[LoopChosen].End);
        // @TODO - verify this works: if we don't have a fragment library, still need to clear the fragment contribution
        smc.pfe.on_restart_from_fragment(smc._conf, LoopChosen, Original_start,
            multiloops[LoopChosen].End, Original_start, Original_start);
    }

} // end anonymous namespace

#ifdef DISGRO_BUILD_ENABLE_CLASH_FREE_BACKBONE

#include "vdw_utils.h"

#ifndef DISGRO_BUILD_ENABLE_SELF_CLASH_CHECKS
    //#define DISGRO_BUILD_ENABLE_CB_CLASH_CHECKS
#endif // DISGRO_BUILD_ENABLE_SELF_CLASH_CHECKS

namespace {

    /**
     * @return 1 if atoms overlap, 0 otherwise
     */
    inline int atoms_overlap(const Atom &a, const Atom &b, const double ofac2) {
        const double a_r_vdw = collision_grid_prot::elem_radius(a);
        const double b_r_vdw = collision_grid_prot::elem_radius(b);
        const double disquare = a.disquare(b);
        double sum_radii_square = a_r_vdw + b_r_vdw;
        sum_radii_square *= sum_radii_square;
        return (disquare < ofac2*sum_radii_square);
    }

    /**
     * WARNING: Hard-coded values based on atomProp2.txt
     * If order of atoms, particularly backbone atoms, changes in "# restype"
     # table, then tables here need to change to reflect new ordering
     */

    // For a mixed residue, as is the case for a residue being grown in
    // grow_one_bb_only(), atoms N and CA actually belong to residue i+1.
    // This table maps the back bone atom to its residue offset.
    const int bb_mixed_residue_offset[NUM_BB_ATOM] = {
        // N  CA  C  O  H  CB
           1, 1,  0, 0, 0, 0
    };

    // 1-3 mapping for backbone atoms within self or adjacent residues.
    // Table is 1 if atoms are within 2-bonds of each other, 0 o/w
    // <atom_name>+1 is for the atom in the adjacent residue
    // e.g. N+1 is for backbone nitrogen of adjacent residue
    const int bb_1_3[NUM_BB_ATOM][2 * NUM_BB_ATOM] = {
        //        N  CA  C  O  H  CB N+1 CA+1 C+1 O+1 H+1 CB+1
        /*N:*/  { 1, 1,  1, 0, 1, 1, 0,  0,   0,  0,  0,  0 },
        /*CA:*/ { 1, 1,  1, 1, 1, 1, 1,  0,   0,  0,  0,  0 },
        /*C:*/  { 1, 1,  1, 1, 1, 1, 1,  1,   0,  0,  0,  0 },
        /*O:*/  { 0, 1,  1, 1, 0, 0, 1,  0,   0,  0,  0,  0 },
        /*H:*/  { 1, 1,  1, 0, 1, 1, 0,  0,   0,  0,  0,  0 },
        /*CB:*/ { 1, 1,  1, 0, 1, 1, 0,  0,   0,  0,  0,  0 },
    };

    /**
     * Enumerates atom_1_3_info parameters
     */
    enum e1_3 {
        // Info index for atom being belonging to candidate residue being
        // grown via grow_one_bb_only() call
        e1_3_gobbo=0,
        // Info index for atom not currently being grown
        e1_3_alt,
        e1_3_num
    };

    /**
     * Method assumes the following:
     *  - infos[e1_3_gobbo] is for an atom within a candidate residue not yet
     *      added to structure and which is missing a side chain. Therefore,
     *      the atom within this info must be a backbone atom.
     * @return 1 if parameter atoms are within two bonds or less of each
     *  other, 0 otherwise
     */
    inline int gobbo_are_1_3_connected(const struct atom_info infos[e1_3_num]) {
        // Make sure infos[] is expected size
        assert(e1_3_num == 2);
        // Make sure gobbo info is for a backbone atom
        assert(infos[e1_3_gobbo].atom_ix >= 0);
        assert(infos[e1_3_gobbo].atom_ix < NUM_BB_ATOM);

        // Determine min, max residue indices
        const int i_max_res = (infos[1].res_ix > infos[0].res_ix);
        const int i_min_res = 1 - i_max_res;

        //****************************************
        // Case: Residues must be self or adjacent
        //****************************************

        const int delta_res = infos[i_max_res].res_ix - infos[i_min_res].res_ix;
        if (delta_res > 1) {
            // Not self or adjacent, cannot be 1-3 connected
            return 0;
        }

        //****************************************
        // Case: Comparing two backbone atoms
        //****************************************

        if (infos[e1_3_alt].atom_ix < NUM_BB_ATOM) {
            // Can use a simple look-up table
            // delta_res in {0,1}: offset for self residue (0) vs adjacent residue (1)
            const int ix_min_bb_atom = infos[i_min_res].atom_ix;
            const int ix_max_bb_atom = infos[i_max_res].atom_ix + (delta_res * NUM_BB_ATOM);
            return bb_1_3[ix_min_bb_atom][ix_max_bb_atom];
        }

        //****************************************
        // Case: adjacent i+1 proline residue
        //****************************************

        // If we've reached this point: since gobbo info is assumed to not
        // have a side chain, then we must be comparing a side chain from an
        // adjacent residue to the gobbo backbone atom. The only side chain
        // atom which could be 1-3 connected to an adjacent backbone atom is
        // the CG2 atom of proline which connects to its own N atom.
        assert(delta_res == 1);
        return
            // Only care if alt residue is proline
            (infos[e1_3_alt].res_type == PRO) &&
            // and it's the (i+1)-th residue relative to gobbo residue
            (e1_3_alt == i_max_res) &&
            // Only atom CG2 (index=7) connects to N atom
            (7  /*CG2*/ == infos[e1_3_alt].atom_ix) &&
            // Only C atom of gobbo residue i is within 2 bonds
            (infos[e1_3_gobbo].atom_ix == ATM_C);
    }

    /**
     * @param atom - the atom belonging to residue currently being grown
     *  via grow_one_bb_only() call. Assumed to be a mixed residue
     *  containing only backbones atoms where N and CA are from residue i+1
     *  and all other backbone atoms are from residue i. Contains no sidechain.
     * @param anchor_res - the start or end residue of a loop region
     * @params ignore_atom_ix_a, ignore_atom_ix_b - For start residue,
     *  N and CA are assumed known. For end residue CA and C are assumed known.
     *  In either case, we ignore clash detection in this phase as broad phase
     *  has already handled it.
     */
    inline int gobbo_overlaps_with_anchor_residue_bb(
        const Atom& atom,
        const Residue &anchor_res,
        struct atom_info infos[e1_3_num],
        const int ignore_atom_ix_a,
        const int ignore_atom_ix_b,
        const CollisionFrontend &cfe) {
    
        infos[e1_3_alt].res_ix = anchor_res._posn;
        infos[e1_3_alt].res_type = anchor_res._type;
        const double ofac2 = cfe.get_ofac2(infos[e1_3_gobbo].res_ix, infos[e1_3_alt].res_ix);
        for (int i_sim_atom = 0; i_sim_atom < NUM_BB_ATOM; ++i_sim_atom)
        {
            const Atom &other_atom = anchor_res._atom[i_sim_atom];
            if ((i_sim_atom == ignore_atom_ix_a) || (i_sim_atom == ignore_atom_ix_b)) {
                continue;
            }
            if (VdwUtils::should_ignore_atom(other_atom)) {
                continue;
            }
            // skip if 1-3 connected (connected within at most 2 bonds)
            infos[e1_3_alt].atom_ix = i_sim_atom;
            if (gobbo_are_1_3_connected(infos)) {
                continue;
            }
            if (atoms_overlap(atom, other_atom, ofac2)) {
                COLLISION_LOGGER_RECORD(cfe.clog, atom, other_atom, infos, "anchor");
                return 1;
            }
        } // end iteration over anchor residue atoms

        // No overlaps detected
        return 0;
    }

    /**
     * @return true if assumptions about a candidate residue
     *  grown via grow_one_bb_only() are met, false o/w
     */
    bool verify_is_gobbo_candidate_res(const Residue &res) {
        // Verify backbone exists
        for (int i = 0; i < NUM_BB_ATOM; ++i) {
            // Skip hydrogen
            if (i == ATM_H) {
                continue;
            }
            // Glycine
            if (i == ATM_CB && res._type == GLY) {
                continue;
            }
            const Atom &atom = res._atom[i];
            if (atom.is_at_origin()) {
                assert(false);
                return false;
            }
        }
        // Verify side chain is missing
        for (int i = NUM_BB_ATOM; i < res._numAtom; ++i) {
            const Atom &atom = res._atom[i];
            if (!atom.is_at_origin()) {
                assert(false);
                return false;
            }
        }
        return true;
    }

    /**
     * Checks assumption that if ATM_N is at origin, then all atoms at
     * parameter residue are also at origin (this is not valid
     * assumption for anchor residues).
     */
    bool verify_gobbo_int_res_missing_assumptions(const Residue &res) {
        if (res._atom[ATM_N].is_at_origin()) {
            for (int i = 0; i < res._numAtom; ++i) {
                const Atom &atom = res._atom[i];
                if (!atom.is_at_origin()) {
                    assert(false);
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Can only be called for an intermediate residue within a simulated loop
     * region grown via grow_one_bb_only (acronym: gobbo) calls. Will check if
     * first ATM_N is missing and assume that remainder of residue is missing
     * as well.
     */
    inline int gobbo_is_int_res_missing(const Residue &res) {
        assert(verify_gobbo_int_res_missing_assumptions(res));
        return res._atom[ATM_N].is_at_origin();
    }

    /**
     * Ignores non-simulated collision checks as fragment backbones
     * are assumed to be clash free.
     * DOES NOT CHECK AGAINST SIMULATED SIDE CHAINS
     * @param LoopChosen - the loop being regrown from a fragment
     * @param frag_idx - the index within the fragment library containing
     *  the fragment to check for collisions
     * @param cfe - collision front end interface
     * @param smc - simulation object
     * @return true if fragment does not collides with any simulated loop
     *  backbone atoms not within the LoopChosen region, false otherwise
     */
    bool frag_is_bb_clash_free(const int LoopChosen,
                               const int frag_idx,
                               const int frag_end,
                               const CollisionFrontend &cfe,
                               const SMC &smc) {

        const size_t num_loops = smc.starts.size();

        // Structure used for detecting whether two atoms are within two bonds or
        // less of each other
        struct atom_info infos[e1_3_num] = { 0 };

        // The overlap factor: If the ratio of distance between two atoms centers
        // to the sum of their atomic radii is less than this value, then the atoms
        // are colliding. Use non-adjacent overlap factor.
        const double ofac2 = cfe.get_ofac2(CollisionFrontend::eIxOFAC_nadj);

        // Iterate over fragment residues
        assert(frag_end <= smc.MultiLooplib[LoopChosen][frag_idx]._numRes);
        for (int i_fres = 1; i_fres <= frag_end; ++i_fres) {

            const Residue &res = smc.MultiLooplib[LoopChosen][frag_idx]._res[i_fres];

            // Iterate over atoms for fragment residue
            for (int i_fatom = 0; i_fatom < NUM_BB_ATOM; ++i_fatom) {

                const Atom &atom = res._atom[i_fatom];

                if (VdwUtils::should_ignore_atom(atom)) {
                    continue;
                }

                // These are not mixed residues, so do not offset residue index!
                infos[e1_3_gobbo].res_ix = res._posn;
                assert(infos[e1_3_gobbo].res_ix >= 1);
                assert(infos[e1_3_gobbo].res_ix <= smc._conf._numRes);
                infos[e1_3_gobbo].atom_ix = i_fatom;
                infos[e1_3_gobbo].res_type = smc._conf._res[infos[e1_3_gobbo].res_ix]._type;

                //****************************************************
                // narrow phase check against simulated protein atoms
                //****************************************************

                // Iterate over all (simulated) loop regions:
                for (size_t i_loop = 0; i_loop < num_loops; ++i_loop) {
                    assert(smc.starts[i_loop] < smc.ends[i_loop]);

                    // Skip self loop checks
                    if (LoopChosen == i_loop) {
                        continue;
                    }

                    //****************************************************
                    // check start residue of loop region
                    //****************************************************

                    if (gobbo_overlaps_with_anchor_residue_bb(
                        atom,
                        smc._conf._res[smc.starts[i_loop]],
                        infos,
                        ATM_N,
                        ATM_CA,
                        cfe)) {
                        return false;
                    }

                    //****************************************************
                    // check end residue of loop region
                    //****************************************************

                    if (gobbo_overlaps_with_anchor_residue_bb(
                        atom,
                        smc._conf._res[smc.ends[i_loop]],
                        infos,
                        ATM_CA,
                        ATM_C,
                        cfe)) {
                        return false;
                    }

                    //****************************************************
                    // check intermediate residues of loop region
                    //****************************************************

                    for (size_t i_res = smc.starts[i_loop] + 1; i_res < smc.ends[i_loop]; ++i_res) {
                        const Residue &other_res = smc._conf._res[i_res];
                        // Assuming consecutive residue growth, we can skip remainder
                        // of loop if anchor N atom is missing. 
                        if (gobbo_is_int_res_missing(other_res)) {
                            break;
                        }
                        infos[e1_3_alt].res_ix = other_res._posn;
                        infos[e1_3_alt].res_type = other_res._type;
                        for (int i_sim_atom = 0; i_sim_atom < NUM_BB_ATOM; ++i_sim_atom) {
                            const Atom &other_atom = other_res._atom[i_sim_atom];
                            if (VdwUtils::should_ignore_atom(other_atom)) {
                                continue;
                            }
                            // No need to check 1-3 connectivity as
                            // these atoms should not be adjacent
                            infos[e1_3_alt].atom_ix = i_sim_atom;
                            if (atoms_overlap(atom, other_atom, ofac2)) {
                                return false;
                            }
                        } // end iteration over intermediate residue atoms
                    } // end iteration over intermediate residues
                } // end iteration over simulated loop regions
            } // end iteration over fragment backbone atoms
        } // end iteration over fragment residues

        // No collision with simulated backbone detected
        return true;
    }

    /**
     * HACK: Ca+1 and carbonyl O always seem to clash, even in native structure,
     * so ignore that specific interaction.
     * @return 1 if parameter atoms are Ca+1 and carbonyl Oxygen, 0 otherwise
     */
    int hack_ignore_Ca1_with_O_clash(const struct atom_info infos[e1_3_num]) {
        // Make sure infos[] is expected size
        assert(e1_3_num == 2);

        // Determine min, max residue indices
        const int i_max_res = (infos[1].res_ix > infos[0].res_ix);
        const int i_min_res = 1 - i_max_res;

        const int delta_res = infos[i_max_res].res_ix - infos[i_min_res].res_ix;
        if (delta_res == 1) {
            return (infos[i_max_res].atom_ix == ATM_CA) &&
                (infos[i_min_res].atom_ix == ATM_O);
        }

        return 0;
    }

} // end of anonymous namespace

/**
 * Initializes broad phase
 * @param markedConf - simulated loop regions must be moved to origin
 * @param params - command line arguments
 */
void CollisionFrontend::init(const Structure &markedConf, const struct Params &params) {
    // Initialize broad phase collision checks for static geometry
    m_broad_phase.init(
          markedConf
        , collision_grid_protein_heuristics::max_elem_diameter()
        , collision_grid_protein_heuristics::find_min_coord_for_grid(markedConf, params.starts, params.ends)
        );
    m_ofac2[eIxOFAC_adj] = params.ofa * params.ofa;
    assert(m_ofac2[eIxOFAC_adj] > 0.0);
    m_ofac2[eIxOFAC_nadj] = params.ofna * params.ofna;
    assert(m_ofac2[eIxOFAC_nadj] > 0.0);
    m_ofac2[eIxOFAC_memb] = params.ofmemb * params.ofmemb;
    assert(m_ofac2[eIxOFAC_memb] > 0.0);
    printf("Collision grid initialized with:\n");
    printf("\tadjacent overlap factor: %f\n", params.ofa);
    printf("\tnon-adjacent overlap factor: %f\n", params.ofna);
    printf("\tmembrane overlap factor: %f\n", params.ofmemb);
    // Initialize membrane
    m_membrane.init(params.memb_file);
}

/**
 * Prints clash statistics against broad phase only
 * @param Conf - the structure to check for clashes at loop regions
 * @param starts - vector of start residue indices for loop regions
 * @param ends - vector of end residue indices for loop regions
 * @param b_ignore_self - if true, then assume that broad phase was
 *  generated by Conf and 'self' collisions are ignored
 */
void CollisionFrontend::print_broad_clash_counts(const class Structure &Conf,
    const std::vector<int> &starts,
    const std::vector<int> &ends,
    const bool b_ignore_self) const {
 
    assert(starts.size() == ends.size());   
    std::cout << "res,bbc,scc,tot,mc\n";
    
    // Structure used for detecting whether two backbone atoms are within two bonds or
    // less of each other
    struct atom_info infos[e1_3_num] = { 0 };

    // The overlap factor: If the ratio of distance between two atoms centers
    // to the sum of their atomic radii is less than this value, then the atoms
    // are colliding.
    double ofac2 = 1.0;
    
    const size_t num_loops = starts.size();

    // Total number of loop backbone atoms that clash with protein
    int total_bb_clashes = 0;
    // Total number of loop side chain atoms that clash with protein
    int total_sc_clashes = 0;
    // Total number of loop atoms that clash with membrane
    int total_memb_clashes = 0;

    for (size_t i_loop = 0; i_loop < num_loops; ++i_loop) {
        for (size_t i_res = starts[i_loop]; i_res <= ends[i_loop]; ++i_res) {
            const Residue &res = Conf._res[i_res];
            assert(res._posn == i_res);

            // Total number of backbone atoms at this residue that clash with protein
            int bb_clashes = 0;
            // Total number of side chain atoms at this residue that clash with protein
            int sc_clashes = 0;
            // Total number of atoms at this residue that clash with membrane
            int memb_clashes = 0;

            for (short i_atom = 0; i_atom < res._numAtom; ++i_atom) {
                // Obtain handle to atom
                const Atom &atom = res._atom[i_atom];
                if (VdwUtils::should_ignore_atom(atom)) {
                    continue;
                }

                // 1-3 collision filter for backbone atoms
                const bool is_bb_atom = (i_atom < NUM_BB_ATOM);
                if (is_bb_atom) {
                    infos[e1_3_gobbo].res_ix = static_cast<int>(i_res);
                    assert(infos[e1_3_gobbo].res_ix >= 1);
                    assert(infos[e1_3_gobbo].res_ix <= Conf._numRes);
                    infos[e1_3_gobbo].atom_ix = i_atom;
                    infos[e1_3_gobbo].res_type = Conf._res[infos[e1_3_gobbo].res_ix]._type;
                }

                // Check against membrane
                if (m_membrane.collides(atom, m_ofac2[eIxOFAC_memb])) {
                    ++memb_clashes;
                }

                // Broad phase
                collision_grid_prot::filter_results_t nearby_atoms;
                m_broad_phase.filter(nearby_atoms,
                                     atom.x, atom.y, atom.z,
                                     collision_grid_prot::elem_radius(atom));

                // Narrow phase
                const size_t num_nearby_atoms = nearby_atoms.size();
                for (size_t i_nearby_atom = 0; i_nearby_atom < num_nearby_atoms; ++i_nearby_atom) {
                    const collision_grid_prot::atom_key &nearby_atom_key = nearby_atoms[i_nearby_atom];
                    assert(nearby_atom_key.res_ix >= 1);
                    assert(nearby_atom_key.res_ix <= static_cast<unsigned int>(Conf._numRes));
                    assert(nearby_atom_key.atom_ix >= 0);
                    assert(nearby_atom_key.atom_ix < static_cast<unsigned int>(Conf._res[nearby_atom_key.res_ix]._numAtom));
                    const Atom &nearby_atom = Conf._res[nearby_atom_key.res_ix]._atom[nearby_atom_key.atom_ix];
                    assert(!VdwUtils::should_ignore_atom(nearby_atom));

                    // Skip self residues
                    if (b_ignore_self && (nearby_atom_key.res_ix == i_res)) {
                        continue;
                    }

                    // Skip if 1-3 connected (connected within at most 2 bonds)
                    // only supported if first atom is backbone atom
                    if (is_bb_atom)
                    {
                        infos[e1_3_alt].res_ix = nearby_atom_key.res_ix;
                        infos[e1_3_alt].atom_ix = nearby_atom_key.atom_ix;
                        infos[e1_3_alt].res_type = Conf._res[nearby_atom_key.res_ix]._type;
                        if (gobbo_are_1_3_connected(infos)) {
                            continue;
                        }

                        // HACK: Ca+1 and carbonyl O always seem to clash,
                        // even in native structure so ignore that specific interaction
                        if (hack_ignore_Ca1_with_O_clash(infos)) {
                            continue;
                        }
                    }

                    // Check for overlap
                    ofac2 = get_ofac2(static_cast<int>(i_res), static_cast<int>(nearby_atom_key.res_ix));
                    if (atoms_overlap(atom, nearby_atom, ofac2)) {
                        if (is_bb_atom) {
                            ++bb_clashes;
                        }
                        else {
                            ++sc_clashes;
                        }
                        // Check next atom for clashes
                        break;
                    }
                } // end narrow phase check

            } // end iteration over atoms

            std::cout << res.get_name_3() << ":" << i_res << ","
                      << bb_clashes << ","
                      << sc_clashes << ","
                      << bb_clashes + sc_clashes << ","
                      << memb_clashes << std::endl;

            total_bb_clashes += bb_clashes;
            total_sc_clashes += sc_clashes;
            total_memb_clashes += memb_clashes;

        } // end iteration over residues
    } // end iteration over loops

    std::cout << "TOTAL BB CLASHES: " << total_bb_clashes << std::endl;
    std::cout << "TOTAL SC CLASHES: " << total_sc_clashes << std::endl;
    std::cout << "TOTAL MEMB CLASHES: " << total_memb_clashes << std::endl;
}

/**
 * gobbo - *g*row *o*ne *b*ack*b*one only
 * This method is only meant to be called from grow_one_bb_only() as
 * as simplifying assumptions made are only valid for this use case.
 * Assumptions:
 *  -param res is a mixed residue containing only backbone atoms
 *      (with no side chain) and where atoms N and CA are for residue i+1
 *      and remaining backbone atoms are for residue i
 *  -simulated loop regions in param Conf can be ignored once the first
 *      non-hydrogen atom placed at origin is encountered
 *  -simulated loop regions have no side chains placed
 * @return true if residue backbone is clash free, false o/w
 */
bool CollisionFrontend::gobbo_is_bb_clash_free(
    const Residue &res,
    const Structure &Conf,
    const SMC &smc) const {
    // Verify parameter residue is candidate from gobbo call
    assert(verify_is_gobbo_candidate_res(res));
    // Verify loop interval arrays are same size
    assert(smc.starts.size() == smc.ends.size());
    // Verify residue type matches expected residue type
    assert(res._type == Conf._res[res._posn]._type);

    const size_t num_loops = smc.starts.size();

    // Structure used for detecting whether two atoms are within two bonds or
    // less of each other
    struct atom_info infos[e1_3_num] = { 0 };

    // The overlap factor: If the ratio of distance between two atoms centers
    // to the sum of their atomic radii is less than this value, then the atoms
    // are colliding.
    double ofac2 = 1.0;

    for (int i_atom = 0; i_atom < NUM_BB_ATOM; ++i_atom) {
        const Atom &atom = res._atom[i_atom];
        if (VdwUtils::should_ignore_atom(atom)) {
            continue;
        }
        infos[e1_3_gobbo].res_ix = res._posn + bb_mixed_residue_offset[i_atom];
        assert(infos[e1_3_gobbo].res_ix >= 1);
        assert(infos[e1_3_gobbo].res_ix <= Conf._numRes);
        infos[e1_3_gobbo].atom_ix = i_atom;
        infos[e1_3_gobbo].res_type = Conf._res[infos[e1_3_gobbo].res_ix]._type;

        //****************************************************
        // broad phase check against unsimulated protein atoms
        //****************************************************

        collision_grid_prot::filter_results_t nearby_atoms;
        m_broad_phase.filter(nearby_atoms,
                             atom.x, atom.y, atom.z,
                             collision_grid_prot::elem_radius(atom));

        //****************************************************
        // narrow phase check against unsimulated protein atoms
        //****************************************************

        const size_t num_nearby_atoms = nearby_atoms.size();
        for (size_t i_unsim_atom = 0; i_unsim_atom < num_nearby_atoms; ++i_unsim_atom) {
            const collision_grid_prot::atom_key &nearby_atom_key = nearby_atoms[i_unsim_atom];
            assert(nearby_atom_key.res_ix >= 1);
            assert(nearby_atom_key.res_ix <= static_cast<unsigned int>(Conf._numRes));
            assert(nearby_atom_key.atom_ix >= 0);
            assert(nearby_atom_key.atom_ix < static_cast<unsigned int>(Conf._res[nearby_atom_key.res_ix]._numAtom));
            const Atom &nearby_atom = Conf._res[nearby_atom_key.res_ix]._atom[nearby_atom_key.atom_ix];
            assert(!VdwUtils::should_ignore_atom(nearby_atom));
            // skip if 1-3 connected (connected within at most 2 bonds)
            infos[e1_3_alt].res_ix = nearby_atom_key.res_ix;
            infos[e1_3_alt].atom_ix = nearby_atom_key.atom_ix;
            infos[e1_3_alt].res_type = Conf._res[nearby_atom_key.res_ix]._type;
            if (gobbo_are_1_3_connected(infos)) {
                continue;
            }
            ofac2 = get_ofac2(infos[e1_3_gobbo].res_ix, infos[e1_3_alt].res_ix);
            if (atoms_overlap(atom, nearby_atom, ofac2)) {
                COLLISION_LOGGER_RECORD(clog, atom, nearby_atom, infos, "static");
                return false;
            }
        } // end narrow phase check

        //****************************************************
        // check against membrane atoms
        //****************************************************

        if (this->m_membrane.collides(atom, m_ofac2[eIxOFAC_memb])) {
            // @TODO - enable more detailed logging of membrane collisions,
            //  perhaps move to membrane method
            COLLISION_LOGGER_RECORD(clog, atom, atom, infos, "memb");
            return false;
        }

        //****************************************************
        // narrow phase check against simulated protein atoms
        //****************************************************

        // Iterate over all (simulated) loop regions:
        for (size_t i_loop = 0; i_loop < num_loops; ++i_loop) {
            assert(smc.starts[i_loop] < smc.ends[i_loop]);

            //****************************************************
            // check start residue of loop region
            //****************************************************

            if (gobbo_overlaps_with_anchor_residue_bb(
                atom,
                Conf._res[smc.starts[i_loop]],
                infos,
                ATM_N,
                ATM_CA,
                *this)) {
                return false;
            }

            //****************************************************
            // check end residue of loop region
            //****************************************************

            if (gobbo_overlaps_with_anchor_residue_bb(
                atom,
                Conf._res[smc.ends[i_loop]],
                infos,
                ATM_CA,
                ATM_C,
                *this)) {
                return false;
            }

            //****************************************************
            // check intermediate residues of loop region
            //****************************************************

            for (size_t i_res = smc.starts[i_loop] + 1; i_res < smc.ends[i_loop]; ++i_res) {
                const Residue &other_res = Conf._res[i_res];
                // Assuming consecutive residue growth, we can skip remainder
                // of loop if anchor N atom is missing. 
                if (gobbo_is_int_res_missing(other_res)) {
                    break;
                }
                infos[e1_3_alt].res_ix = other_res._posn;
                infos[e1_3_alt].res_type = other_res._type;
                ofac2 = get_ofac2(infos[e1_3_gobbo].res_ix, infos[e1_3_alt].res_ix);
                for (int i_sim_atom = 0; i_sim_atom < NUM_BB_ATOM; ++i_sim_atom) {
                    const Atom &other_atom = other_res._atom[i_sim_atom];
                    if (VdwUtils::should_ignore_atom(other_atom)) {
                        continue;
                    }
                    // skip if 1-3 connected (connected within at most 2 bonds)
                    infos[e1_3_alt].atom_ix = i_sim_atom;
                    if (gobbo_are_1_3_connected(infos)) {
                        continue;
                    }
                    if (atoms_overlap(atom, other_atom, ofac2)) {
                        COLLISION_LOGGER_RECORD(clog, atom, other_atom, infos, "sim");
                        return false;
                    }
                } // end iteration over intermediate residue atoms
            } // end iteration over intermediate residues

        } // end iteration over simulated loops

#ifdef DISGRO_BUILD_ENABLE_SELF_CLASH_CHECKS
        //****************************************************
        // check against atoms within candidate mixed residue
        //****************************************************

        for (int i_self_atom = i_atom + 1; i_self_atom < NUM_BB_ATOM; ++i_self_atom) {
            const Atom &self_atom = res._atom[i_self_atom];
            if (VdwUtils::should_ignore_atom(self_atom)) {
                continue;
            }
            // skip if 1-3 connected (connected within at most 2 bonds)
            infos[e1_3_alt].res_ix = res._posn + bb_mixed_residue_offset[i_self_atom];
            assert(infos[e1_3_alt].res_ix >= 1);
            assert(infos[e1_3_alt].res_ix <= Conf._numRes);
            infos[e1_3_alt].atom_ix = i_self_atom;
            infos[e1_3_alt].res_type = Conf._res[infos[e1_3_alt].res_ix]._type;
            if (gobbo_are_1_3_connected(infos)) {
                continue;
            }
            if (atoms_overlap(atom, self_atom, m_ofac2[eIxOFAC_adj])) {
                COLLISION_LOGGER_RECORD(clog, atom, self_atom, infos, "self");
                return false;
            }
        } // end self check
#endif // DISGRO_BUILD_ENABLE_SELF_CLASH_CHECKS
    } // end iteration over atoms within residue

    return true;
}

/**
 * Only meant to be called during multi-loop growth -
 * restarts a loop region for regrowth. If fragment libraries
 * are enabled, then a fragment is attempted to be used
 * @param LoopChosen - index of loop to regrow
 * @param multiloops - (start, end) residue indices for simulated loop regions
 * @param min_frag_len - minimum length of loop region to allow fragment
 *  library usage (if region is smaller than this value, fragment library usage
 *  is disabled)
 * @param tail_len - The amount of residues from end of fragment to not use
 * @param smc - the simulation object
 */
void CollisionFrontend::restart_loop(const int LoopChosen,
                                     vector<loop_info> &multiloops,
                                     const vector<int> &multilooplength,
                                     const int min_frag_len,
                                     const int tail_len,
                                     SMC &smc) const {

    assert(min_frag_len >= 0);
    assert(tail_len >= 0);
    const int Original_start = multiloops[LoopChosen].Start;
    // CASE 1: We have a fragment library
    if (smc.use_multiloop_lib && (smc.end - Original_start >= min_frag_len)) {
        assert(smc.MultiLooplib[LoopChosen].size() > 0);
        // "Uniformly" select a fragment
        // @TODO - use C++11 random or boost random
        const int num_frags = static_cast<int>(smc.MultiLooplib[LoopChosen].size());
        const int first_frag_idx = rand() % num_frags;
        const int frag_end = multilooplength[LoopChosen] - tail_len;

        // Assuming libraries have uniform geometry class representation and
        // that classes are uniformly shuffled within library, then when
        // regrowing from a fragment library - can pick an arbitrary cut point
        // and proceed search for a clash-free fragment from that point

        // @TODO - consider following feature:
        // - min_frag_len: The minimum clash-free fragment length. A fragment
        //  is selected if at least the first min_frag_len residues are clash
        //  free. Note, the max fragment len is loop_len - tail_len 

        int frag_checks = 0;
        int frag_idx = first_frag_idx;
        do {

            // Cap the search span
            if (++frag_checks > smc.max_frag_checks) {
                break;
            }

            // Check if fragment collides
            if (frag_is_bb_clash_free(LoopChosen, frag_idx, frag_end, *this, smc)) {
                // No collision! Apply fragment to loop
                frag_apply(frag_idx, LoopChosen, multiloops, multilooplength, tail_len, smc);
                return;
            }

            frag_idx = (frag_idx + 1) % num_frags;

        } while (frag_idx != first_frag_idx);

        // No fragments were clash free - restart from scratch
        restart_loop_no_frag(LoopChosen, multiloops, smc);

    } // CASE 2: We don't have a fragment library
    else {
        restart_loop_no_frag(LoopChosen, multiloops, smc);
    }
}

#else // null implementation


/**
 * Does nothing
 */
void CollisionFrontend::init(const SMC &smc, const struct Params &params) {}

/**
 * @return true no matter what (does no checking)
 */
bool CollisionFrontend::gobbo_is_bb_clash_free(const Residue &res,
                                               const Structure &Conf,
                                               const SMC &smc) const { return true; }

/**
 * Only meant to be called during multi-loop growth -
 * restarts a loop region for regrowth. If fragment libraries
 * are enabled, then a fragment is attempted to be used
 * @param LoopChosen - index of loop to regrow
 * @param multiloops - (start, end) residue indices for simulated loop regions
 * @param min_frag_len - minimum length of loop region to allow fragment
 *  library usage (if region is smaller than this value, fragment library usage
 *  is disabled)
 * @param tail_len - The amount of residues from end of fragment to not use
 * @param smc - the simulation object
 */
void CollisionFrontend::restart_loop(
    const int LoopChosen,
    vector<loop_info> &multiloops,
    const vector<int> &multilooplength,
    const int min_frag_len,
    const int tail_len,
    SMC &smc) const {

    assert(min_frag_len >= 0);
    assert(tail_len >= 0);
    const int Original_start = multiloops[LoopChosen].Start;
    // CASE 1: We have a fragment library
    if (smc.use_multiloop_lib && (smc.end - Original_start >= min_frag_len)) {
        assert(smc.MultiLooplib[LoopChosen].size() > 0);
        // Randomly select a fragment
        const int select_lib_idx = rand() % smc.MultiLooplib[LoopChosen].size();
        frag_apply(select_lib_idx,
                   LoopChosen,
                   multiloops,
                   multilooplength,
                   tail_len,
                   smc);
    } // CASE 2: We don't have a fragment library
    else {
        restart_loop_no_frag(LoopChosen, multiloops, smc);
    }
}

#endif // DISGRO_BUILD_ENABLE_CLASH_FREE_BACKBONE
