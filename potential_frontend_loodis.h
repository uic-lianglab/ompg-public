#ifndef PotentialFrontendLoodis_h
#define PotentialFrontendLoodis_h

/** 
 * Interface for computing potential energy calculations. Allows easier hooking
 * if we need to modify/add potentials. This implementation only uses the
 * LOODIS statistical potential. This class used to be called
 * CompositePotentialCalculator (see history in older repositories)
 */
class PotentialFrontend {
public:

    // Must be called once a candidate residue has chosen:
    // MUST BE CALLED BEFORE CHOSEN RESIDUE IS ADDED TO SAMPLE
    // @param Conf - the protein structure containing candidate residue
    // @param chosen_candidate_res - the candidate residue selected to be added to Conf
    // @param chosen_candidate_key - The index identifying the chosen candidate
    // @param sTart - vector of same length as eNd, contains start residue
    //  index of each simulated loop region
    // @param eNd - vector of same length as sTart, contains end residue
    //  index of each simulated loop region
    // @param target_res_ix - index within Conf at which chosen_candidate_res is to be added to the sample
    // @param LoopChosen - the index of the loop that residue is being added to
    void on_grow_one_res_chosen(const Structure &Conf,
                                const Residue &chosen_candidate_res,
                                const int chosen_candidate_key,
                                const std::vector<int> &sTart,
                                const std::vector<int> &eNd,
                                const int target_res_ix,
                                const int LoopChosen) { /*do nothing!*/ }

    // Must be called once a candidate residue has chosen:
    // MUST BE CALLED BEFORE CHOSEN RESIDUE IS ADDED TO SAMPLE
    // @param Conf - the protein structure containing candidate residue
    // @param chosen_candidate_res - the candidate residue selected to be added to Conf
    // @param chosen_candidate_key - The index identifying the chosen candidate
    // @param sTart - vector of same length as eNd, contains start residue
    //  index of each simulated loop region
    // @param eNd - vector of same length as sTart, contains end residue
    //  index of each simulated loop region
    // @param target_res_ix - index within Conf at which chosen_candidate_res is to be added to the sample
    // @param LoopChosen - the index of the loop that residue is being added to - necessary
    //  for keeping track of fragment density contributions.
    void on_grow_sc_res_chosen(const Structure &Conf,
                               const Residue &chosen_candidate_res,
                               const int chosen_candidate_key,
                               const std::vector<int> &sTart,
                               const std::vector<int> &eNd,
                               const int target_res_ix,
                               const int LoopChosen) { /*do nothing*/ }

    // Must be called after a sample is completed
    void on_sample_finish() { /*do nothing!*/ }

    // MUST BE CALLED AFTER THE FRAGMENT IS ADDED TO SAMPLE
    // Must be called after a loop is restarted from a fragment library.
    // @param Conf - a whole protein structure
    // @param LoopChosen - index of loop to be replaced
    // @param loop_start_ix - the start residue index of the loop region within Conf
    // @param loop_end_ix - the end residue index of the entire loop region within Conf
    // @param frag_start_ix - the start index of the fragment region to regrow from
    // @param frag_end_ix - the end index of the fragment region to regrow from (must be less
    //  than fragment._numRes as terminal residue is special case (CA and C atoms known) which
    //  is currently not being handled. (see on_analytic_closure() for how to handle this
    //  case should need arise)
    // @param core - contains cached indices and types for ionizable residues/atoms
    void on_restart_from_fragment(const Structure &Conf,
                                  const int LoopChosen,
                                  const int loop_start_ix,
                                  const int loop_end_ix,
                                  const int frag_start_ix,
                                  const int frag_end_ix) { /*do nothing!*/ }

    // Callback for when a loop is successfully closed.
    // @param Conf - a protein structure - already modified with closed atoms
    // @param LoopChosen - the loop index which has just been closed
    // @param loop_end_ix - the index into Conf of the loop end/terminal residue
    void on_analytic_closure(const Structure &Conf,
                             const int LoopChosen,
                             const int loop_end_ix) { /*do nothing!*/ }

    // @param Conf - the protein structure
    // @param res - the candidate residue
    // @param List - indices of residues to check for in statistical potential
    // @param List_size - the size of (number of elements of) param List
    // @param res_ix - the index of param Res within Conf
    // @param candidate_key - the index of the candidate being grown
    // @param sTart - vector of same length as eNd, contains start residue
    //  index of each simulated loop region
    // @param eNd - vector of same length as sTart, contains end residue
    //  index of each simulated loop region
    // @param b_is_hybrid_res - When calling from grow_one(), residue is actually mix of atoms
    //  from res_ix and res_ix + 1
    // @param LoopChosen - the loop selected to be grown
    // @return potential energy
    double one_res_en_all_list(const Structure &Conf,
                               const Residue &res,
            // loodis parameters:
                               const std::vector<int> &List,
            // (possible) electrostatic parameters:
                               const int res_ix,
                               const int candidate_key,
                               const std::vector<int> &sTart,
                               const std::vector<int> &eNd,
                               const bool b_is_hybrid_res,
                               const int LoopChosen) {
        // Statistical energy calculation
        extern double one_res_en_loodis_all_list(const Structure &Conf, const Residue &res,
                                                 const std::vector<int> &List);
        return one_res_en_loodis_all_list(Conf, res, List);
    }

    // @param Conf - the protein structure
    // @param res - the candidate residue
    // @param res_ix - the index of param Res within Conf
    // @param candidate_key - the index of the candidate being grown
    // @param sTart - vector of same length as eNd, contains start residue
    //  index of each simulated loop region
    // @param eNd - vector of same length as sTart, contains end residue
    //  index of each simulated loop region
    // @param LoopChosen - the loop selected to be grown
    // @return potential energy
    double one_res_en_sc(const Structure &Conf,
                         const Residue &res,
                         const int res_ix,
                         const int candidate_key,
                         const std::vector<int> &sTart,
                         const std::vector<int> &eNd,
                         const int LoopChosen) {

        // Statistical energy function
        extern double one_res_en_loodis_sc(const Structure &conf, const Residue &res, const int position,
                                           const int Start, const int End);
        // Compute statistical potential energy
        return one_res_en_loodis_sc(Conf, res, res_ix, sTart[LoopChosen], eNd[LoopChosen]);
    }

    // Computes potential energy for entire structure
    // @param Conf - the protein structure to compute potential energy of
    // @param sTart - vector of same length as eNd, contains start residue
    //  index of each simulated loop region
    // @param eNd - vector of same length as sTart, contains end residue
    //  index of each simulated loop region
    double energy(Structure &Conf,
                  const vector<int> &sTart,
                  const vector<int> &eNd) {
        // Statistical energy function
        extern double loodis_e(const Structure &Conf, const vector<int> &sTart, const vector<int> &eNd, bool type);
        memset(&Conf._enStats, 0, sizeof(EnergyStats));
        // Compute statistical potential energy
        Conf._enStats.loodis_term = loodis_e(Conf, sTart, eNd, 1 /*type = 1 -> use side chains*/);
        return Conf._enStats.loodis_term;
    }
};

#endif // PotentialFrontendLoodis_h
