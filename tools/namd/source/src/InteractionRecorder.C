/**
 * Records aggregate pairwise interaction energies between residues. Utility maps
 * atom identifiers to residues.
 *
 * Only electrostatic/GBIS interaction energies supported for
 * a 0-step, single-core simulation.
 */
 
#include "InteractionRecorder.h"

// PDB atom records not compatible with MEM_OPT_VERSION
#if (defined(ENABLE_INTERACTION_RECORDER) && !defined(MEM_OPT_VERSION))

#include "Node.h"
#include "PDB.h"

#include <algorithm>
#include <vector>

/**
 * @TODO - Everything is printf, would be nice to use actual
 *  NAMD logging system
 */
 
// Store singleton data in unnamed namespace
namespace
{

    /**
     * Cached pointer to atom and residue information
     */
    static PDB *pdb = NULL;
    
    /**
     * 2-d matrix for each slot (energy type) being recorded. There is no
     * attempt at memory efficiency - we really only need a symmetric,
     * triangular matrix but we are allocating the full matrix. Also, residue
     * indices are 1-based to better match PDB records which means there are
     * dummy records inserted just for padding (this padding may be especially
     * egregious if the first PDB residue identifier does not actually start
     * at 1 but rather at some offset >> 1.) It assumed that PDB residue
     * identifiers will never be less than 1.
     */
    static std::vector<std::vector<BigReal> > records[InteractionRecorder::slot_Num];
    
    /**
     * 2-d matrix slot which marks which interactions have actually been
     * recorded. This is to minimize the log spew in the print() method as
     * only marked interactions will be output. No attempt has been made at
     * memory efficiency (we only need a triangle matrix but are allocating
     * a full matrix with dummy entries for padding) as this makes things
     * simpler for now.
     */
    static std::vector<std::vector<unsigned int> > marks;
    
    /**
     * Utility to allocate a 'size' number of vector of vectors, where each
     * internal vector also has 'size' number of elements. In other words,
     * just like a square matrix.
     */
    template <typename t>
    static void init_square_matrix(std::vector<std::vector< t > > &m,
                                   const int size,
                                   const t default_value)
    {
        const std::vector< t > def_vector(size, default_value);
        m.clear();
        m.resize(size, def_vector);
    }  
}
 
/**
 * Can only be called after Node::Object()->saveMolDataPointers(state)
 * as it relies on singleton Node object to have a valid pdb member
 * pointer.
 *
 * Initializes internal atomID to residueID mapping
 */
void InteractionRecorder::init()
{
    // PDB structure contains the atom to residue mapping.
    // Simply cache the pointer
    pdb = Node::Object()->pdb;
    
    // Determine max residue identifier and allocate storage accordingly
    
    const int num_atoms = pdb->num_atoms();
    if (num_atoms < 0)
    {
        printf("INTERACTION LOGGER: Warning, no atoms found during initialization.\n");
        return;
    }
    
    PDBAtom *last_atom = pdb->atom_unsafe(num_atoms-1);
    const int max_residue_id = last_atom->residueseq();
    
    for (int itr_slot=0; itr_slot<InteractionRecorder::slot_Num; ++itr_slot) {
        // Make 1-based indexing by adding 1 element padding
        init_square_matrix(records[itr_slot], max_residue_id + 1, 0.0);
    }
    
    init_square_matrix(marks, max_residue_id + 1, (unsigned int)0);
    
    printf("INTERACTION LOGGER: Initialized recording for %d x %d residues.\n",
           max_residue_id, max_residue_id);
}

/**
 * Can only be called after init() has been called as assumes atomID to
 * residueID mapping is valid.
 *
 * Incrementally updates the interaction energy for the residue pair
 * denoted by parameter atom identifiers. If interacting pair has not
 * been encountered, then the pair's interaction energy is set to
 * parameter value.
 *
 * @param atomID1 - first atom of interaction pair. A 0-based atom index,
 *  assumed to be in same order as listing of atoms within PDB file.
 * @param atomID2 - second atom of interaction pair. A 0-based atom index,
 *  assumed to be in same order as listing of atoms within PDB file.
 * @param slot - specifies which interaction type is being logged
 * @param value - the value to add to current interaction energy value.
 */
void InteractionRecorder::inc_update(const int atomID1,
                                     const int atomID2,
                                     const enum slot s,
                                     const BigReal value)
{
    const int resID1 = pdb->atom_unsafe(atomID1)->residueseq();
    const int resID2 = pdb->atom_unsafe(atomID2)->residueseq();
    const int min_res_id = std::min(resID1, resID2);
    const int max_res_id = std::max(resID1, resID2);
    marks[min_res_id][max_res_id] = 1;
    records[s][min_res_id][max_res_id] += value;
}
                           
/**
 * Prints to stdout all recorded interactions
 */
void InteractionRecorder::print() 
{
    const unsigned int max_res = (unsigned int)marks.size();
    printf("BEGIN_INTERACTION_RECORDS\n");
    printf("RES1,RES2,Electrostatic,GBIS\n");
    for (unsigned int i=0; i<max_res; ++i) {
        for (unsigned int j=i; j<max_res; ++j) {
            if (marks[i][j]) {
                const BigReal es = records[slot_Electrostatic][i][j];
                const BigReal gb = records[slot_GBIS][i][j];
                printf("%u,%u,%f,%f\n", i, j, (es + gb), gb);
            }
        }
    }
    printf("END_INTERACTION_RECORDS\n");
}

#endif // (defined(ENABLE_INTERACTION_RECORDER) && !defined(MEM_OPT_VERSION))

