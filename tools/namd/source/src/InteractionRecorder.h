/**
 * Records aggregate pairwise interaction energies between residues. Utility maps
 * atom identifiers to residues.
 *
 * Only electrostatic/GBIS interaction energies supported for
 * a 0-step, single-core simulation.
 */

#ifndef _INTERACTION_RECORDER_H
#define _INTERACTION_RECORDER_H

// Master switch to toggle detailed non-bonded interaction recording
#define ENABLE_INTERACTION_RECORDER

// PDB records do not contain residue information in MEM_OPT_VERSION, hence we
// cannot map interactions back to residue -> disable recording in this case
#if (defined(ENABLE_INTERACTION_RECORDER) && !defined(MEM_OPT_VERSION))
 
#include "common.h"
 
/**
 * Singleton data structure - public interface
 */
namespace InteractionRecorder {
 
    /**
     * Used for specifying the type of interaction being recorded
     */
    enum slot {
        // Electrostatic (Coulomb potential) 
        slot_Electrostatic = 0,
        // Energy from Generalized Born Model
        slot_GBIS,
        slot_Num 
    };
 
    /**
     * Can only be called after Node::Object()->saveMolDataPointers(state)
     * as it relies on singleton Node object to have a valid pdb member
     * pointer.
     *
     * Initializes internal atomID to residueID mapping
     */
    extern void init();
    
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
    extern void inc_update(const int atomID1,
                           const int atomID2,
                           const enum slot s,
                           const BigReal value);
                           
    /**
     * Prints to stdout all recorded interactions
     */
    extern void print();
}
 
#define INT_REC_INIT InteractionRecorder::init
#define INT_REC_INC_UPDATE InteractionRecorder::inc_update
#define INT_REC_PRINT InteractionRecorder::print
 
#else
 
#define INT_REC_INIT()
#define INT_REC_INC_UPDATE(...)
#define INT_REC_PRINT()

#endif // (defined(ENABLE_INTERACTION_RECORDER) && !defined(MEM_OPT_VERSION))

#endif // _INTERACTION_RECORDER_H

