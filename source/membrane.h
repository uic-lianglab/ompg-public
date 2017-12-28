#ifndef MEMBRANE_H
#define MEMBRANE_H

#include "point.h"
#include "vdw_utils.h"
#include "collision_grid.h"

#include <string>
#include <vector>

/**
 * Membrane (lipid bilayer) representation
 * Utilities for loading and interacting with a membrane
 */
class Membrane {
public:

    /**
     * A membrane atom
     */
    struct atom_t : public Point {
        /**
         * The simplified atom type
         */
        VdwUtils::AtomFlavor flavor;
    };

    /**
     * The data structure for storing membrane atoms
     */
    typedef std::vector<atom_t> atom_list_t;

    /**
     * Initializes membrane from PDB file. Will wipe internal state before
     * initializing even if parameters are invalid.
     * @param pdb_path - path to PDB file containing membrane data
     */
    void init(const std::string& pdb_path);

    /**
     * @param atom - a protein atom to check collisions against membrane
     * @param ofac2 - the squared overlap factor, if ratio of squared distance
     *  to a membrane atom and the squared sum of vdw radii is less than this
     *  value, then the atoms are colliding
     * @return 1 if param atom collides with membrane atoms, 0 o/w
     */
    int collides(const class Atom &atom, const double ofac2) const;

    /**
     * Collection of atoms defining the membrane
     */
    atom_list_t m_atoms;

    /**
     * Volumetric grid for broad phase collisions. Grid cells store
     * indices into atoms array.
     */
    collision_grid m_grid;
};

#endif // MEMBRANE_H
