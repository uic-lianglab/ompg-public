//****************************************************************************
// collision_grid_atom.h
//****************************************************************************

/**
 * Implementation for a collision grid at atom-scale resolution
 * (as opposed to residue level resolution).
 */

#ifndef COLLISION_GRID_PROT_H
#define COLLISION_GRID_PROT_H

#include "atom.h"
#include "collision_grid.h"
#include "vdw_utils.h"

#include <assert.h>
#include <limits>

/**
 * Collision grid that uses atoms as grid elements
 */
class collision_grid_prot {
public:

    /**
     * Uniquely identifies an atom within a protein structure
     */
    struct atom_key {
        // Residue index within Structure::_res
        unsigned int res_ix;
        // Atom index within Residue::_atom
        unsigned int atom_ix;
    };

    /**
     * Structure for storing result of filter() call
     */
    typedef std::vector<atom_key> filter_results_t;

    /**
     * Default constructor
     */
    collision_grid_prot() {}

    /**
     * Initializes cubic grid structure with minimum corner at
     * (x,y,z) = (min_coord, min_coord, min_coord)
     * @param cell_diameter - diameter of a grid cell
     *  note: cell_diameter can be found using max_elem_diameter() member
     * @param min_coord - minimum corner of grid
     *  note: min_coord can be found using find_min_coord_for_grid() member
     */
    void init(const double cell_diameter, const double min_coord) {
        m_impl.init(cell_diameter, min_coord);
    }

    /**
     * Initializes cubic grid structure for parameter Structure
     * @param Conf - a protein structure
     * @param cell_diameter - diameter of a grid cell
     *  note: cell_diameter can be found using max_elem_diameter() member
     * @param min_coord - minimum corner of grid
     *  note: min_coord can be found using find_min_coord_for_grid() member
     */
    void init(const class Structure &Conf, const double cell_diameter, const double min_coord);

    /**
     * Resets to default uninitialized state
     */
    void clear() {
        m_impl.clear();
    }

    /**
     * Clears elements only
     */
    void clear_elements() {
        m_impl.clear_elements();
    }

    /**
     * Broad phase: Spatial subdivision acceleration to avoid colliding
     *  against all elements
     * @param out_results - set of element ids that need further narrow
     *  phase consideration
     * @param elem_center (x,y,z) - centroid of element
     * @param elem_radius - the radius of element
     */
    inline void filter(filter_results_t &out_results,
                       const double elem_center_x,
                       const double elem_center_y,
                       const double elem_center_z,
                       const double elem_radius) const {
        out_results.clear();
        collision_grid::element_set_t encoded_results;
        m_impl.filter(
                encoded_results,
                elem_center_x,
                elem_center_y,
                elem_center_z,
                elem_radius
        );
        // Decode results from key to (res_ix, atom_ix) pair
        if (!encoded_results.empty()) {
            out_results.resize(encoded_results.size());
            int ix = 0;
            for (collision_grid::element_set_t::const_iterator i = encoded_results.cbegin();
                 i != encoded_results.cend();
                 ++i) {
                this->decode(out_results[ix].res_ix, out_results[ix].atom_ix, *i);
                ++ix;
            }
        }
    }

    /**
     * Adds a new element to the broad phase collision state. It is
     * assumed that elemental volume is less than or equal to a single
     * cell volume.
     * @param res_ix - The index of the residue atom is within
     * @param atom_ix - The atom offset within the residue
     * @param atom - The atom coordinates
     */
    inline void add(const unsigned int res_ix,
                    const unsigned int atom_ix,
                    const Atom &atom) {
        m_impl.add(this->encode(res_ix, atom_ix),
                   atom.x,
                   atom.y,
                   atom.z,
                   this->elem_radius(atom));
    }

    
    /**
     * @param atom - atom to determine radius for
     * @return radius of parameter atom for use in grid
     */
    inline static double elem_radius(const Atom &atom) {
        return DG_ATOM_TYPE_TO_VDW_RADIUS(atom._type);
    }

private:

    /**
     * The number of bits to encode an atom index
     */
    enum {
        NUM_ATOM_IX_BITS = 4
    };

    /**
     * The number of bits to encode a residue index
     */
    enum {
        NUM_RES_IX_BITS = sizeof(unsigned int) * 8 - NUM_ATOM_IX_BITS
    };

    /**
     * Munges two integer indices into a single integer key. This works
     * because each residue only has a few atoms and the largest known human
     * protein is only tens of thousands of residues (therefore, we should
     * be able to fit all indices within the precision of an integer).
     * @param res_ix - the residue index within Structure::_res
     * @param atom_ix - the atom index within Residue::_atom
     * @return A unique integer key assuming the pair (res_ix, atom_ix) is unique
     */
    inline static unsigned int encode(const unsigned int res_ix, const unsigned int atom_ix) {
        assert((atom_ix >= 0) && (atom_ix < (0x1 << NUM_ATOM_IX_BITS)));
        assert((res_ix >= 1) && (res_ix < (0x1 << NUM_RES_IX_BITS)));
        const unsigned int elem_id = ((res_ix << NUM_ATOM_IX_BITS) | atom_ix);
        return elem_id;
    }

    /**
     * Converts from a single integer key back to a pair (res_ix, atom_ix)
     * @param res_ix - the residue index into Structure::_res
     * @param atom_ix - the atom index into Residue::_atom
     * @param elem_id - the key to decompose into (res_ix, atom_ix)
     */
    inline static void decode(unsigned int &res_ix, unsigned int &atom_ix, const unsigned int elem_id) {
        // The following only works if NUM_ATOM_IX_BITS is 4
        assert(NUM_ATOM_IX_BITS == 4);
        atom_ix = elem_id & 0xF;
        res_ix = elem_id >> NUM_ATOM_IX_BITS;
    }

    // We are implemented in terms of a collision grid
    collision_grid m_impl;
};

/**
 * Utilities for initializing a collision grid for
 * use with protein structures
 */
class collision_grid_protein_heuristics {
public:

    /**
     * Utility for finding minimum x, y, or z coordinate of a protein
     * among all atoms within protein. Needed to initialize
     * minimum extreme corner of grid bounding box.
     * @param Conf - a protein structure
     */
    static double find_min_coord_for_protein(const class Structure &Conf);

    /**
     * Heuristic: determines minimum coordinate among unsimulated atoms
     * then pads by a "worst-case" length of the simulated region
     * @param Conf - a protein structure
     * @param num_simulated_residues - the number of residues with as yet unknown positions
     * @return the minimum corner of a simulation bounding grid
     */
    static double find_min_coord_for_grid(const class Structure &Conf,
                                          const unsigned int num_simulated_residues);

    /**
     * Heuristic: determines minimum coordinate among unsimulated atoms
     * then pads by a "worst-case" length of the simulated region
     * @param Conf - a protein structure
     * @param sTart - vector of same length as eNd, contains start residue
     *  index of each simulated loop region
     * @param eNd - vector of same length as sTart, contains end residue
     *  index of each simulated loop region
     * @return the minimum corner of a simulation bounding grid
     */
    static double find_min_coord_for_grid(const class Structure &Conf,
                                          const std::vector<int> &sTart,
                                          const std::vector<int> &eNd);

    /**
     * @return diameter of max_element based on atomic level resolution
     *  (as opposed to residue level resolution)
     */
    static double max_elem_diameter();

private:
    // No instances of this class allowed!
    collision_grid_protein_heuristics();
    collision_grid_protein_heuristics(const collision_grid_protein_heuristics&);
    collision_grid_protein_heuristics& operator=(const collision_grid_protein_heuristics&);
};

#endif // COLLISION_GRID_PROT_H
