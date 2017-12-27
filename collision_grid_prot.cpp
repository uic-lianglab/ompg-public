//
// Created by apr on 10/13/15.
//

#include "collision_grid_prot.h"
#include "structure.h"

#include <algorithm>

/**
 * Initializes cubic grid structure for parameter Structure
 * @param Conf - a protein structure
 * @param cell_diameter - diameter of a grid cell
 *  note: cell_diameter can be found using max_elem_diameter() member
 * @param min_coord - minimum corner of grid
 *  note: min_coord can be found using find_min_coord_for_grid() member
 */
void collision_grid_prot::init(const Structure &Conf, const double cell_diameter, const double min_coord) {
    // Reset state
    this->clear();
    // Initialize grid implementation
    m_impl.init(cell_diameter, min_coord);
    // Add atoms to grid
    // Residue indices are 1-based
    const unsigned int numRes = static_cast<unsigned int>(Conf._numRes);
    for (unsigned int res_ix = 1; res_ix <= numRes; ++res_ix) {
        const unsigned int numAtom = static_cast<unsigned int>(Conf._res[res_ix]._numAtom);
        for (unsigned int atom_ix = 0; atom_ix < numAtom; ++atom_ix) {
            const Atom &atom = Conf._res[res_ix]._atom[atom_ix];
            if (VdwUtils::should_ignore_atom(atom)) {
                continue;
            }
            const unsigned int elem_id = this->encode(res_ix, atom_ix);
            m_impl.add(
                    elem_id,
                    atom.x,
                    atom.y,
                    atom.z,
                    this->elem_radius(atom)
            );
        } // end iteration over atoms
    } // end iteration over residues
}

/**
 * Utility for finding minimum x, y, or z coordinate of a protein
 * among all atoms within protein. Needed to initialize
 * minimum extreme corner of grid bounding box.
 * @param Conf - a protein structure
 */
double collision_grid_protein_heuristics::find_min_coord_for_protein(const Structure &Conf) {
    double min_coord = std::numeric_limits<double>::max();
    // Residue indices are 1-based
    for (int res_ix = 1; res_ix <= Conf._numRes; ++res_ix) {
        for (int atom_ix = 0; atom_ix < Conf._res[res_ix]._numAtom; ++atom_ix) {
            const Atom &atom = Conf._res[res_ix]._atom[atom_ix];
            if (VdwUtils::should_ignore_atom(atom)) {
                continue;
            }
            min_coord = std::min(std::min(std::min(min_coord, atom.x), atom.y), atom.z);
        } // end iteration over atoms
    } // end iteration over residues

    // Pad min coord by largest element diameter
    min_coord -= max_elem_diameter();
    return min_coord;
}

/**
 * Heuristic: determines minimum coordinate among unsimulated atoms
 * then pads by a "worst-case" length of the simulated region
 * @param Conf - a protein structure
 * @param num_simulated_residues - the number of residues with as yet unknown positions
 * @return the minimum corner of a simulation bounding grid
 */
double collision_grid_protein_heuristics::find_min_coord_for_grid(
    const Structure &Conf,
    const unsigned int num_simulated_residues) {

    assert(num_simulated_residues > 0);
    double min_coord = find_min_coord_for_protein(Conf);
    // @HACK - pad the min_coord by a fudge factor that we expect simulated loop regions
    // to be contained within
    min_coord -= ((num_simulated_residues * Residue::size[ARG]) + CUB_SIZE);
    return min_coord;
}

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
double collision_grid_protein_heuristics::find_min_coord_for_grid(
    const Structure &Conf,
    const std::vector<int> &sTart,
    const std::vector<int> &eNd) {

    unsigned int num_simulated_residues = 0;
    assert(sTart.size() == eNd.size());
    const size_t num_loops = sTart.size();
    for (size_t i=0; i<num_loops; ++i) {
        assert(eNd[i] >= sTart[i] );
        num_simulated_residues += eNd[i] - sTart[i] + 1;
    }
    return find_min_coord_for_grid(Conf, num_simulated_residues);
}

/**
* @return diameter of max_element based on atomic level resolution
*  (as opposed to residue level resolution)
*/
double collision_grid_protein_heuristics::max_elem_diameter() {
    return VdwUtils::max_atom_diameter();
}
