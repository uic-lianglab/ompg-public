#include "membrane.h"
#include "util.h"
#include "atom.h"
#include "pdb_reader.h"

#include <assert.h>
#include <fstream>
#include <iostream>

namespace {

    /**
     * @return true if line represent a PDB ATOM record from a membrane segment, false o/w
     */
    bool is_pdb_membrane_atom_rec(const std::string &line) {
        if (Pdb::is_atom_rec(line)) {
            std::string token;
            Pdb::get_atom_rec_token(token, line, PdbAtomSeg);
            return token == "MEMB";
        }
        return false;
    }

    /**
     * Map PDB atom name to an atom flavor
     */
    VdwUtils::AtomFlavor parse_pdb_atom_flavor(const std::string &atom_rec_line) {
        // Determine atom flavor
        std::string token;
        Pdb::get_atom_rec_token(token, atom_rec_line, PdbAtomName);
        if (token.empty()) {
            std::cout << "Warning: missing atom name for line:\n\t"
                      << atom_rec_line << std::endl;
            return VdwUtils::AF_UNDEF;
        }
        switch (token[0]) {
        case 'C': return VdwUtils::AF_C;
        case 'H': return VdwUtils::AF_H;
        case 'O': return VdwUtils::AF_N;
        case 'N': return VdwUtils::AF_N;
        case 'S': return VdwUtils::AF_S;
        case 'P': return VdwUtils::AF_P;
        default:
            std::cout << "Warning: unrecognized membrane atom type: " << token << std::endl
                      << "\tfor line: " << atom_rec_line << std::endl;
            return VdwUtils::AF_UNDEF;
        };
    }

    /**
     * Utility for finding minimum x, y, or z coordinate of a membrane
     * among all atoms within membrane. Needed to initialize
     * minimum extreme corner of grid bounding box.
     * @param atoms - list of membrane atoms
     */
    double find_min_coord_for_membrane(const Membrane::atom_list_t &atoms) {
        double min_coord = std::numeric_limits<double>::max();
        for (Membrane::atom_list_t::const_iterator atom_it = atoms.cbegin();
             atom_it != atoms.cend(); ++atom_it) {

            min_coord = std::min(std::min(std::min(min_coord, atom_it->x), atom_it->y), atom_it->z);
        }
        // Pad min coord by largest element diameter
        min_coord -= VdwUtils::max_atom_diameter();
        return min_coord;
    }

    /**
     * If membrane atoms exist, will initialize a broad phase collision grid
     */
    void conditional_init_membrane_grid(collision_grid &grid, Membrane::atom_list_t &atoms) {
        // Early out if no atoms
        if (atoms.empty()) {
            return;
        }

        // Initialize grid dimensions
        grid.init(VdwUtils::max_atom_diameter(), find_min_coord_for_membrane(atoms));

        // Add membrane atoms to grid
        const collision_grid::uint_t num_atoms =
            static_cast<collision_grid::uint_t>(atoms.size());
        for (collision_grid::uint_t i=0; i<num_atoms; ++i) {
            assert(atoms[i].flavor != VdwUtils::AF_UNDEF);
            assert(atoms[i].flavor < VdwUtils::AF_NUM);
            grid.add(
                i,
                atoms[i].x,
                atoms[i].y,
                atoms[i].z,
                VdwUtils::atom_flavor_to_vdw_radius[atoms[i].flavor]);
        }
    }

    /**
     * @return 1 if atoms overlap, 0 otherwise
     */
    inline int atoms_overlap(const Atom &a, const Membrane::atom_t &b, const double ofac2) {
        const double a_r_vdw = DG_ATOM_TYPE_TO_VDW_RADIUS(a._type);
        const double b_r_vdw = VdwUtils::atom_flavor_to_vdw_radius[b.flavor];
        const double disquare = a.disquare(b);
        double sum_radii_square = a_r_vdw + b_r_vdw;
        sum_radii_square *= sum_radii_square;
        return (disquare < ofac2*sum_radii_square);
    }
} // end of anonymous namespace

/**
 * Initializes membrane from PDB file. Will wipe internal state before
 * initializing even if parameters are invalid.
 * @param pdb_path - path to PDB file containing membrane data
 */
void Membrane::init(const std::string& pdb_path) {
    // Clear internal state
    m_atoms.clear();
    m_grid.clear();

    // Early out if no PDB specified
    if (pdb_path.empty()) {
        return;
    }

    // Early out if PDB path is invalid
    if (!FileExists(pdb_path)) {
        std::cout << "Warning: membrane PDB file not found: " << pdb_path << std::endl;
        return;
    }

    // Parse PDB file line by line
    atom_t atom;
    std::ifstream fpdb(pdb_path);
    for (std::string line; std::getline(fpdb, line);) {

        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Skip non-atom records and non-membrane atoms
        if (!is_pdb_membrane_atom_rec(line)) {
            continue;
        }

        // Parse atom flavor
        atom.flavor = parse_pdb_atom_flavor(line);

        // Skip if unrecognized atom
        if (atom.flavor == VdwUtils::AF_UNDEF) {
            continue;
        }

        // Parse atom coordinates
        Pdb::parse_atom_coords(atom, line);

        // Add atom to membrane
        m_atoms.push_back(atom);

    } // end iteration over lines of PDB file

    // Close file handle
    fpdb.close();

    // Add membrane atoms to bounding grid
    conditional_init_membrane_grid(m_grid, m_atoms);

    std::cout << "Initialized lipid membrane from " 
              << pdb_path << "\n\twith " << m_atoms.size() << " atoms.\n";
}

/**
 * @param atom - a protein atom to check collisions against membrane
 * @param ofac2 - the squared overlap factor, if ratio of squared distance
 *  to a membrane atom and the squared sum of vdw radii is less than this
 *  value, then the atoms are colliding
 * @return 1 if param atom collides with membrane atoms, 0 o/w
 */
int Membrane::collides(const Atom &atom, const double ofac2) const {
    // Early out if no membrane loaded
    if (m_atoms.empty()) {
        return 0;
    }

    // Broad phase
    collision_grid::element_set_t nearby_atoms;
    m_grid.filter(nearby_atoms,
        atom.x, atom.y, atom.z,
        DG_ATOM_TYPE_TO_VDW_RADIUS(atom._type));

    // Narrow phase
    for (collision_grid::element_set_t::const_iterator i = nearby_atoms.cbegin();
         i != nearby_atoms.cend(); ++i) {
        // Bounds check
        assert(*i >= 0);
        assert(*i < static_cast<collision_grid::uint_t>(this->m_atoms.size()));
        const atom_t &nearby_atom = this->m_atoms[*i];
        assert(nearby_atom.flavor != VdwUtils::AF_UNDEF);

        if (atoms_overlap(atom, nearby_atom, ofac2)) {
            // Collision detected
            return 1;
        }
    }

    // No collisions detected
    return 0;
}
