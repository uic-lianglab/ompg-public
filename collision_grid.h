//****************************************************************************
// collision_grid.h
//****************************************************************************

/**
 * Collision detection is divided into broad and narrow phases. The purpose of
 * the broad phase is to perform efficient ruling out of large portions of the
 * potential collidable objects. Those objects which are not filtered via the
 * broad phase must then undergo detailed collision checking in the narrow
 * phase.
 *
 * The data structure here is a 3-D regular grid. All grid cells are of the
 * same size and are cubic.
 *
 * This specific implementation uses a hash table to hash 3D positions to a
 * set of element identifiers. The set defines which elements have bounding
 * volumes which overlap the grid cell volume.
 */

#ifndef COLLISION_GRID_H
#define COLLISION_GRID_H

//****************************************************************************
// Includes
//****************************************************************************

#include <assert.h>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>

//****************************************************************************
// Class
//****************************************************************************

/**
 * A uniform grid implementation
 */
class collision_grid {
private:

    /**
     * Assumed offsets into an extents array
     */
    enum {
        ixMinX = 0,
        ixMaxX,
        ixMinY,
        ixMaxY,
        ixMinZ,
        ixMaxZ,
        NumExtents
    };

public:

    /**
     * Unsigned integer type
     */
    typedef unsigned int uint_t;

    /**
     * Element id type
     */
    typedef uint_t elem_id_t;

    /**
     * Our set type - a set is the collection of element identifiers
     * at a single grid cell (voxel)
     */
    typedef std::unordered_set<elem_id_t> element_set_t;

    /**
     * Floating point type
     */
    typedef double real_t;

    /**
     * Default constructor
     */
    collision_grid()
            : m_inv_cell_diameter(static_cast<real_t>(0.0)), m_min_coord(static_cast<real_t>(0.0)) { }

    /**
     * Construct from cell diameter and coordinate
     * @param cell_diameter - diameter of a grid cell
     * @param min_coord - minimum corner of grid
     *  (x,y,z) = (min_coord, min_coord, min_coord)
     */
    collision_grid(const real_t cell_diameter, const real_t min_coord) {
        init(cell_diameter, min_coord);
    }

    /**
     * Initializes cubic grid structure with minimum corner at
     * (x,y,z) = (min_coord, min_coord, min_coord)
     * @param cell_diameter - diameter of a grid cell
     * @param min_coord - minimum corner of grid
     */
    void init(const real_t cell_diameter, const real_t min_coord) {
        assert(std::isfinite(cell_diameter));
        assert(cell_diameter > static_cast<real_t>(0.0));
        this->clear();
        m_inv_cell_diameter = static_cast<real_t>(1.0) / cell_diameter;
        // Add a fudge factor to make sure everything stays within
        // the grid bounding box
        assert(std::isfinite(min_coord));
        m_min_coord = min_coord -
                      (static_cast<real_t>(2.0) * std::numeric_limits<real_t>::epsilon());
    }

    /**
     * Resets to default uninitialized state
     */
    void clear() {
        this->clear_elements();
        m_inv_cell_diameter = static_cast<real_t>(0.0);
        m_min_coord = static_cast<real_t>(0.0);
    }

    /**
     * Clears elements only
     */
    void clear_elements() {
        m_grid.clear();
    }

    /**
     * Broad phase: Spatial subdivision acceleration to avoid colliding
     *  against all elements
     * @param out_results - set of element ids that need further narrow
     *  phase consideration
     * @param elem_center (x,y,z) - centroid of element
     * @param elem_radius - the radius of element
     */
    inline void filter(element_set_t &out_results,
                       const real_t elem_center_x,
                       const real_t elem_center_y,
                       const real_t elem_center_z,
                       const real_t elem_radius) const {
        // Reset results
        out_results.clear();
        this->set_reserve(out_results);

        // Determine which grid cells elemental volume intersects
        uint_t extents[NumExtents];
        this->quantize(&(extents[0]),
                       elem_center_x,
                       elem_center_y,
                       elem_center_z,
                       elem_radius);

        // See if there are any elements in the overlapping grid cells
        uint_t x, y, z;
        uint_t key;
        for (x = extents[ixMinX]; x <= extents[ixMaxX]; ++x) {
            for (y = extents[ixMinY]; y <= extents[ixMaxY]; ++y) {
                for (z = extents[ixMinZ]; z <= extents[ixMaxZ]; ++z) {
                    // @TODO - investigate moving x and y encoding out of innermost loop
                    // for possible better performance
                    key = this->encode(x, y, z);

                    grid_t::const_iterator cell_itr(m_grid.find(key));
                    if (cell_itr != m_grid.end()) {
                        this->set_union(out_results, cell_itr->second);
                    }
                }
            }
        }
    }

    /**
     * Adds a new element to the broad phase collision state. It is
     * assumed that elemental volume is less than or equal to a single
     * cell volume.
     * @param elem_id - a unique element identifier
     * @param elem_center (x,y,z) - centroid of element
     * @param elem_radius - the radius of element
     */
    inline void add(const uint_t elem_id,
                    const real_t elem_center_x,
                    const real_t elem_center_y,
                    const real_t elem_center_z,
                    const real_t elem_radius) {
        // Determine which grid cells elemental volume intersects
        uint_t extents[NumExtents];
        this->quantize(&(extents[0]),
                       elem_center_x,
                       elem_center_y,
                       elem_center_z,
                       elem_radius);

        // Add element to each overlapping grid cell
        uint_t x, y, z;
        uint_t key;
        for (x = extents[ixMinX]; x <= extents[ixMaxX]; ++x) {
            for (y = extents[ixMinY]; y <= extents[ixMaxY]; ++y) {
                for (z = extents[ixMinZ]; z <= extents[ixMaxZ]; ++z) {
                    // @TODO - investigate moving x and y encoding out of innermost loop
                    // for possible better performance
                    key = this->encode(x, y, z);
                    this->set_union(m_grid[key], elem_id);
                } // end iteration over z
            } // end iteration over y
        } // end iteration over x
    } // end function add()

private:

    /**
     * Converts a 1-d coordinate to a discretized grid cell coordinate
     * @param elem_coord - a real 1-d coordinate > m_min_coord
     * @return integral grid coordinate >= 0
     */
    inline uint_t quantize_scalar(const real_t elem_coord) const {
        assert(std::isfinite(elem_coord));
        assert(std::isfinite(m_inv_cell_diameter));
        assert(std::isfinite(m_min_coord));
        assert(elem_coord > m_min_coord);
        assert(m_inv_cell_diameter > static_cast<real_t>(0.0));
        assert(std::isfinite(((elem_coord - m_min_coord) * m_inv_cell_diameter)));
        assert(((elem_coord - m_min_coord) * m_inv_cell_diameter) >= static_cast<real_t>(0.0));
        return static_cast<uint_t>((elem_coord - m_min_coord) * m_inv_cell_diameter );
    }

    /**
     * Discretizes a 3-D real positioned cubic bounding box into the
     * integral grid cell coordinates which overlap it.
     * @param extents - the output buffer for the grid cell indices
     * @param elem_center (x,y,z) - centroid of element
     * @param elem_radius - the radius of element
     */
    inline void quantize(uint_t *const extents,
                         const real_t elem_center_x,
                         const real_t elem_center_y,
                         const real_t elem_center_z,
                         const real_t elem_radius) const {
        assert(extents != NULL);
        assert(std::isfinite(elem_radius));
        assert(elem_radius > static_cast<real_t>(0.0));
        extents[ixMinX] = this->quantize_scalar(elem_center_x - elem_radius);
        extents[ixMaxX] = this->quantize_scalar(elem_center_x + elem_radius);
        extents[ixMinY] = this->quantize_scalar(elem_center_y - elem_radius);
        extents[ixMaxY] = this->quantize_scalar(elem_center_y + elem_radius);
        extents[ixMinZ] = this->quantize_scalar(elem_center_z - elem_radius);
        extents[ixMaxZ] = this->quantize_scalar(elem_center_z + elem_radius);
    }

    /**
     * For hashing-based grid implementations
     * @param x - a quantized x coordinate
     * @param y - a quantized y coordinate
     * @param z - a quantized z coordinate
     * @return a single key representing the 3-D cell coordinate
     */
    static inline uint_t encode(const uint_t x,
                                const uint_t y,
                                const uint_t z) {
        /**
         * Bit shifts for encoding 3-D cell coordinate as a single integer
         * Based on http://http.developer.nvidia.com/GPUGems3/gpugems3_ch32.html
         */
        enum {
            // x-coord starts at bit 0
            xBitShift = 0,
            // y-coord starts at bit 10 (32-bit) or bit 21 (64-bit)
            yBitShift = ((8 * sizeof(uint_t)) / 3),
            // z-coord starts at bit 21 (32-bit) or bit 42 (64-bit)
            zBitShift = ((16 * sizeof(uint_t)) / 3)
        };
        // Assume yBitshift is median (else subsequent asserts don't make sense)
        assert(yBitShift > xBitShift);
        assert(yBitShift < zBitShift);
        // If any of these asserts trip, we need to increase the sizeof(key)
        assert(x < (1 << yBitShift));
        assert(y < (1 << yBitShift));
        assert(z < (1 << yBitShift));
        // @TODO - investing creating macros for encoding each dimension
        // separately (using & masks). This would allow triple nested loops to
        // not redo some of the shifts on each iteration.
        return ((x << xBitShift) | (y << yBitShift) | (z << zBitShift));
    }

    /**
     * Hook for hinting at a set size
     * Meant to help avoid expensive intermediate resize operations
     */
    inline static void set_reserve(element_set_t &out_set) {
#ifndef __INTEL_COMPILER // intel does not define reserve method
        enum { n_reserve_elements = 32 };
        out_set.reserve(n_reserve_elements);
#endif // __INTEL_COMPILER
    }

    /**
     * Hook for unioning two sets
     * @param out_set - the resulting union will be stored in this set
     * @param other_set - the elements in this set will be added to out_set
     */
    inline static void set_union(element_set_t &out_set, const element_set_t &other_set) {
        out_set.insert(other_set.begin(), other_set.end());
    }

    /**
     * Hook for adding a new element to an element set
     * @param out_set - the set to append the new element id to
     * @param elem_id - the element identifier
     */
    inline static void set_union(element_set_t &out_set, const elem_id_t elem_id) {
        out_set.insert(elem_id);
    }

    /**
     *  1.0 / (diameter of a grid cell)
     *  (multiplication is supposedly faster than division)
     */
    real_t m_inv_cell_diameter;

    /**
     *  The location of the min extrema corner of the grid
     *  (x, y, z) = (m_min_coord, m_min_coord, m_min_coord)
     */
    real_t m_min_coord;

    // Sparse grid structure - map from a quantized grid position
    // to the set of element identifiers at that position. A "quantized"
    // position is a 3-d position that has been mapped to a unique
    // hash key: f(real x, real y, real z) -> uint key. The key is
    // more than likely to be something along the lines:
    // f1(real x, real y, real z) -> uint grid_x, uint grid_y, uint grid_z
    // f2(grid_x, grid_y, grid_z) -> key := grid_x | grid_y | grid_z (with
    // some bitshifts to smoosh everything together). In other words,
    // we discretize a real 3-D vector into grid coordinates. Then we
    // smoosh the 3-D integer grid coordinates into a single integer key
    // which we can then use for our hash table grid representation.
    typedef std::unordered_map<unsigned int, element_set_t> grid_t;

    /**
     * Mapping from cell coordinate to element identifiers at that cell
     */
    grid_t m_grid;
};

#endif // COLLISION_GRID_H
