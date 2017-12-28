/**
 * Utilities for converting from internal DiSGro atom types to simple
 * C, H, O, N, or S types
 */

#ifndef VdwUtils_h
#define VdwUtils_h

//****************************************************************************
// Includes
//****************************************************************************

#include "atom.h"

#include <assert.h>

//****************************************************************************
// Class
//****************************************************************************

#define DISGRO_MAX_ATOM_TYPES MAX_ATOM_TYPE

class VdwUtils {
public:

    enum AtomFlavor {
        AF_C = 0, // carbonyl/carboxyl sp2 carbon
        AF_CA,    // aliphatic sp3 carbon
        AF_H,     // hydrogen
        AF_O,     // oxygen
        AF_N,     // nitrogen
        AF_S,     // sulfur
        AF_P,     // phosphorus
        AF_UNDEF, // undefined
        AF_NUM
    };

    /**
     * Mapping from DiSGro Atom::_type to simplified AtomFlavor (C, H, O, N, S)
     */
    static const AtomFlavor disgro_atom_type_to_CHONS[DISGRO_MAX_ATOM_TYPES];

    /**
     * Mapping from AtomFlavor to van der Waals radius
     */
    static const double atom_flavor_to_vdw_radius[AF_NUM];

    /**
     * @return 1 (TRUE) if atom should be ignored by all grid methods, 0 (FALSE) o/w
     */
    inline static unsigned int should_ignore_atom(const Atom &atom) {
        // Ignore undefined or hydrogen atoms
        // Ignore uninitialized atoms which have been placed at origin
        assert(atom.is_at_origin() || (atom._type > 0));
        assert(atom.is_at_origin() || (atom._type < DISGRO_MAX_ATOM_TYPES));
        if (atom.is_at_origin() ||
            disgro_atom_type_to_CHONS[atom._type] == VdwUtils::AF_H ||
            disgro_atom_type_to_CHONS[atom._type] == VdwUtils::AF_UNDEF) {
            return 1;
        }
        return 0;
    }

    /**
     * @return Max possible atom radius based on atom flavors
     */
    static double max_atom_radius() {
        const unsigned int n = sizeof(atom_flavor_to_vdw_radius) / sizeof(double);
        double max_radius = VdwUtils::atom_flavor_to_vdw_radius[0];
        for (unsigned int i = 1; i < n; ++i) {
            if (atom_flavor_to_vdw_radius[i] > max_radius) {
                max_radius = atom_flavor_to_vdw_radius[i];
            }
        }
        return max_radius;
    }

    /**
     * @return Max possible atom diameter based on atom flavors
     */
    static double max_atom_diameter() {
        return 2.0 * max_atom_radius();
    }

private:
    // No instances of this class allowed!
    VdwUtils();
    VdwUtils(const VdwUtils&);
    VdwUtils& operator=(const VdwUtils&);
};

//****************************************************************************
// Macros
//****************************************************************************

// Do not use this macro directly
// Use DG_ATOM_TYPE_TO_CHONS instead as it provides
// a checked version if asserts are enabled. This macro is
// defined to adhere to DRY (do-not-repeat-yourself) idiom
#define DO_DG_ATOM_TYPE_TO_CHONS(dg_atom_type) \
    VdwUtils::disgro_atom_type_to_CHONS[dg_atom_type]

// Do not use this macro directly
// Use DG_ATOM_TYPE_TO_ES_VDW_RADIUS instead as it provides
// a checked version if asserts are enabled.
#define DO_DG_ATOM_TYPE_TO_VDW_RADIUS(dg_atom_type) \
    VdwUtils::atom_flavor_to_vdw_radius[DG_ATOM_TYPE_TO_CHONS(dg_atom_type)]

// Provide checked version of look-up macro if asserts enabled
#ifdef NDEBUG
#   define DG_ATOM_TYPE_TO_CHONS(dg_atom_type) \
        DO_DG_ATOM_TYPE_TO_CHONS(dg_atom_type)
#   define DG_ATOM_TYPE_TO_VDW_RADIUS(dg_atom_type) \
        DO_DG_ATOM_TYPE_TO_VDW_RADIUS(dg_atom_type)
#else

/**
 * Converts a DiSGro atom type (Atom::_type) to the simplified
 * AtomFlavor enumeration (C, H, O, N, S, undefined). This simplified
 * type can be used to set a common Van der Waals radius.
 * @param dg_atom_type - the DiSGro atom type
 * @return the enumerated AtomFlavor mapping to parameter type
 */
inline VdwUtils::AtomFlavor DG_ATOM_TYPE_TO_CHONS(const short dg_atom_type) {
    assert(dg_atom_type > 0);
    assert(dg_atom_type < DISGRO_MAX_ATOM_TYPES);
    return DO_DG_ATOM_TYPE_TO_CHONS(dg_atom_type);
}

/**
 * @return Van der Waal radius for parameter DiSGro atom type (Atom::_type).
 * This radius is what was used during training of the electrostatic
 * density model and it may differ from the radius defined in atomProp2.txt
 */
inline double DG_ATOM_TYPE_TO_VDW_RADIUS(const short dg_atom_type) {
    assert(dg_atom_type > 0);
    assert(dg_atom_type < DISGRO_MAX_ATOM_TYPES);
    assert(VdwUtils::disgro_atom_type_to_CHONS[dg_atom_type] >= 0);
    assert(VdwUtils::disgro_atom_type_to_CHONS[dg_atom_type] < VdwUtils::AF_NUM);
    return DO_DG_ATOM_TYPE_TO_VDW_RADIUS(dg_atom_type);
}

#endif // NDEBUG

#endif // VdwUtils_h
