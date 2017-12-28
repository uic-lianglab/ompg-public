#include "pdb_reader.h"

/**
 * Table storing PDB ATOM record token character offsets
 */
const PdbAtomRecTokenInfo Pdb::AtomRecTokenInfos[PdbAtomRec_NUM] = {
    // Record name
    { 0, 6 },
    // Atom serial number
    { 6, (11 - 6) },
    // Atom name
    { 12, (16 - 12) },
    // Alternate location
    { 16, 1 },
    // Residue name
    { 17, (20 - 17) },
    // Chain identifier
    { 21, 1 },
    // Residue sequence number
    { 22, (26 - 22) },
    // Residue insertion code
    { 26, 1 },
    // Atom x-coordinate
    { 30, (38 - 30) },
    // Atom y-coordinate
    { 38, (46 - 38) },
    // Atom z-coordinate
    { 46, (54 - 46) },
    // Atom occupancy
    { 54, (60 - 54) },
    // Temperature factor
    { 60, (66 - 60) },
    // Segment identifier
    { 72, (76 - 72) },
    // Element symbol
    { 76, (78 - 76) },
    // Charge
    { 78, (80 - 78) }
};
