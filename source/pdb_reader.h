#ifndef PDB_READER_H
#define PDB_READER_H

#include "point.h"
#include "util.h"

#include <assert.h>
#include <string>

/**
 * Enumerated tokens within a PDB ATOM record
 */
enum EPdbAtomRecToken {
    // Record name
    PdbRecName = 0,
    // Atom serial number
    PdbAtomSerialNo,
    // Atom name
    PdbAtomName,
    // Alternate location
    PdbAtomAltLoc,
    // Residue name
    PdbAtomResName,
    // Chain identifier
    PdbAtomChain,
    // Residue sequence number
    PdbAtomResSeqNo,
    // Residue insertion code
    PdbAtomResInsCode,
    // Atom x-coordinate
    PdbAtomXCoord,
    // Atom y-coordinate
    PdbAtomYCoord,
    // Atom z-coordinate
    PdbAtomZCoord,
    // Atom occupancy
    PdbAtomOcc,
    // Temperature factor
    PdbAtomTempFact,
    // Segment identifier
    PdbAtomSeg,
    // Element symbol
    PdbAtomElem,
    // Charge
    PdbAtomCharge,
    // Number of tokens in a PDB ATOM record
    PdbAtomRec_NUM
};

/**
 * Stores offset and length for a token within an PDB ATOM record
 */
struct PdbAtomRecTokenInfo {
    // Start index of token
    unsigned int start;
    // Number of characters defining token
    unsigned int len;
};

/**
 * Effectively a namespace for PDB parsing utilities
 */
class Pdb {
public:

    /**
     * Obtains a token from a string containing a PDB ATOM record.
     * The token is trimmed for leading and trailing whitespace
     * @param token - stores value of trimmed token
     * @param atom_rec_line - line from PDB file defining ATOM record
     * @param id - specifies which token to extract
     */
    static void get_atom_rec_token(
        std::string& token,
        const std::string &atom_rec_line,
        const EPdbAtomRecToken id) {
        assert(PdbAtomRec_NUM > id);
        assert(PdbAtomRec_NUM == (sizeof(AtomRecTokenInfos) / sizeof(struct PdbAtomRecTokenInfo)));
        const struct PdbAtomRecTokenInfo& nfo = AtomRecTokenInfos[id];
        if (nfo.start >= atom_rec_line.size()) {
            // Return empty if line is too short
            token.clear();
        }
        else {
            // Else parse line for token
            token = atom_rec_line.substr(nfo.start, nfo.len);
            // Trim ends of whitespace
            trim(token);
        }
    }

    /**
     * @return true if line represents an PDB ATOM record, false o/w
     */
    static bool is_atom_rec(const std::string &line) {
        std::string token;
        get_atom_rec_token(token, line, PdbRecName);
        return token == "ATOM";
    }

    /**
     * Parses x,y,z coordinates from PDB ATOM record
     */
    static void parse_atom_coords(Point &atom, const std::string &atom_rec_line) {
        std::string token;
        // Parse atom coordinates
        get_atom_rec_token(token, atom_rec_line, PdbAtomXCoord);
        atom.x = atof(token.c_str());
        get_atom_rec_token(token, atom_rec_line, PdbAtomYCoord);
        atom.y = atof(token.c_str());
        get_atom_rec_token(token, atom_rec_line, PdbAtomZCoord);
        atom.z = atof(token.c_str());
    }

private:

    /**
     * Table storing PDB ATOM record token character offsets
     */
    static const PdbAtomRecTokenInfo AtomRecTokenInfos[PdbAtomRec_NUM];
};

#endif // PDB_READER_H
