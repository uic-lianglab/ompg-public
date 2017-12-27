//****************************************************************************
// mutation_manager.h
//
// Class for reading and applying mutation options from .ini
//****************************************************************************

//****************************************************************************
// Includes
//****************************************************************************

#include "mutation_manager.h"
#include "structure.h"
#include "util.h"

#include <iostream>
#include <ctype.h>

// Loads the mutations from .ini file
// @param ini_fpath - string with path to .ini file
// @return TRUE if mutations read properly, FALSE otherwise
void MutationManager::load(const std::string &ini_fpath) {
    std::cout << "Reading mutation file: " << ini_fpath << std::endl;
    this->m_cfg_file.load(ini_fpath);
}

// Modifies parameter structure with configured mutations.
// Note: all atoms of mutated residues are moved to origin.
// @param Conf - protein structure to apply mutations to
void MutationManager::apply(Structure &Conf) const {

    for (auto i = m_cfg_file.getContents().cbegin();
         i != m_cfg_file.getContents().cend();
         ++i) {

        const int mut_idx = stoi(i->first);

        Residue *pRes = NULL;
        for (int j=1; j<=Conf._numRes; ++j) {
            if (Conf._res[j]._pdbIndex == mut_idx) {
                pRes = &(Conf._res[j]);
                break;
            }
        }

        if (pRes == NULL) {
            std::cout << "ERROR: mutation residue position not found at: " << mut_idx << std::endl;
            exit(0);
        }

        // Target amino acid type must be 1-letter AA code
        if (Residue::AIMap.find(i->second) == Residue::AIMap.end()) {
            std::cout << "ERROR: mutation target amino acid has unrecognized type: " << i->second << std::endl;
            exit(0);
        }

        Residue &res = *pRes;
        std::cout << "Mutating residue " << mut_idx
                  << " (" << res._posn << ")"
                  << " from " << res.get_name_1() << " to " << i->second << std::endl;

        res._type = Residue::AIMap.find(i->second)->second;
        res._numAtom = Residue::numAtom[res._type];

        // Move mutated residue atoms to origin
        for (int k = 0; k < res._numAtom; ++k) {
            Atom &atom = res._atom[k];
            atom._posn = k;
            atom._name = Residue::cType[res._type][k];
            atom._type = Residue::vdwType[res._type][k];
            atom.move_to_origin();
        }
    }

    // @TODO - is this necessary?
    if (!m_cfg_file.getContents().empty()) {
        Conf.calCenter(1, Conf._numRes);
    }
}

// Exits program if any mutation is not within simulated region
// Must be called prior to correcting to internal residue indices
void MutationManager::verify_all_in_sim_region(const std::vector<int> &sTart, const std::vector<int> &eNd) const {

    for (auto i = m_cfg_file.getContents().cbegin();
         i != m_cfg_file.getContents().cend();
         ++i) {

        const int mut_idx = stoi(i->first);

        if (!is_in_simulated_region(mut_idx, sTart, eNd)) {
            std::cout << "ERROR - mutation " << mut_idx << " is not within simulated region.\n";
            exit(0);
        }
    }
}
