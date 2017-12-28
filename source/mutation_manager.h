//****************************************************************************
// mutation_manager.h
//
// Class for reading and applying mutation options from .ini
//****************************************************************************

#ifndef mutation_manager_h
#define mutation_manager_h

//****************************************************************************
// Includes
//****************************************************************************

// For parsing .ini files
#include "config/ConfigFile.h"

#include <string>
#include <vector>

/**
 * Reads and applies mutation configuration parameters
 */
class MutationManager
{
public:

    // Loads the mutations from .ini file
    // @param ini_fpath - string with path to .ini file
    // @return TRUE if mutations read properly, FALSE otherwise
    void load(const std::string &ini_fpath);

    // Modifies parameter structure with configured mutations.
    // Note: all atoms of mutated residues are moved to origin.
    // @param Conf - proteint structure to apply mutations to
    void apply(class Structure &Conf) const;

    // Exits program if any mutation is not within simulated region
    // Must be called prior to correcting to internal residue indices
    void verify_all_in_sim_region(const std::vector<int> &sTart, const std::vector<int> &eNd) const;

private:

    ConfigFile m_cfg_file;
};

#endif // mutation_manager_h
