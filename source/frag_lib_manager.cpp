//
// Created by apr on 11/4/15.
//

#include "frag_lib_manager.h"
#include "util.h"

#include <assert.h>

// Maximum fragment library size
#define FRAG_MAX_LIB_SIZE 10E6

// Minimum loop size in order to have a fragment library
#define FRAG_MIN_LOOP_SIZE 5

// Loads the mutations from .ini file
// Loops are assumed to be in same order as those in the loop region/interval file (mulist)
// @param MultiLooplib - stores the the fragment libraries read from disk
// @param ini_fpath - .ini file with names of fragment libraries for each individual loop modeled
//  the library is assumed to reside within lib_dir
// @param lib_dir - the directory containing all fragment libraries
// @param sTart - the loop start indices
// @param eNd - the loop end indices
void FragLibManager::load(std::vector<std::vector<Structure> > &MultiLooplib, const std::string &ini_fpath,
                          const std::string &lib_dir, const vector<int> &starts, const vector<int> &ends) {
    std::cout << "Reading frag lib config file: " << ini_fpath << std::endl;
    ConfigFile cfg;
    cfg.load(ini_fpath);

    // Currently, size of multilooplib and number of loops has to match - even if lib is empty
    const size_t n_loops = starts.size();
    assert(starts.size() == ends.size());
    MultiLooplib.clear();
    MultiLooplib.resize(n_loops);

    // Used for storing the read in library
    vector<Structure> TmpVec;
    const SSET EMPTY_SET;

    int i_loop = 0;

    // Mapping of loop_id -> frag_lib_file_name
    auto records = cfg.getContents();
    // Vector of loop_ids in order specified by ini
    auto ordered_loop_ids = cfg.getReadOrderedKeys();
    assert(records.size() == ordered_loop_ids.size());
    
    // Make sure correct number of loops specified
    if (records.size() != n_loops) {
        std::cout << "ERROR: loop count mismatch in " << ini_fpath << std::endl;
        std::cout << "\tExpected " << n_loops << " loops but " << records.size()
                  << " loops found in file.\nExiting." << std::endl;
        exit(0);
    }
    
    for (size_t i_loop = 0; i_loop < n_loops; ++i_loop) {

        // Skip past loops that are too small
        const int loop_sz = ends[i_loop] - starts[i_loop];
        if (loop_sz >= FRAG_MIN_LOOP_SIZE) {
            vector<vector<string> > TMPCC;

            const string &loop_id = ordered_loop_ids[i_loop];
            assert(records.count(loop_id));
            const string &libpdb_rec = records[loop_id]; 
            
            const string libpdb_fname = lib_dir + "/" + libpdb_rec;
            cout << "Loading fragment library for loop " << starts[i_loop] << " to " << ends[i_loop] 
                 << " from file:\n\t" << libpdb_fname << endl;

            // Parse the library pdb to populate TMPCC structure
            Loop_DecoyPDB_read(libpdb_fname, TMPCC, FRAG_MAX_LIB_SIZE);

            // Now load each individual fragment
            Structure Mconf(MAX_NUM_RES);
            cout << "\tLoop size: " << loop_sz + 1 << endl;
            cout << "\tFragment library size: " << TMPCC.size() << endl;
            for (int j = 0; j < TMPCC.size(); j++) {
                Mconf.readPdb(TMPCC[j], EMPTY_SET);
                TmpVec.push_back(Mconf);
            }

            // Store read fragment
            MultiLooplib[i_loop] = TmpVec;
            // Clear buffer so that we can read next library
            TmpVec.clear();
        } else {
            cout << "Skipping fragment library loading for loop " << starts[i_loop] << " to " << ends[i_loop]
                 << "\n\t(length " << loop_sz + 1 << " smaller than " << FRAG_MIN_LOOP_SIZE + 1 << ")" << endl;
        }

    }
}

// Loads a single fragment library located at target path into Looplib
void FragLibManager::load(std::vector<Structure> &Looplib, const std::string &frag_path) {
    cout << "Loading fragment library " << frag_path << endl;

    // Used for storing the read in library
    Looplib.clear();
    const SSET EMPTY_SET;

    // Parse the library pdb to populate TMPCC structure
    vector<vector<string> > TMPCC;
    Loop_DecoyPDB_read(frag_path, TMPCC, FRAG_MAX_LIB_SIZE);

    // Now load each individual fragment
    Structure Mconf(MAX_NUM_RES);
    cout << "\tFragment library size: " << TMPCC.size() << endl;
    for (int j = 0; j < TMPCC.size(); j++) {
        Mconf.readPdb(TMPCC[j], EMPTY_SET);
        Looplib.push_back(Mconf);
    }
}

// Loads a fragment library into a SlimStruct data structure
void FragLibManager::load(std::vector<SlimStruct> &Looplib, const std::string &frag_path, const bool ignore_H_atoms,
                          const bool backbone_only) {
    cout << "Loading fragment library " << frag_path << endl;

    Looplib.clear();
    const SSET EMPTY_SET;

    // Parse the library pdb to populate TMPCC structure
    vector<vector<string> > TMPCC;
    Loop_DecoyPDB_read(frag_path, TMPCC, FRAG_MAX_LIB_SIZE);

    // Now load each individual fragment
    Structure Mconf(MAX_NUM_RES);
    SlimStruct conf;
    cout << "\tFragment library size: " << TMPCC.size() << endl;
    for (int j = 0; j < TMPCC.size(); j++) {
        conf.clear();
        Mconf.readPdb(TMPCC[j], EMPTY_SET);
        // Convert to SlimStruct format
        for (int i_res = 1; i_res <= Mconf._numRes; ++i_res) {
            const Residue &res = Mconf._res[i_res];
            const int num_atoms = backbone_only ? std::min(res._numAtom, (short)NUM_BB_ATOM) : res._numAtom;
            for (int i_atom = 0; i_atom < num_atoms; ++i_atom) {
                if (ignore_H_atoms && i_atom == ATM_H) {
                    continue;
                }
                const Point &atom = res._atom[i_atom];
                conf.push_back(atom);
            }
        }

        Looplib.push_back(conf);
    }
}
