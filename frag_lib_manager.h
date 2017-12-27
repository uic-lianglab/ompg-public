//
// Created by apr on 11/4/15.
//

#ifndef MDISGRO_FRAG_LIB_MANAGER_H
#define MDISGRO_FRAG_LIB_MANAGER_H

// For parsing .ini files
#include "config/ConfigFile.h"
#include "structure.h"

#include <vector>
#include <string>

// A minimalist structure definition
typedef std::vector<Point> SlimStruct;

// Reads fragment libraries from disk
class FragLibManager {
public:

    // Loads the mutations from .ini file
    // Loops are assumed to be in same order as those in the loop region/interval file (mulist)
    // @param MultiLooplib - stores the the fragment libraries read from disk
    // @param ini_fpath - .ini file with names of fragment libraries for each individual loop modeled
    //  the library is assumed to reside within lib_dir
    // @param lib_dir - the directory containing all fragment libraries
    // @param sTart - the loop start indices
    // @param eNd - the loop end indices
    static void load(std::vector<std::vector<Structure> > &MultiLooplib, const std::string &ini_fpath,
                     const std::string &lib_dir, const vector<int> &sTart, const vector<int> &eNd);

    // Loads a single fragment library located at target path into Looplib
    static void load(std::vector<Structure> &Looplib, const std::string &frag_path);

    // Loads a fragment library into a SlimStruct data structure
    static void load(std::vector<SlimStruct> &Looplib, const std::string &frag_path, const bool ignore_H_atoms, const bool backbone_only);
};


#endif //MDISGRO_FRAG_LIB_MANAGER_H
