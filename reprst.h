// dsm.h
//header file for discrete state representation.

#ifndef _REPRST_
#define _REPRST_

#include "atom.h"

// Data files
#define FILE_SCTORSION2 "data/SCT_PF.txt"

// Parameters
#define NUM_RES_TP 20 // number of amino acids

// side chain representation
class SCR {
public:
    static FIMAP SCRMAP[NUM_RES_TP];   // same as BBRMAP
    static void InitSCAng(const char *scFile);
}; // SCR

#endif
