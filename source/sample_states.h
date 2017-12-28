// sample_angles.h

#ifndef _SAMPLE_ANGLES_
#define _SAMPLE_ANGLES_

#include "structure.h"
#include "rotamer.h"
#include "reprst.h"
#include "util.h"
#include "vla.h"

#include <iostream>

// @HACK - Visual studio does not support variable length arrays and underneath
// the hood they are standard vectors. This causes function signature issues
// as sometimes a compile time C-style array is actually passed to this method.
// In visual studio, can just use a template method to handle both of these
// situations. GCC is less forgiving, so we will have to declare the
// degenerate pointer array prototype + implementation in this case.

#if _MSC_VER

template<typename t_>
void sample_sc_angles(const Residue &res, t_ &angles,
                      int numStates, int scType) {
    int i, j, resType, numRot;

    if (scType == 1 || scType == 2) { // continuous representation
        int key;
        double r;
        FIMAP::iterator fiItr;
        resType = res._type;
        if (resType == -1) resType = 1;
        numRot = Rotamer::numRotBond[resType];
        for (i = 0; i < numStates; ++i) {
            r = frand(0, 1);
            if ((fiItr = SCR::SCRMAP[resType].upper_bound(r)) != SCR::SCRMAP[resType].end()) {
                key = fiItr->second;
                for (j = (numRot - 1); j >= 0; --j) {
                    angles[i][j] = (key % 100 * SC_T_INT + frand(0, SC_T_INT) - 180) * PI / 180;
                    key = key / 100;
                }

            }
            else {
                std::cout << "Warning: SCRMAP out of bound, resample the angle. " << res.get_name_3() << " " << r << " "
                          << SCR::SCRMAP[resType].size() << std::endl;
                --i;
            }
        }
    }
    else if (scType == 3) { // sample the +- 60 degree of the native torsion angles
        double r, a, b;
        resType = res._type;
        numRot = Rotamer::numRotBond[resType];
        for (i = 0; i < numStates; ++i) {
            for (j = 0; j < numRot; ++j) {
                a = res._scChi[j] - 10;
                b = res._scChi[j] + 10;
                r = frand(a, b);
                if (r < -180) r = r + 360;
                else if (r > 180) r = r - 360;
                angles[i][j] = r * PI / 180;
            }
        }
    }
    else {
        std::cout << "Wrong Sidechain sampling type!!!" << std::endl;
        exit(0);
    }
}

#else

void sample_sc_angles(const Residue &res, double(*angles)[6],
                      int numStates, int scType);

#endif // __MSC_VER

#endif // _SAMPLE_ANGLES_
