// rotamer.cpp

#include "rotamer.h"
#include "residue.h"
#include <new>

const int Rotamer::numRotBond[20] = {
        0, 1, 2, 3, 2,    // A C D E F
        0, 2, 2, 4, 2,    // G H I K L
        3, 2, 2, 3, 4,    // M N P Q R
        1, 1, 1, 2, 2     // S T V W Y
};
