// potential.cpp

#include "potential.h"
#include "atom.h"
#include "structure.h"
#include "rotamer.h"
#include "util.h"
#include "cal_energy.h"
#include "reprst.h"
#include <iostream>
#include <sstream>

bool PF::cal[ENERGY_MODES];
double PF::LOODIS[20][20][LOODIS_DIS_BIN];
vector<vector<int> > PF::VDW_AA[20];
map<string, IDMAP> PF::Parameter;

void PF::InitPar(string ParFile) {
    int j, k;
    // Set up maps to fold parameters
    IDMAP EMap, AMap, BMap, HMap, SMap;
    Parameter["E"] = EMap;
    Parameter["A"] = AMap;
    Parameter["B"] = BMap;
    Parameter["H"] = HMap;
    Parameter["S"] = SMap;
    for (int i = 0; i < ENERGY_TYPES; i++)
        PF::Parameter["E"][i] = 0.0;
    // Read parameters and load into maps
    if (ParFile != "none") {
        vector<string> ParFileV = FileLines(ParFile);
        for (int i = 0; i < ParFileV.size(); i++) {
            if (ParFileV[i][0] == '#') continue;
            SVEC bob;
            split(ParFileV[i], ':', bob);
            string id = ParFileV[i].substr(0, 1);
            bob[0].erase(0, 1);
            int key = atoi(bob[0].c_str());
            double value = atof(bob[1].c_str());
            if (value != 0)
                Parameter[id][key] = value;
        }
    }

    if (cal[EM_VDW]) {
        // initialize VDW_AA for VDW energy within an amino acids
        vector<int> tmpV;
        // set number of atoms for each AA that need to calculate VDW
        int VDW_CNT[20] = {
                0, 1, 3, 4, 3,
                0, 3, 3, 5, 3,
                3, 3, 0, 4, 6,
                1, 1, 1, 4, 3
        };
        for (int i = 0; i < 20; ++i) {
            // push 7 vector<int> to each VDW_AA[i], most of them will be empty
            // VDW_AA[i][ATM_N] will store the atoms that need to calculate
            // VDW energy with ATM_N of residue i
            for (j = 0; j < 7; ++j) {
                VDW_AA[i].push_back(tmpV);
            }
            if (VDW_CNT[i] >= 1) {
                for (k = NUM_BB_ATOM; k < Residue::numAtom[i]; k++)
                    VDW_AA[i][ATM_O].push_back(k);
            }
            if (VDW_CNT[i] >= 3) {
                if (i == 7) {
                    // isoleucine
                    for (k = (NUM_BB_ATOM + 2); k < Residue::numAtom[i]; ++k) {
                        VDW_AA[i][ATM_N].push_back(k);
                        VDW_AA[i][ATM_C].push_back(k);
                    }
                }
                else {
                    for (k = (NUM_BB_ATOM + 1); k < Residue::numAtom[i]; ++k) {
                        VDW_AA[i][ATM_N].push_back(k);
                        VDW_AA[i][ATM_C].push_back(k);
                    }
                }
            }
            if (VDW_CNT[i] >= 4) {
                for (k = (NUM_BB_ATOM + 2); k < Residue::numAtom[i]; ++k) {
                    VDW_AA[i][ATM_CA].push_back(k);
                }
            }
            else if (VDW_CNT[i] >= 5) {
                for (k = (NUM_BB_ATOM + 3); k < Residue::numAtom[i]; ++k) {
                    VDW_AA[i][ATM_CB].push_back(k);
                }
            }
            else {
                for (k = (NUM_BB_ATOM + 4); k < Residue::numAtom[i]; k++)
                    VDW_AA[i][6].push_back(k);
            }
        }
    }
}

//initialize potential function LooDis
void PF::initLOODIS(string filename) {
    int i, j, k;
    double pt;
    char tmpLine[500];
    string tmpStr;
    vector<string> strVec;       // used in reading files
    string tmp;
    ifstream inFile;

    inFile.open(filename.c_str(), ios::in);

    if (!inFile.is_open()) {
        std::cout << "Error!!! Cannot open LOODIS potential file! " << filename << endl;
        exit(0);
    }

    while (!inFile.eof()) {
        inFile.getline(tmpLine, 500);
        tmpStr = tmpLine;
        if (tmpStr[0] == '#') continue;
        if (tmpStr.length() < 3) continue;
        strVec.clear();
        split(tmpStr, ' ', strVec);
        if (strVec.size() == 6) {
            i = atoi(strVec[1].c_str());
            j = atoi(strVec[2].c_str());
            k = atoi(strVec[0].c_str());
            pt = atof(strVec[5].c_str());
            LOODIS[i][j][k] = pt;
            LOODIS[j][i][k] = pt;
        }

    }
    inFile.close();
}
