// reprst.cpp

#include "reprst.h"
#include "util.h"
#include "residue.h"
#include "rotamer.h"
#include "potential.h"

#include <iostream>

FIMAP SCR::SCRMAP[NUM_RES_TP];

// initialize side chain angle distribution, which is used in sampling
// side chain torsion angles
void SCR::InitSCAng(const char *scFile) {
    int i, j, aaIndex, numRot;
    char tmpLine[500];
    ifstream inFile;
    string tmpStr;
    vector<string> strVec;
    // read side chain torsion angles from scFile
    inFile.open(scFile, ios::in);
    if (!inFile.is_open()) {
        std::cout << "cannot open side chain torsion angle file at " << scFile << endl;
        return;
    }
    int key;
    double prob, cumProb, count, totalCnt[NUM_RES_TP];
    vector<int> scDataKey[NUM_RES_TP];
    vector<double> scDataCnt[NUM_RES_TP];

    memset(&(totalCnt[0]), 0, sizeof(double) * NUM_RES_TP);

    while (!inFile.eof()) {
        inFile.getline(tmpLine, 500);
        tmpStr = tmpLine;
        if (tmpStr.empty()) continue;
        if (tmpStr[0] == '#') continue;
        strVec.clear();
        split(tmpStr, ' ', strVec);
        aaIndex = atoi(strVec[0].c_str());
        numRot = Rotamer::numRotBond[aaIndex];
        if (strVec.size() < (numRot + 3)) continue;
        key = 0;
        for (j = 1; j <= numRot; ++j)
            key = key * 100 + atoi(strVec[j].c_str());
        count = atof(strVec[numRot + 2].c_str());
        scDataKey[aaIndex].push_back(key);
        scDataCnt[aaIndex].push_back(count);
        totalCnt[aaIndex] += count;
    }
    for (i = 0; i < NUM_RES_TP; ++i) {
        if (scDataKey[i].size() == 0) continue;
        cumProb = 0;
        for (j = 0; j < scDataKey[i].size(); ++j) {
            prob = scDataCnt[i][j] / totalCnt[i];
            cumProb += prob;
            SCRMAP[i].insert(FIMAP::value_type(cumProb, scDataKey[i][j]));
        }

        // make sure the last one has cumulative probability equals to 1
        if (cumProb < 1)
            SCRMAP[i].insert(FIMAP::value_type(1,
                                               scDataKey[i][scDataKey[i].size() - 1]));
    }
}

