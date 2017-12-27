// cal_energy.h
// header file for cal_energy.cpp
// Dec 20, 2005

#ifndef _CALE_
#define _CALE_

#include <vector>

using namespace std;

#include "structure.h"
#include "util.h"
#include "potential.h"
#include "rotamer.h"

// During side chain clash resolution, torsion angles are incremented according
// to this value (currently hard-coded to 5 degrees)
// @TODO - consider making this command line configurable
#define SCE_MIN_DELTA_ROT (PI / 36.0)

class Structure;

extern double one_res_en_loodis_bb2all(const Structure &conf, const Residue &res, const int start, const int end,
                                       const int Start, const int End, const int type = 0);

// Computes potential of param residue backbone atoms against all atoms in param List
extern double one_res_en_loodis_bb2all_list(const Structure &conf, const Residue &res,
                                            const int Start, const int End, const int type, const vector<int> &List);

extern double one_res_en_loodis_sc(const Structure &conf, const Residue &res, const int position, const int Start,
                                   const int End);

extern double loodis_e_list(const Structure &Conf, const vector<int> &starts, const vector<int> &ends,
                            const vector<int> &List, const int type);

extern double loodis_e(const Structure &Conf, const vector<int> &starts, const vector<int> &ends, const bool type);

extern void Clash_detection_list(const Structure &conf, const int Start, const int End, vector<int> &ResIdx,
                                 vector<int> &ClashNum, const vector<int> &List);

extern void Clash_detection_list(const Structure &conf, const vector<int> &starts, const vector<int> &ends, vector<int> &ResIdx,
                                 vector<int> &ClashNum, const vector<int> &List);

extern void BBClash_detection_list(const Structure &conf, const vector<int> &starts, const vector<int> &ends,
                                   vector<int> &ResIdx, vector<int> &ClashNum, const vector<int> &List);

extern void BBClash_detection_list(const Structure &conf, const int Start, const int End, vector<int> &ResIdx,
                                   vector<int> &ClashNum, const vector<int> &List);

extern int Res_clash_detection_list(const Structure &conf, const Residue &res,
                                    const int Start, const int End, const vector<int> &List);

extern int Res_clash_detection_list(const Structure &conf, const Residue &res, const vector<int> &starts,
                                    const vector<int> &ends, const int loopidx, const vector<int> &List);

extern void SCE_Minimization_list(const Structure &conf, const int start, const int end, const vector<int> &ResIdx, vector<int> &ClashNum,
                                  const std::vector<int> &List, const double basic_rot);

extern void SCE_Minimization_list(const Structure &conf, const vector<int> &start, const vector<int> &end, const vector<int> &ResIdx,
                                  vector<int> &ClashNum, const std::vector<int> &List, const double basic_rot);

#endif // _CALE_
