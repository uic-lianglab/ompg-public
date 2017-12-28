// sc_util.cpp
// useful functions for side-chain modeling program

#include "util.h"
#include "structure.h"
#include "matrix.h"

#include <iostream>

// angle formed by ab and bc
double angle(const Point a, const Point b, const Point c) {
    double d12 = b.dis(a);
    double d23 = b.dis(c);
    double dot, angle;
    dot = (a - b).dot(b - c);
    angle = 180 * acos(dot / (d12 * d23)) / PI;
    return (180 - angle);
}

// torsion angle formed by four points a,b,c,and d.
// @TODO - should Point parameters be passed by reference?
// @return - torsion angle in radians [-PI,PI]
double torsion(const Point a, const Point b, const Point c, const Point d, const eTorsionRadians angtype) {
    double xij, yij, zij,
           xkj, ykj, zkj,
           xkl, ykl, zkl,
           dxi, dyi, dzi,
           gxi, gyi, gzi,
           bi, bk, ct,
           z1, z2, ap, s;
    double xi = a.x, yi = a.y, zi = a.z;
    double xj = b.x, yj = b.y, zj = b.z;
    double xk = c.x, yk = c.y, zk = c.z;
    double xl = d.x, yl = d.y, zl = d.z;

    /* Calculate the vectors C,B,C                                       */
    xij = xi - xj;
    yij = yi - yj;
    zij = zi - zj;
    xkj = xk - xj;
    ykj = yk - yj;
    zkj = zk - zj;
    xkl = xk - xl;
    ykl = yk - yl;
    zkl = zk - zl;

    /* Calculate the normals to the two planes n1 and n2
       this is given as the cross products:
       AB x BC
       --------- = n1
       |AB x BC|

       BC x CD
       --------- = n2
       |BC x CD|
    */
    dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
    dyi = zij * xkj - xij * zkj;
    dzi = xij * ykj - yij * xkj;
    gxi = zkj * ykl - ykj * zkl;     /* Normal to plane 2                */
    gyi = xkj * zkl - zkj * xkl;
    gzi = ykj * xkl - xkj * ykl;

    /* Calculate the length of the two normals                           */
    bi = dxi * dxi + dyi * dyi + dzi * dzi;
    bk = gxi * gxi + gyi * gyi + gzi * gzi;
    ct = dxi * gxi + dyi * gyi + dzi * gzi;

    bi = (double) sqrt((double) bi);
    bk = (double) sqrt((double) bk);

    z1 = 1. / bi;
    z2 = 1. / bk;
    ct = ct * z1 * z2;
    if (ct > 1.0) ct = 1.0;
    if (ct < (-1.0)) ct = -1.0;
    ap = acos(ct);

    s = xkj * (dzi * gyi - dyi * gzi)
        + ykj * (dxi * gzi - dzi * gxi)
        + zkj * (dyi * gxi - dxi * gyi);

    if (s < 0.0) ap = -ap;

    ap = (ap > 0.0) ? PI - ap : -(PI + ap);

    return ap;
}

// torsion angle formed by four points a,b,c,and d.
// @TODO - should Point parameters be passed by reference?
// @return - torsion angle in degrees [-180,180]
double torsion(const Point a, const Point b, const Point c, const Point d, const eTorsionDegrees angtype) {
    return TO_DEG(torsion(a, b, c, d, torsion_radians));
}

int split(const string &str, const char delim, vector<string> &results) {
    size_t start = str.find_first_not_of(delim), end = start;
    while (start != std::string::npos){
        // Find next occurrence of delimiter
        end = str.find(delim, start);
        // Push back the token found into vector
        results.push_back(str.substr(start, end - start));
        // Skip all occurrences of the delimiter to find new start
        start = str.find_first_not_of(delim, end);
    }
    return static_cast<int>(results.size());
}

string itoa(const int num) {
    stringstream converter;
    converter << num;
    return converter.str();
}

string ftoa(const double num) {
    stringstream converter;
    converter << num;
    return converter.str();
}

// calculate position of n_res given the previous positions of three atoms and the bond length, bond angle and torsion angle
void calCo(const Point *prev_atoms, double length, double bAngle, double tAngle,
           Point &n_res) {
    Point su, u3, _SvdV;
    Point _u, _v, _w, cur_res, last_res, bfl_res;
    double d, dis2;
    cur_res = prev_atoms[2];
    last_res = prev_atoms[1];
    bfl_res = prev_atoms[0];
    _SvdV = bfl_res - last_res;
    su = cur_res - last_res;
    d = sqrt(su.x * su.x + su.y * su.y + su.z * su.z);
    _u = su / d;
    d = _SvdV.dot(_u);
    u3 = _u * d;
    dis2 = _SvdV.dis(u3);
    _v = (_SvdV - u3) / dis2;
    _w = _u.cross(_v);
    n_res = cur_res + _u * length * cos(PI - bAngle)
            + _v * length * sin(PI - bAngle) * cos(tAngle)
            + _w * length * sin(PI - bAngle) * sin(tAngle);
}

template<class T>
void quicksort(T a[], const int &leftarg,
               const int &rightarg) {
    if (leftarg < rightarg) {

        T pivotvalue = a[leftarg];
        int left = leftarg - 1;
        int right = rightarg + 1;

        for (; ;) {

            while (a[--right] > pivotvalue);
            while (a[++left] < pivotvalue);

            if (left >= right) break;

            T temp = a[right];
            a[right] = a[left];
            a[left] = temp;
        }

        int pivot = right;
        quicksort(a, leftarg, pivot);
        quicksort(a, pivot + 1, rightarg);
    }
}

int SampleOne(const vector<double> &prob) {
    double total = 0.0;
    const int nprob = (int) prob.size();
    for (int i = 0; i < nprob; ++i)
        total += prob[i];
    double chosentotal = drand(0, total);
    total = 0.0;
    for (int i = 0; i < nprob; ++i) {
        total += prob[i];
        if (total > chosentotal) {
            return (i);
        }
    }
    // We reach here because drand adds +.5 which causes overflow issues
    // @TODO - fix this using actual random library in C++11
    // @HACK - for now just return last element
    return nprob - 1;
}

vector<string> FileLines(const string &inFile) {
    // Returns a vector with strings representing the lines of a file
    vector<string> inFileV;
    ifstream inFileS;
    inFileS.open(inFile.c_str(), ios::in);
    inFileS.clear();
    inFileS.seekg(0);
    while (!inFileS.eof()) {
        string line;
        getline(inFileS, line);
        if (line.size() > 0)
            inFileV.push_back(line);
    }
    return inFileV;
}

string File2ProtName(const string &fname) {
    vector<string> fnameTok;
    split(fname, '/', fnameTok);
    string lastTok = fnameTok[fnameTok.size() - 1];
    int end = 0;
    for (end = 0; end <= 6; end++)
        if (lastTok.substr(end, 1) == ".")
            break;
    return (lastTok.substr(0, end));
}

bool FileExists(const string &fname) {
    std::ifstream fin(fname);
    // destructor call should auto close stream
    return fin.good();
}

void box_MullerNsample_single(double &normal_Svalue, const double mean, const double sigma) {
    double U, V;
    U = frand(0, 1);
    V = frand(0, 1);
    normal_Svalue = sqrt((-2) * log(U)) * cos(2 * PI * V) * sigma + mean;
    if (normal_Svalue > PI) {
        normal_Svalue = normal_Svalue - mean - mean;
    }
}

bool sortfun_E(combiEIndex x, combiEIndex y) {
    return (x.energy < y.energy);
}

// To get the atom name in each line in PDB file
string PDB_atomname_Extraction(const string &tmpStr) {
    string atomName;
    if (tmpStr.substr(12, 1) != " ")
        atomName += tmpStr.substr(12, 4);
    else if (tmpStr.substr(14, 1) == " ")
        atomName += tmpStr.substr(13, 1);
    else if (tmpStr.substr(15, 1) == " ")
        atomName += tmpStr.substr(13, 2);
    else atomName += tmpStr.substr(13, 3);
    return atomName;
}

// Read the PDB file generated by myself only store loop region
void Loop_DecoyPDB_read(const string &PathparPDB, vector<vector<string> > &TMPCC, int conf_Upbound) {
    vector<string> TMPC;
    string tmpStr;
    int confCount = 0;
    ifstream INPUTFILE;
    INPUTFILE.open(PathparPDB.c_str(), ios::in);
    if (!INPUTFILE.is_open()) {
        cout << "cannot open pdb file " << PathparPDB << endl;
        exit(0);
    }

    while (!INPUTFILE.eof()) {
        getline(INPUTFILE, tmpStr);
        if (tmpStr.substr(0, 6) == "MODEL ") {
            if (confCount < conf_Upbound) {
                confCount++;
            }
            else
                break;
        }
        else if (tmpStr.substr(0, 6) == "ATOM  ") {
            TMPC.push_back(tmpStr);
        }
        else if (tmpStr.substr(0, 6) == "ENDMDL") {
            TMPCC.push_back(TMPC);
            TMPC.clear();

        }
        else
            continue;
    }

    INPUTFILE.close();
}

// For multi-loops
void SCE_Minimization_list(const Structure &conf, const vector<int> &start, const vector<int> &end, const vector<int> &ResIdx,
                           vector<int> &ClashNum, const std::vector<int> &List, const double basic_rot)
{
    int i, j, k, p, direction, posn, resType, collnum, rot_per_unit;
    int minCollnum = 10000;
    int size = static_cast<int>(ResIdx.size());
    int loopidx = 0;
    double pri_x, pri_y, pri_z, R, change_rot = 0;
    Atom rot_axis, tmp;
    Residue tmpRes, tmpRes2;
    double **trans_matrix;
    trans_matrix = matrix(1, 3, 1, 3);

    for (i = 0; i < size; i++) {
        minCollnum = 10000;
        if (ClashNum[i] == 0)
            continue;
        posn = ResIdx[i];
        for (p = 0; p < start.size(); p++)
            if (posn >= start[p] && posn <= end[p]) {
                loopidx = p;
                break;
            }
        resType = conf._res[posn]._type;
        tmpRes = conf._res[posn];
        if (resType == PRO)   // Proline has special treatment
        {
            rot_axis.x = conf._res[posn]._atom[ATM_CA].x - conf._res[posn]._atom[ATM_N].x;
            rot_axis.y = conf._res[posn]._atom[ATM_CA].y - conf._res[posn]._atom[ATM_N].y;
            rot_axis.z = conf._res[posn]._atom[ATM_CA].z - conf._res[posn]._atom[ATM_N].z;
            R = sqrt(rot_axis.x * rot_axis.x + rot_axis.y * rot_axis.y + rot_axis.z * rot_axis.z);
            rot_axis = rot_axis / R;
            for (k = 1; k <= 35; k++) {
                for (direction = -1; direction <= 1; direction += 2) {
                    rot_per_unit = k * direction;
                    change_rot = k * basic_rot * direction;
                    RotMatrix(trans_matrix, (double) rot_axis.x, (double) rot_axis.y, (double) rot_axis.z,
                              (double) change_rot);
                    for (j = NUM_BB_ATOM - 1; j < conf._res[posn]._numAtom; ++j) {
                        if (conf._res[posn]._atom[j]._type == UNDEF || conf._res[posn]._atom[j]._type >= 22) continue;
                        else {
                            tmp = conf._res[posn]._atom[j] - conf._res[posn]._atom[ATM_CA];
                            pri_x = tmp.x;
                            pri_y = tmp.y;
                            pri_z = tmp.z;
                            tmp.x = (double) (trans_matrix[1][1] * pri_x + trans_matrix[1][2] * pri_y +
                                              trans_matrix[1][3] * pri_z);
                            tmp.y = (double) (trans_matrix[2][1] * pri_x + trans_matrix[2][2] * pri_y +
                                              trans_matrix[2][3] * pri_z);
                            tmp.z = (double) (trans_matrix[3][1] * pri_x + trans_matrix[3][2] * pri_y +
                                              trans_matrix[3][3] * pri_z);
                            tmpRes._atom[j] = conf._res[posn]._atom[ATM_CA] + tmp;
                        }
                    }
                    collnum = Res_clash_detection_list(conf, tmpRes, start, end, loopidx, List);

                    if (collnum == 0) {
                        minCollnum = 0;
                        tmpRes2 = tmpRes;
                        break;
                    }
                    else {
                        if (collnum < minCollnum) {
                            minCollnum = collnum;
                            tmpRes2 = tmpRes;
                        }
                    }
                }
                if (minCollnum == 0)
                    break;
            }
            conf._res[posn] = tmpRes2;
            ClashNum[i] = minCollnum;
        }
        else {
            rot_axis.x = conf._res[posn]._atom[ATM_CB].x - conf._res[posn]._atom[ATM_CA].x;
            rot_axis.y = conf._res[posn]._atom[ATM_CB].y - conf._res[posn]._atom[ATM_CA].y;
            rot_axis.z = conf._res[posn]._atom[ATM_CB].z - conf._res[posn]._atom[ATM_CA].z;
            R = sqrt(rot_axis.x * rot_axis.x + rot_axis.y * rot_axis.y + rot_axis.z * rot_axis.z);
            rot_axis = rot_axis / R;
            for (k = 1; k <= 35; k++) {
                for (direction = -1; direction <= 1; direction += 2) {
                    rot_per_unit = k * direction;
                    change_rot = k * basic_rot * direction;
                    RotMatrix(trans_matrix, (double) rot_axis.x, (double) rot_axis.y, (double) rot_axis.z,
                              (double) change_rot);
                    for (j = NUM_BB_ATOM; j < conf._res[posn]._numAtom; ++j) {
                        if (conf._res[posn]._atom[j]._type == UNDEF || conf._res[posn]._atom[j]._type >= 22) continue;
                        else {
                            tmp = conf._res[posn]._atom[j] - conf._res[posn]._atom[ATM_CB];
                            pri_x = tmp.x;
                            pri_y = tmp.y;
                            pri_z = tmp.z;
                            tmp.x = (double) (trans_matrix[1][1] * pri_x + trans_matrix[1][2] * pri_y +
                                              trans_matrix[1][3] * pri_z);
                            tmp.y = (double) (trans_matrix[2][1] * pri_x + trans_matrix[2][2] * pri_y +
                                              trans_matrix[2][3] * pri_z);
                            tmp.z = (double) (trans_matrix[3][1] * pri_x + trans_matrix[3][2] * pri_y +
                                              trans_matrix[3][3] * pri_z);
                            tmpRes._atom[j] = conf._res[posn]._atom[ATM_CB] + tmp;
                        }
                    }
                    collnum = Res_clash_detection_list(conf, tmpRes, start, end, loopidx, List);
                    if (collnum == 0) {
                        minCollnum = 0;
                        tmpRes2 = tmpRes;
                        break;
                    }
                    else {
                        if (collnum < minCollnum) {
                            minCollnum = collnum;
                            tmpRes2 = tmpRes;
                        }
                    }
                }
                if (minCollnum == 0)
                    break;
            }
            conf._res[posn] = tmpRes2;
            ClashNum[i] = minCollnum;
        }
    }

    free_matrix(trans_matrix, 1, 3, 1, 3);
}

void SCE_Minimization_list(const Structure &conf, const int start, const int end, const vector<int> &ResIdx, vector<int> &ClashNum,
                           const std::vector<int> &List, const double basic_rot)
{
    int i, j, k, direction, posn, resType, collnum, minIdx, rot_per_unit;
    int minCollnum = 10000;
    int size = static_cast<int>(ResIdx.size());
    double pri_x, pri_y, pri_z, R, change_rot = 0;
    Atom rot_axis, tmp;
    Residue tmpRes, tmpRes2;
    double **trans_matrix;
    trans_matrix = matrix(1, 3, 1, 3);

    for (i = 0; i < size; i++) {
        minCollnum = 10000;
        if (ClashNum[i] == 0)
            continue;
        posn = ResIdx[i];
        resType = conf._res[posn]._type;
        tmpRes = conf._res[posn];
        if (resType == 12)   // Proline has special treatment
        {
            rot_axis.x = conf._res[posn]._atom[ATM_CA].x - conf._res[posn]._atom[ATM_N].x;
            rot_axis.y = conf._res[posn]._atom[ATM_CA].y - conf._res[posn]._atom[ATM_N].y;
            rot_axis.z = conf._res[posn]._atom[ATM_CA].z - conf._res[posn]._atom[ATM_N].z;
            R = sqrt(rot_axis.x * rot_axis.x + rot_axis.y * rot_axis.y + rot_axis.z * rot_axis.z);
            rot_axis = rot_axis / R;
            for (k = 1; k <= 35; k++) {
                for (direction = -1; direction <= 1; direction++) {
                    if (direction == 0)
                        continue;
                    rot_per_unit = k * direction;
                    change_rot = k * basic_rot * direction;
                    RotMatrix(trans_matrix, (double) rot_axis.x, (double) rot_axis.y, (double) rot_axis.z,
                              (double) change_rot);
                    for (j = NUM_BB_ATOM - 1; j < conf._res[posn]._numAtom; ++j) {
                        if (conf._res[posn]._atom[j]._type == UNDEF || conf._res[posn]._atom[j]._type >= 22) continue;
                        else {
                            tmp = conf._res[posn]._atom[j] - conf._res[posn]._atom[ATM_CA];
                            pri_x = tmp.x;
                            pri_y = tmp.y;
                            pri_z = tmp.z;
                            tmp.x = (double) (trans_matrix[1][1] * pri_x + trans_matrix[1][2] * pri_y +
                                              trans_matrix[1][3] * pri_z);
                            tmp.y = (double) (trans_matrix[2][1] * pri_x + trans_matrix[2][2] * pri_y +
                                              trans_matrix[2][3] * pri_z);
                            tmp.z = (double) (trans_matrix[3][1] * pri_x + trans_matrix[3][2] * pri_y +
                                              trans_matrix[3][3] * pri_z);
                            tmpRes._atom[j] = conf._res[posn]._atom[ATM_CA] + tmp;
                        }
                    }
                    collnum = Res_clash_detection_list(conf, tmpRes, start, end, List);
                    if (collnum == 0) {
                        minCollnum = 0;
                        tmpRes2 = tmpRes;
                        minIdx = rot_per_unit;
                        break;
                    }
                    else {
                        if (collnum < minCollnum) {
                            minCollnum = collnum;
                            tmpRes2 = tmpRes;
                            minIdx = rot_per_unit;
                        }
                    }
                }
                if (minCollnum == 0)
                    break;
            }
            conf._res[posn] = tmpRes2;
            ClashNum[i] = minCollnum;
        }
        else {
            rot_axis.x = conf._res[posn]._atom[ATM_CB].x - conf._res[posn]._atom[ATM_CA].x;
            rot_axis.y = conf._res[posn]._atom[ATM_CB].y - conf._res[posn]._atom[ATM_CA].y;
            rot_axis.z = conf._res[posn]._atom[ATM_CB].z - conf._res[posn]._atom[ATM_CA].z;
            R = sqrt(rot_axis.x * rot_axis.x + rot_axis.y * rot_axis.y + rot_axis.z * rot_axis.z);
            rot_axis = rot_axis / R;
            for (k = 1; k <= 35; k++) {
                for (direction = -1; direction <= 1; direction++) {
                    if (direction == 0)
                        continue;
                    rot_per_unit = k * direction;
                    change_rot = k * basic_rot * direction;
                    RotMatrix(trans_matrix, (double) rot_axis.x, (double) rot_axis.y, (double) rot_axis.z,
                              (double) change_rot);
                    for (j = NUM_BB_ATOM; j < conf._res[posn]._numAtom; ++j) {
                        if (conf._res[posn]._atom[j]._type == UNDEF || conf._res[posn]._atom[j]._type >= 22) continue;
                        else {
                            tmp = conf._res[posn]._atom[j] - conf._res[posn]._atom[ATM_CB];
                            pri_x = tmp.x;
                            pri_y = tmp.y;
                            pri_z = tmp.z;
                            tmp.x = (double) (trans_matrix[1][1] * pri_x + trans_matrix[1][2] * pri_y +
                                              trans_matrix[1][3] * pri_z);
                            tmp.y = (double) (trans_matrix[2][1] * pri_x + trans_matrix[2][2] * pri_y +
                                              trans_matrix[2][3] * pri_z);
                            tmp.z = (double) (trans_matrix[3][1] * pri_x + trans_matrix[3][2] * pri_y +
                                              trans_matrix[3][3] * pri_z);
                            tmpRes._atom[j] = conf._res[posn]._atom[ATM_CB] + tmp;
                        }
                    }
                    collnum = Res_clash_detection_list(conf, tmpRes, start, end, List);
                    if (collnum == 0) {
                        minCollnum = 0;
                        tmpRes2 = tmpRes;
                        minIdx = rot_per_unit;
                        break;
                    }
                    else {
                        if (collnum < minCollnum) {
                            minCollnum = collnum;
                            tmpRes2 = tmpRes;
                            minIdx = rot_per_unit;
                        }
                    }
                }
                if (minCollnum == 0)
                    break;
            }
            conf._res[posn] = tmpRes2;
            ClashNum[i] = minCollnum;
        }
    }

    free_matrix(trans_matrix, 1, 3, 1, 3);
}

double Ellipsoid_Detect(Structure &Conf, int Start, int End, bool type) {

    // determine 2a
    double ConstL = Residue::bond_length[Conf._res[Start]._type][ATM_C] + 3 * (End - Start);  // 3A is the C-C distance
    double ConstLtotal = 0;
    if (type == false)
        ConstLtotal = ConstL;
    else {
        int MAX_SC_SIZE = 7;
        double dist = Conf._res[Start]._atom[ATM_CA].dis(Conf._res[End]._atom[ATM_C]);

        ConstLtotal = sqrt((sqrt((ConstL / 2) * (ConstL / 2) - (dist / 2) * (dist / 2)) + MAX_SC_SIZE) *
                           (sqrt((ConstL / 2) * (ConstL / 2) - (dist / 2) * (dist / 2)) + MAX_SC_SIZE) +
                           (dist / 2) * (dist / 2)) * 2;
    }
    return ConstLtotal;
}
