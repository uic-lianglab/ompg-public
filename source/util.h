// util.h
// header file for util.cpp

#ifndef _UTIL_
#define _UTIL_

#include <sstream>
#include <limits>
#include <assert.h>
#include <algorithm>

#include "atom.h"
#include "structure.h"
#include "cal_energy.h"

struct combiEIndex {
    double energy;
    int index;
};

double angle(const Point, const Point, const Point);

enum eTorsionRadians {torsion_radians};
enum eTorsionDegrees {torsion_degrees};
double torsion(const Point a, const Point b, const Point c, const Point d, const eTorsionRadians angtype);
double torsion(const Point a, const Point b, const Point c, const Point d, const eTorsionDegrees angtype);

int split(const string &input, const char delim, vector<string> &results);

string itoa(const int num);

string ftoa(const double num);

void calCo(const Point *prev_atoms, double length, double bAngle, double tAngle, Point &n_res);

template<class T>
void quicksort(T a[], const int &leftarg, const int &rightarg);

int SampleOne(const vector<double> &prob);

vector<string> FileLines(const string &inFile);

string File2ProtName(const string &fname);

bool FileExists(const string &fname);

void box_MullerNsample_single(double &normal_Svalue, const double mean, const double sigma);

void Loop_DecoyPDB_read(const string &PathparPDB, vector<vector<string> > &TMPCC, int conf_Upbound);

double Ellipsoid_Detect(Structure &Conf, int Start, int End, bool type);

bool sortfun_E(combiEIndex x, combiEIndex y);

template<typename T>
std::string to_string(T const &value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}

inline bool check_real_eq(const double a, const double b) {
    return (std::abs(a - b) < ((2.0) * std::numeric_limits<float>::epsilon()));
}

inline bool check_real_gte(const double a, const double b) {
    return (a > b) || check_real_eq(a, b);
}

// @return TRUE if res_ix is in a simulated loop region, FALSE o/w
inline bool is_in_simulated_region(const int res_ix,
                                   const std::vector<int> &starts,
                                   const std::vector<int> &ends) {
    // residue indices start 1
    assert(res_ix >= 1);
    assert(starts.size() == ends.size());
    bool is_in_loop_region = false;
    const int num_loops = (int) starts.size();
    for (int i = 0; i < num_loops; ++i) {
        assert(starts[i] < ends[i]);
        if (res_ix >= starts[i] && res_ix <= ends[i]) {
            is_in_loop_region = true;
            break;
        }
    }
    return is_in_loop_region;
}

// Crops string to keep 'precision' number of digits past decimal point
inline std::string ftoa(double real, const size_t precision) {
    std::string real_str = ftoa(real);
    const size_t dec_point = real_str.find(".");
    size_t str_len = real_str.size();
    if (dec_point != real_str.npos) {
        str_len = std::min(dec_point + precision + 1, real_str.size());
    }
    return real_str.substr(0, str_len);
}

// Remove leading and trailing instances of parameter characters
inline void trim(std::string& s, const std::string &char_set) {
    s.erase(0, s.find_first_not_of(char_set));
    s.erase(s.find_last_not_of(char_set) + 1U);
}

// Remove leading and trailing whitespace
inline void trim(std::string& s) {
    static const char whitespace[] = " \n\t\v\r\f";
    s.erase(0, s.find_first_not_of(whitespace));
    s.erase(s.find_last_not_of(whitespace) + 1U);
}

#endif
