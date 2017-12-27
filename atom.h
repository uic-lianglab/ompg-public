/* atom.h
The basic header file for all programs.
Class Atom inherited from Point.

*/

#ifndef _ATOM_
#define _ATOM_

#include "point.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cctype>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <new>
#include <cstring>
#include <math.h>

using namespace std;
class Residue;

#define MAX_ATOM_TYPE 30  // maximum number of atom types allowed
#define MAX_NUM_ATOM_RES 30
#define NUM_RES_TYPE 30
#define MAX_NUM_RES 1500
#define PI 3.14159265358979
#define EXPO 2.718281828
#define UNDEF -12345   // undefined value
#define intrand(a,b) ((int)(a+(b-a)*((double)rand())/RAND_MAX)) // [a, b-1]
#define frand(a,b) ((rand()+.5)/RAND_MAX*(b-a)+a)
#define drand(a,b) (((double)rand()+.5)/RAND_MAX*(b-a)+a)
#define BBTbinSize 5 // degrees
#define TORBIN 360/BBTbinSize

// macros for angle conversion
#define TO_DEG(x) (x*180/PI)
#define TO_RAD(x) (x*PI/180)

// defined types
typedef set <string> SSET;               // set of strings
typedef vector <int> IVEC;               // vector of int
typedef vector <string> SVEC;            // vector of string
typedef map <int, double> IFMAP;         // map of int and double
typedef map <int, double> IDMAP;         // map of int and double
typedef map <double, int> FIMAP;         // double to integer map
typedef map <long, string> ISMAP;        // long to string map
typedef map <string, string> SSMAP;      // string to string map
typedef map <string, int> SIMAP;         // string to int map

class Atom: public Point {
public:
    int       _type;    // for van der wall radii
    string    _name;    // name of the atom
    double    _Bfactor; // store the temperature factor in pdb files
    int       _posn;    // position in the parent residue
    short int _state;   // This attribute can be used for multiple purposes
    //0: default
    //1: this atom is missing from native

    // common atom parameters
    static double radius[MAX_ATOM_TYPE];
    static double welldepth[MAX_ATOM_TYPE];
    static double s_volume[MAX_ATOM_TYPE];
    static double s_lambda[MAX_ATOM_TYPE];
    static double s_dgfree[MAX_ATOM_TYPE];
    static short acceptor[MAX_ATOM_TYPE];
    static short donor[MAX_ATOM_TYPE];
    static short hbondH[MAX_ATOM_TYPE];
    static double R_SC[NUM_RES_TYPE];

    explicit Atom(double X = 0., double Y = 0., double Z = 0., int I = UNDEF);
    explicit Atom(const Point& P);
    Atom(const Atom& A);
    void init(double X = 0., double Y = 0., double Z = 0., int I = UNDEF);
    Atom* operator=(const Atom& A);
    Atom* operator=(const Point& P);
    void CopyPos(const Atom& A);
    void reset();
    static void InitPar(const char*, const string&);
    void out() const;

    // @HACK Atoms placed at origin are considered "uninitialized"
    // and are ignored during potential energy computations
    inline int is_at_origin() const { return (x == 0. && y == 0. && z == 0.); }
    inline void move_to_origin() { x = y = z = 0.; }
};

#endif // _ATOM_
