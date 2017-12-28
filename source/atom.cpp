/* atom.cpp
*/

#include "atom.h"

#include <iostream>

// Atom class//
double Atom::radius[MAX_ATOM_TYPE];
double Atom::welldepth[MAX_ATOM_TYPE];
double Atom::s_volume[MAX_ATOM_TYPE];
double Atom::s_lambda[MAX_ATOM_TYPE];
double Atom::s_dgfree[MAX_ATOM_TYPE];
short Atom::acceptor[MAX_ATOM_TYPE];
short Atom::donor[MAX_ATOM_TYPE];
short Atom::hbondH[MAX_ATOM_TYPE];
double Atom::R_SC[NUM_RES_TYPE] = {1.53, 2.14, 2.49, 3.16, 3.42, 0.0, 3.17, 2.35, 3.57, 2.65, 3.00, 2.52, 1.88, 3.15, 4.17, 1.95, 1.95, 1.97, 3.89, 3.78, 0,0,0,0,0,0,0,0,0,0};  // new assignment calculated recently

// Simple initialization
Atom::Atom(double X, double Y, double Z, int I) {
    init(X, Y, Z, I);
}
Atom::Atom(const Point& P) {
    init(P.x, P.y, P.z);
}
Atom::Atom(const Atom& A) {
    init();
    *this = A;
}
void Atom::init(double X, double Y, double Z, int I) {
    x = X;
    y = Y;
    z = Z;
    _type = I;
    _Bfactor = -1;
    _state = 0;
}

// This is the full assignment operator, copying all members
Atom* Atom::operator=(const Atom& A) {
    x = A.x;
    y = A.y;
    z = A.z;
    _type = A._type;
    _name = A._name;
    _Bfactor = A._Bfactor;
    _posn = A._posn;
    _state = A._state;
    return this;
}
Atom* Atom::operator=(const Point& P) {
    x = P.x;
    y = P.y;
    z = P.z;
    return this;
}

// This is a limited assignment operator, copying only the spatial position
void Atom::CopyPos(const Atom& A) {
    x = A.x;
    y = A.y;
    z = A.z;
}

void Atom::reset() {
    init();
}

void Atom::out() const {
    std::cout << "(" << x << ", " << y << ", " << z << ", " << _type << ") ";
}

void Atom::InitPar(const char* parFile, const string& dir) {
    int i,numLine;
    ifstream inFile;
    char tempCh[200],tmpLine[500];
    cout<<"Reading Atom and Residue Parameters in "<<parFile<<endl;
    inFile.open(parFile,ios::in);
    if(!inFile.is_open()) {
        strcpy(tempCh,dir.c_str());
        strcat(tempCh,parFile);
        inFile.open(tempCh,ios::in);
        if(!inFile.is_open()) {
            cout<<"cannot open parameter file "<<tempCh<<endl;
            exit(0);
        }
    }

    // read parameters on atoms
    inFile.clear();
    inFile.seekg(0,ios::beg);
    while(!inFile.eof()) {
        inFile.getline(tmpLine,1000);
        sscanf(tmpLine,"%*s %s %d",tempCh,&numLine);
        if(strcmp(tempCh,"atom_parameters")==0) { // read vdw parameters
            for(i=0; i<numLine; ++i) {
                inFile.getline(tmpLine,1000);
                sscanf(tmpLine, "%lf, %lf, %lf, %lf, %lf, %d %d %d",
                       &Atom::radius[i+1],
                       &Atom::welldepth[i+1],
                       &Atom::s_volume[i+1],
                       &Atom::s_lambda[i+1],
                       &Atom::s_dgfree[i+1],
                       (int*)&Atom::acceptor[i+1],
                       (int*)&Atom::donor[i+1],
                       (int*)&Atom::hbondH[i+1]);
            }
        }
        else continue;
    }
    inFile.close();
}
