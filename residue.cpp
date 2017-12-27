/*
residue.cpp
*/

#include "residue.h"
#include "util.h"

#include <iostream>

string Residue::Name1[25];  // one letter name including bases
string Residue::Name3[25];  // three letter name
SIMAP Residue::AIMap;  // Amino acid to integer map
SIMAP Residue::SIMap;  // Secondary structure to integer map
SSMAP Residue::AAMap;  // amino acid three letter name and one letter name map
SIMAP Residue::AtomMap;
SIMAP Residue::AtomIndexMap;
ISMAP Residue::ResMap;

string Residue::cType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
int Residue::vdwType[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
int Residue::prev_atom[NUM_RES_TYPE][MAX_NUM_ATOM_RES][3];
// Bond lengths and bond angles can be found:
// Engh, Richard A., and Robert Huber. "Accurate bond and angle parameters for X-ray protein structure refinement." Acta Crystallographica Section A : Foundations of Crystallography 47.4 (1991) : 392 - 400.
// Some helpful tables which cite above resource:
// http://www.ruppweb.org/Xray/tutorial/engtabl.html
// http://www.ruppweb.org/Xray/tutorial/protein_structure.htm
// Also:
// http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
double Residue::bond_length[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
double Residue::bond_angle[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
double Residue::torsion[NUM_RES_TYPE][MAX_NUM_ATOM_RES];
double Residue::size[NUM_RES_TYPE];
double Residue::sc_size[NUM_RES_TYPE];
double Residue::bb_size;
int Residue::numAtom[NUM_RES_TYPE];
IVEC Residue::funcAtom[NUM_RES_TYPE];

Residue::Residue(int n) {
    init(n);
}

Residue::Residue(const Residue &R) {
    init(MAX_NUM_ATOM_RES);
    *this = R;
}

void Residue::init(int n) {
    _atom = new Atom[n];
    for (int i = 0; i < n; i++) {
        _atom[i]._posn = i;
    }
    _type = -1;
    _numAtom = n;
    _scState = -1;
    for (int i = 0; i < 5; i++)
        _scChi[i] = 999.99;
}

Residue *Residue::operator=(const Residue &R) {
// Copy the atom's from this residue
    _posn = R._posn;
    _type = R._type;
    _numAtom = R._numAtom;
    _phi = R._phi;
    _psi = R._psi;
    _omega = R._omega;
    _scState = R._scState;
    _pdbIndex = R._pdbIndex;
    _bbc = R._bbc;
    _scc = R._scc;
    _center = R._center;
    for (int i = 0; i < R._numAtom; i++) {
        _atom[i] = R._atom[i];
    }
    for (int i = R._numAtom; i < MAX_NUM_ATOM_RES; i++)
        _atom[i].init();

    for (int i = 0; i < 5; i++)
        _scChi[i] = R._scChi[i];
    return this;
}

Residue::~Residue() {
    Destruct();
}

void Residue::Destruct() {
    delete[] _atom;
}

void Residue::InitMap() {
    AAMap.insert(SSMAP::value_type("Ala", "A"));
    AAMap.insert(SSMAP::value_type("Cys", "C"));
    AAMap.insert(SSMAP::value_type("Asp", "D"));
    AAMap.insert(SSMAP::value_type("Glu", "E"));
    AAMap.insert(SSMAP::value_type("Phe", "F"));
    AAMap.insert(SSMAP::value_type("Gly", "G"));
    AAMap.insert(SSMAP::value_type("His", "H"));
    AAMap.insert(SSMAP::value_type("Ile", "I"));
    AAMap.insert(SSMAP::value_type("Lys", "K"));
    AAMap.insert(SSMAP::value_type("Leu", "L"));
    AAMap.insert(SSMAP::value_type("Met", "M"));
    AAMap.insert(SSMAP::value_type("Asn", "N"));
    AAMap.insert(SSMAP::value_type("Pro", "P"));
    AAMap.insert(SSMAP::value_type("Gln", "Q"));
    AAMap.insert(SSMAP::value_type("Arg", "R"));
    AAMap.insert(SSMAP::value_type("Ser", "S"));
    AAMap.insert(SSMAP::value_type("Thr", "T"));
    AAMap.insert(SSMAP::value_type("Val", "V"));
    AAMap.insert(SSMAP::value_type("Trp", "W"));
    AAMap.insert(SSMAP::value_type("Tyr", "Y"));
    AAMap.insert(SSMAP::value_type("ALA", "A"));
    AAMap.insert(SSMAP::value_type("CYS", "C"));
    AAMap.insert(SSMAP::value_type("ASP", "D"));
    AAMap.insert(SSMAP::value_type("GLU", "E"));
    AAMap.insert(SSMAP::value_type("PHE", "F"));
    AAMap.insert(SSMAP::value_type("GLY", "G"));
    AAMap.insert(SSMAP::value_type("HIS", "H"));
    AAMap.insert(SSMAP::value_type("ILE", "I"));
    AAMap.insert(SSMAP::value_type("LYS", "K"));
    AAMap.insert(SSMAP::value_type("LEU", "L"));
    AAMap.insert(SSMAP::value_type("MET", "M"));
    AAMap.insert(SSMAP::value_type("ASN", "N"));
    AAMap.insert(SSMAP::value_type("PRO", "P"));
    AAMap.insert(SSMAP::value_type("GLN", "Q"));
    AAMap.insert(SSMAP::value_type("ARG", "R"));
    AAMap.insert(SSMAP::value_type("SER", "S"));
    AAMap.insert(SSMAP::value_type("THR", "T"));
    AAMap.insert(SSMAP::value_type("VAL", "V"));
    AAMap.insert(SSMAP::value_type("TRP", "W"));
    AAMap.insert(SSMAP::value_type("TYR", "Y"));
    AAMap.insert(SSMAP::value_type("UNK", "U"));
    AAMap.insert(SSMAP::value_type("  C", "c"));
    AAMap.insert(SSMAP::value_type("  T", "t"));
    AAMap.insert(SSMAP::value_type("  A", "a"));
    AAMap.insert(SSMAP::value_type("  G", "g"));
    AAMap.insert(SSMAP::value_type("  U", "u"));

    AAMap.insert(SSMAP::value_type("A", "ALA"));
    AAMap.insert(SSMAP::value_type("C", "CYS"));
    AAMap.insert(SSMAP::value_type("D", "ASP"));
    AAMap.insert(SSMAP::value_type("E", "GLU"));
    AAMap.insert(SSMAP::value_type("F", "PHE"));
    AAMap.insert(SSMAP::value_type("G", "GLY"));
    AAMap.insert(SSMAP::value_type("H", "HIS"));
    AAMap.insert(SSMAP::value_type("I", "ILE"));
    AAMap.insert(SSMAP::value_type("K", "LYS"));
    AAMap.insert(SSMAP::value_type("L", "LEU"));
    AAMap.insert(SSMAP::value_type("M", "MET"));
    AAMap.insert(SSMAP::value_type("N", "ASN"));
    AAMap.insert(SSMAP::value_type("P", "PRO"));
    AAMap.insert(SSMAP::value_type("Q", "GLN"));
    AAMap.insert(SSMAP::value_type("R", "ARG"));
    AAMap.insert(SSMAP::value_type("S", "SER"));
    AAMap.insert(SSMAP::value_type("T", "THR"));
    AAMap.insert(SSMAP::value_type("V", "VAL"));
    AAMap.insert(SSMAP::value_type("W", "TRP"));
    AAMap.insert(SSMAP::value_type("Y", "TYR"));
    AAMap.insert(SSMAP::value_type("U", "UNK"));
    AAMap.insert(SSMAP::value_type("c", "  C"));
    AAMap.insert(SSMAP::value_type("t", "  T"));
    AAMap.insert(SSMAP::value_type("a", "  A"));
    AAMap.insert(SSMAP::value_type("g", "  G"));
    AAMap.insert(SSMAP::value_type("u", "  U"));

    AAMap.insert(SSMAP::value_type("Alanine", "A"));
    AAMap.insert(SSMAP::value_type("Cysteine", "C"));
    AAMap.insert(SSMAP::value_type("Aspartate", "D"));
    AAMap.insert(SSMAP::value_type("Glutamate", "E"));
    AAMap.insert(SSMAP::value_type("Phenylalanine", "F"));
    AAMap.insert(SSMAP::value_type("Glycine", "G"));
    AAMap.insert(SSMAP::value_type("Histidine", "H"));
    AAMap.insert(SSMAP::value_type("Isoleucine", "I"));
    AAMap.insert(SSMAP::value_type("Lysine", "K"));
    AAMap.insert(SSMAP::value_type("Leucine", "L"));
    AAMap.insert(SSMAP::value_type("Methionine", "M"));
    AAMap.insert(SSMAP::value_type("Asparagine", "N"));
    AAMap.insert(SSMAP::value_type("Proline", "P"));
    AAMap.insert(SSMAP::value_type("Glutamine", "Q"));
    AAMap.insert(SSMAP::value_type("Arginine", "R"));
    AAMap.insert(SSMAP::value_type("Serine", "S"));
    AAMap.insert(SSMAP::value_type("Threonine", "T"));
    AAMap.insert(SSMAP::value_type("Valine", "V"));
    AAMap.insert(SSMAP::value_type("Tryptophan", "W"));
    AAMap.insert(SSMAP::value_type("Tyrosine", "Y"));
    AAMap.insert(SSMAP::value_type("Disulfide", "Z"));
    AAMap.insert(SSMAP::value_type("Unknown", "U"));
    AAMap.insert(SSMAP::value_type("Cytosine", "c"));
    AAMap.insert(SSMAP::value_type("Thymine", "t"));
    AAMap.insert(SSMAP::value_type("Adenine", "a"));
    AAMap.insert(SSMAP::value_type("Guanine", "g"));
    AAMap.insert(SSMAP::value_type("Uracil", "u"));

    // AIMap is the map for amino acid name and their integer representation
    AIMap.insert(SIMAP::value_type("A", 0));
    AIMap.insert(SIMAP::value_type("C", 1));
    AIMap.insert(SIMAP::value_type("D", 2));
    AIMap.insert(SIMAP::value_type("E", 3));
    AIMap.insert(SIMAP::value_type("F", 4));
    AIMap.insert(SIMAP::value_type("G", 5));
    AIMap.insert(SIMAP::value_type("H", 6));
    AIMap.insert(SIMAP::value_type("I", 7));
    AIMap.insert(SIMAP::value_type("K", 8));
    AIMap.insert(SIMAP::value_type("L", 9));
    AIMap.insert(SIMAP::value_type("M", 10));
    AIMap.insert(SIMAP::value_type("N", 11));
    AIMap.insert(SIMAP::value_type("P", 12));
    AIMap.insert(SIMAP::value_type("Q", 13));
    AIMap.insert(SIMAP::value_type("R", 14));
    AIMap.insert(SIMAP::value_type("S", 15));
    AIMap.insert(SIMAP::value_type("T", 16));
    AIMap.insert(SIMAP::value_type("V", 17));
    AIMap.insert(SIMAP::value_type("W", 18));
    AIMap.insert(SIMAP::value_type("Y", 19));
    AIMap.insert(SIMAP::value_type("a", 20));
    AIMap.insert(SIMAP::value_type("c", 21));
    AIMap.insert(SIMAP::value_type("g", 22));
    AIMap.insert(SIMAP::value_type("t", 23));
    AIMap.insert(SIMAP::value_type("u", 24));
    AIMap.insert(SIMAP::value_type("Z", 25));
    AIMap.insert(SIMAP::value_type("U", 26));

    AIMap.insert(SIMAP::value_type("ALA", 0));
    AIMap.insert(SIMAP::value_type("CYS", 1));
    AIMap.insert(SIMAP::value_type("ASP", 2));
    AIMap.insert(SIMAP::value_type("GLU", 3));
    AIMap.insert(SIMAP::value_type("PHE", 4));
    AIMap.insert(SIMAP::value_type("GLY", 5));
    AIMap.insert(SIMAP::value_type("HIS", 6));
    AIMap.insert(SIMAP::value_type("ILE", 7));
    AIMap.insert(SIMAP::value_type("LYS", 8));
    AIMap.insert(SIMAP::value_type("LEU", 9));
    AIMap.insert(SIMAP::value_type("MET", 10));
    AIMap.insert(SIMAP::value_type("ASN", 11));
    AIMap.insert(SIMAP::value_type("PRO", 12));
    AIMap.insert(SIMAP::value_type("GLN", 13));
    AIMap.insert(SIMAP::value_type("ARG", 14));
    AIMap.insert(SIMAP::value_type("SER", 15));
    AIMap.insert(SIMAP::value_type("THR", 16));
    AIMap.insert(SIMAP::value_type("VAL", 17));
    AIMap.insert(SIMAP::value_type("TRP", 18));
    AIMap.insert(SIMAP::value_type("TYR", 19));
    AIMap.insert(SIMAP::value_type("  A", 20));
    AIMap.insert(SIMAP::value_type("  C", 21));
    AIMap.insert(SIMAP::value_type("  G", 22));
    AIMap.insert(SIMAP::value_type("  T", 23));
    AIMap.insert(SIMAP::value_type("  U", 24));

    // SIMap is the map for secondary structure letters and their integers
    SIMap.insert(SIMAP::value_type("H", 0));
    SIMap.insert(SIMAP::value_type("E", 1));
    SIMap.insert(SIMAP::value_type("C", 2));

    //ISMap is the map for integers and 3-letters residue type

    ResMap.insert(ISMAP::value_type(0, "ALA"));
    ResMap.insert(ISMAP::value_type(1, "CYS"));
    ResMap.insert(ISMAP::value_type(2, "ASP"));
    ResMap.insert(ISMAP::value_type(3, "GLU"));
    ResMap.insert(ISMAP::value_type(4, "PHE"));
    ResMap.insert(ISMAP::value_type(5, "GLY"));
    ResMap.insert(ISMAP::value_type(6, "HIS"));
    ResMap.insert(ISMAP::value_type(7, "ILE"));
    ResMap.insert(ISMAP::value_type(8, "LYS"));
    ResMap.insert(ISMAP::value_type(9, "LEU"));
    ResMap.insert(ISMAP::value_type(10, "MET"));
    ResMap.insert(ISMAP::value_type(11, "ASN"));
    ResMap.insert(ISMAP::value_type(12, "PRO"));
    ResMap.insert(ISMAP::value_type(13, "GLN"));
    ResMap.insert(ISMAP::value_type(14, "ARG"));
    ResMap.insert(ISMAP::value_type(15, "SER"));
    ResMap.insert(ISMAP::value_type(16, "THR"));
    ResMap.insert(ISMAP::value_type(17, "VAL"));
    ResMap.insert(ISMAP::value_type(18, "TRP"));
    ResMap.insert(ISMAP::value_type(19, "TYR"));
    ResMap.insert(ISMAP::value_type(20, "  A"));
    ResMap.insert(ISMAP::value_type(21, "  C"));
    ResMap.insert(ISMAP::value_type(22, "  G"));
    ResMap.insert(ISMAP::value_type(23, "  T"));
    ResMap.insert(ISMAP::value_type(24, "  U"));
}

void Residue::InitPar(char *parFile, string &dir) {
    int i, j, tmpInt, numLine;
    ifstream inFile;
    string tmpStr;
    char tempCh[200], tmpLine[500];
    vector<string> strVec;       // used in reading files
    SSMAP::iterator aaItr;
    SIMAP::iterator siItr;

    inFile.open(parFile, ios::in);
    if (!inFile.is_open()) {
        strcpy(tempCh, dir.c_str());
        strcat(tempCh, parFile);
        inFile.open(tempCh, ios::in);
        if (!inFile.is_open()) {
            cout << "cannot open parameter file " << tempCh << endl;
            exit(0);
        }
    }

    // read parameters
    inFile.clear();
    inFile.seekg(0, ios::beg);
    while (!inFile.eof()) {
        inFile.getline(tmpLine, 1000);
        sscanf(tmpLine, "%*s %s %d", tempCh, &numLine);
        if (strcmp(tempCh, "restype") == 0) { // read residue types
            int AtomIndex = 0;
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                Residue::numAtom[i] = (split(tmpStr, ' ', strVec) - 2) / 2;
                Name1[i] = strVec[0].c_str();
                tmpStr = Name1[i];
                if ((aaItr = AAMap.find(tmpStr)) == AAMap.end()) {
                    cout << "file reading error, wrong AA type " << tmpStr << endl;
                    exit(0);
                }
                else {
                    Name3[i] = aaItr->second;
                }
                tmpInt = 2;
                for (j = 0; j < Residue::numAtom[i]; ++j) {
                    Residue::cType[i][j] = strVec[tmpInt++];
                    Residue::vdwType[i][j] = atoi(strVec[tmpInt++].c_str());
                    tmpStr = Name1[i];
                    tmpStr += Residue::cType[i][j];
                    Residue::AtomMap.insert(SIMAP::value_type(tmpStr, j));
                    Residue::AtomIndexMap.insert(SIMAP::value_type(tmpStr, AtomIndex++));
                }
            }
        }
        else if (strcmp(tempCh, "atomDependency") == 0) {
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                split(tmpStr, ' ', strVec);
                if (strVec[1] != Residue::Name1[i]) {
                    cout << "atomDependency error, wrong AA type " << i << " " << strVec[1][0] << " " <<
                    Residue::Name1[i] << endl;
                    exit(0);
                }
                // read the atom dependency data for each atom
                tmpInt = 2;
                for (j = 0; j < Residue::numAtom[i]; ++j) {
                    if (strVec[tmpInt++] != Residue::cType[i][j]) {
                        cout << "atomDependency error, wrong atom type " << strVec[tmpInt] << " " <<
                        Residue::cType[i][j] << " " << i << endl;
                        exit(0);
                    }
                    // read the next three integers
                    Residue::prev_atom[i][j][0] = atoi(strVec[tmpInt++].c_str());
                    Residue::prev_atom[i][j][1] = atoi(strVec[tmpInt++].c_str());
                    Residue::prev_atom[i][j][2] = atoi(strVec[tmpInt++].c_str());
                }
            }
        }
        else if (strcmp(tempCh, "bondParameter") == 0) {
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                split(tmpStr, ' ', strVec);
                if (strVec[1] != Name1[i]) {
                    cout << "bondParameter error, wrong AA type " << i << " " << strVec[1][0] << endl;
                    exit(0);
                }
                // read the bond parameters
                tmpInt = 2;
                for (j = 0; j < numAtom[i]; ++j) {
                    if (strVec[tmpInt++] != cType[i][j]) {
                        cout << "bondParamter error, wrong atom type " << tmpInt << " " << strVec[tmpInt - 1] << " " <<
                        cType[i][j] << endl;
                        exit(0);
                    }
                    // read the next two double numbers
                    bond_length[i][j] = atof(strVec[tmpInt++].c_str());
                    bond_angle[i][j] = (PI / 180) * atof(strVec[tmpInt++].c_str());
                }
            }
        }
        else if (strcmp(tempCh, "residueSize") == 0) {
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                split(tmpStr, ' ', strVec);
                size[i] = atof(strVec[2].c_str());
                sc_size[i] = atof(strVec[4].c_str());
            }
        }
        else if (strcmp(tempCh, "torsionAngle") == 0) {
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                split(tmpStr, ' ', strVec);
                tmpInt = 2;
                for (j = 0; j < numAtom[i]; ++j) {
                    if (strVec[tmpInt++] != cType[i][j]) {
                        std::cout << "torsionAngle error, wrong atom type " << strVec[tmpInt - 1] << " " << cType[i][j] <<
                        endl;
                        exit(0);
                    }
                    if (atof(strVec[tmpInt].c_str()) == -1234) {
                        torsion[i][j] = -1234;
                        tmpInt++;
                    }
                    else
                        torsion[i][j] = (PI / 180) * atof(strVec[tmpInt++].c_str());
                }
            }
        }
        else if (strcmp(tempCh, "functionalAtoms") == 0) {
            for (i = 0; i < numLine; ++i) {
                inFile.getline(tmpLine, 1000);
                tmpStr = tmpLine;
                strVec.clear();
                split(tmpStr, ' ', strVec);
                tmpInt = 2;
                for (j = 2; j < strVec.size(); ++j) {
                    tmpStr = strVec[0];
                    tmpStr = tmpStr + strVec[j];
                    siItr = AtomMap.find(tmpStr);
                    if (siItr == AtomMap.end()) {
                        cout << "atom type error in reading functional atoms! " << i << " " << j << " " << tmpStr <<
                        endl;
                        exit(0);
                    }
                    tmpInt = siItr->second;
                    funcAtom[i].push_back(tmpInt);
                }
            }
        }
        else {
            continue;
        }
    }

    // set backbone size
    bb_size = 1.93;
    // initialize functional group atoms

}

// calculate side chain center
void Residue::cal_scc() {
    Atom tmpAtm;
    for (int i = NUM_BB_ATOM; i < _numAtom; ++i) {
        tmpAtm = tmpAtm + _atom[i];
    }
    tmpAtm = tmpAtm / (_numAtom - NUM_BB_ATOM);
    _scc = tmpAtm;
}

void Residue::out() {
    cout << _type << " " << Name1[_type] << " " << _numAtom << ": ";
    for (int i = 0; i < _numAtom; ++i) {
        cout << i << " " << _atom[i]._type << ", ";
    }
    cout << endl;
    cout << "(" << _atom[0].x << ", " << _atom[0].y << ", " << _atom[0].z
    << ") ";
}

// debug function to print atom x,y,z coordinates
void Residue::debug_print() const {
    std::cout << "--------------------------\n";
    std::cout << "RESIDUE STATE:\n";
    std::cout << "\tname: " << this->get_name_3() << "\n";
    std::cout << "\tposition: " << this->_posn << "\n";
    std::cout << "\tnum atoms: " << this->_numAtom << "\n";

    for (int i=0; i<this->_numAtom; ++i) {
        const Atom &atom = this->_atom[i];
        std::cout << "\t\t" << i << ": "
                  << atom._name << ": "
                  << atom._type <<": ("
                  << atom.x << ", " << atom.y << ", " << atom.z << ")\n";
    }
}
