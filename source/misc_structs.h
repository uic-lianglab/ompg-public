/**
 * Header for small POD-type structures to avoid inclusion of
 * large bulky class definitions
 */

#ifndef MiscStructs_h
#define MiscStructs_h

/**
 * Minimalist structure for passing around atom information
 */
struct atom_info {
    int res_ix;
    int atom_ix;
    short res_type;
};

/**
 * Information about a simulation loop region
 */
struct loop_info {
    int Start;
    int End;
    int CurrentPos;
};

#endif // MiscStructs_h
