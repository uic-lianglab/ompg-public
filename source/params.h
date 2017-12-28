#ifndef MDISGRO_PARAMS_H
#define MDISGRO_PARAMS_H

#include "mutation_manager.h"

#include <string>
#include <vector>

struct Params {
    //**********************/
    // Integer parameters
    //**********************/

    /**
     * Max number of conformations to attempt to generate
     * Usage: -n <positive integer>
     */
    int num_conf;
    /**
     * Number of completed samples to maintain in memory
     * Only the top conf_keep most energetically stable conformations
     * are maintained in memory at any point of execution.
     * Usage: -confkeep <positive integer>
     */
    int conf_keep;
    /**
     * Maximum number of successfully closed conformations to write to disk
     * This is a subset of the conformations that were kept in memory via conf_keep
     * Usage: -pdbout <positive integer>
     */
    int output_pdb;
    /**
     * For single loop/fragment library simulation.
     * Start residue identifier of loop to simulate
     * Usage: -start <positive integer>
     */
    int start;
    /**
     * For single loop/fragment library simulation.
     * End residue identifier of loop to simulate
     * Usage: -end <positive integer>
     */
    int end;
    /**
     * Use parameter fixed, random seed rather than random time-based seed
     * Usage: -rs <non-negative integer>
     */
    unsigned int rand_seed;
    /**
     * Torsion angle type - appears to have affect side chain growth (grow_type param for grow_sc())
     * @TODO - remove this - this probably should never be set and should always be 0
     *  or bad things will likely happen...
     * Usage: -angt <0: grow side chains consecutively from start to end region,
     *               1: grow from custom vector of possibly non-consecutive residue indices> 
     */
    int ang_type;
    /**
     * Number of backbone candidate distance states to sample per simulated residue
     * One will be selected stochastically based on potential energy
     * Usage: -nds <positive integer>
     */
    int num_dist_states;
    /**
     * Number of side chain candidate states to sample per simulated residue
     * One will be selected stochastically based on potential energy
     * Usage: -nscc <positive integer>
     */
    int num_sc_states;
    /**
     * A sample will be discarded if the backbone clash count exceeds this value
     * Usage: -MaxBBClashes <non-negative integer>
     */
    int max_num_bb_clashes;
    /**
     * A sample will be discarded if the side chain clash count exceeds this value
     * Usage: -MaxScClashes <non-negative integer>
     */
    int max_num_sc_clashes;
    /**
     * Maximum number of restarts allowed per multi-loop sample. If the total
     * count of loop restarts exceeds this value, the sample fails and a new
     * sample is started from scratch.
     * Usage: -murestart <non-negative integer>
     */
    int max_muloop_restart;
    /**
     * Maximum number of fragment libraries to check during a multiloop restart
     * Usage: -MaxFragChecks <non-negative integer>
     */
    int max_frag_checks;
    /**
     * If switch present, the program will attempt to refine side chain clashes for the input
     * template PDB over the parameter simulation region. No growth is attempted. This is
     * meant to be a post-processing utility for any generated samples that need further
     * side chain refinement due to exceeding a clash threshold. To set output path of refined
     * protein, set -pdb_outdir and/or -ProtName. Note: if iteration count is set to 0, then
     * this is simply prints the number of side chain clashes.
     * Usage: -PostScRefine <num_iterations >= 0>
     */
    int post_sc_refine;

    //**********************/
    // Boolean parameters
    //**********************/

    /**
     * Input PDB is checked for clashes according to overlap factor parameters.
     * The program then exits.
     * Usage: -ClashCheckOnly
     */
    bool clash_check_only;
    /**
     * Use analytic closure for last 3 loop residues
     * Usage: -close <0|1>
     */
    bool should_close;
    /**
     * Use fragment library during multi-loop growth
     * Not set via command line. Implicitly set when fragment library directory
     * directories and ini are specified.
     */
    bool use_multiloop_lib;
    /**
     * If switch present, then program will compute energy score of parameter pdb file (-f)
     * then exit program. If switch is missing from command line, then this value is false.
     * Usage: -score
     */
    bool score_only;
    /**
     * If switch is present, disable use of loop-specific rotamers. This may be useful if
     * sparsity of allowed (phi, psi) combinations prevents clash free loops from being generated.
     * If switch is not present, then rotamer libraries will be used for backbone growth.
     * Usage: -norot
     */
    bool use_rot_lib;

    //**********************/
    // Float parameters
    //**********************/

    /********************************************************************************
     * Overlap factors
     * If the ratio of distance between two atoms centers to the sum of their
     * atomic radii is less than the overlap factor (ofac) value, then the
     * atoms are colliding.
     *
     * Jacobson, Matthew P., David L.Pincus, Chaya S.Rapp, Tyler JF Day, Barry Honig,
     * David E.Shaw, and Richard A.Friesner.
     * "A hierarchical approach to all-atom protein loop prediction."
     * Proteins: Structure, Function, and Bioinformatics 55, no. 2 (2004) : 351 - 367.
     *********************************************************************************/

    /**
     * Overlap factor for collisions with adjacent residues
     * Usage: -ofa <positive real>
     */
    double ofa;
    /**
     * Overlap factor for collisions non-adjacent residues
     * Usage: -ofna <positive real>
     */
    double ofna;
    /**
     * Overlap factor for collisions with membrane atoms
     * Usage: -ofmemb <positive real>
     */
    double ofmemb;

    //**********************/
    // String parameters
    //**********************/

    /**
     * Path to input PDB file to simulate loop regions of
     * Usage: -f <path to PDB file>
     */
    std::string prot_file;
    /**
     * Protein name (PDB ID) implicitly extracted from PDB name
     * specified by prot_file. To overload explicitly, then:
     * Usage: -ProtName <StringIdNoExtension>
     */
    std::string prot_name;
    /**
     * Effective working directory for most parameter files.
     * Defaults to "./"
     * Usage: -dir <working directory path>
     */
    std::string dir;
    /**
     * Path to input parameters - in this case only thing set are
     * the energy coefficients.
     * @TODO - remove this option as it does nothing
     * Usage: -param|-ParameterFile <path to parameter file>
     */
    std::string param_file;
    /**
     * Path to text file containing begin and end residues indices for each
     * loop to be simulated. These regions are used to populate starts and
     * ends vectors.
     * Usage: -multiloops <path to text file>
     */
    std::string muloop_file;
    /**
     * All fragment libraries are assumed to reside in same folder
     * Usage: -multilooplib_dir <path to directory containing all fragment libraries>
     */
    std::string mulib_frag_dir;
    /**
     * Path to ini file specifying the names (i.e. identifiers) for the fragment
     * libraries to use for each simulated loop region. All fragment libraries
     * are assumed to reside in same folder as specified by -multilooplib_dir
     * Usage: -multilooplib_ini <path to ini file>
     */
    std::string mulib_frag_ini_path;
    /**
     * Used to avoid name clashes for output files from different runs of the
     * program.
     * Usage: -job_prefix <unique job identifier string>
     */
    std::string job_prefix;
    /**
     * The folder to write generated PDB structures to
     * Usage: -pdb_outdir <path to folder to write PDB data>
     */
    std::string pdb_outdir;
    /**
     * The folder to write captured statistics for this run of the
     * program (typically energy values for generated loops)
     * Usage: -stats_outdir <path to stats directory>
     */
    std::string stats_outdir;
    /**
     * Optional parameter to a membrane PDB file as generated by CHARMM-GUI.
     * The membrane is a lipid-bilayer in which the protein is embedded.
     * The membrane atoms are used for collision checking to make sure that
     * loops do not egregiously grow into the membrane region.
     * Usage: -memb <path to membrane PDB file>
     */
    std::string memb_file;
    /**
     * Path to text file containing being and end residues indices for each
     * mask region. These regions will ignored during energy and collision
     * calculations. This is to avoid any bias these regions may induce. In
     * particular, when generating fragment libraries for a future multi-
     * loop simulation, the other loop regions should likely be ignored.
     * Usage: -mask <path to text file>
     */
    std::string mask_file;

    //*****************************************************************/
    // Parameters for exporting a fragment library to CSV format used
    // by clustering tools
    //*****************************************************************/

    /**
     * All these parameters are set together by a single command line switch.
     * This utility will convert a PDB formatted fragment library to a CSV file format
     * usable by the clustering tools.
     * Usage: -exp_cl_csv_frags <path_to_input_pdb_frag_lib> <path_to_output_csv_frag_lib>
     */
    
    /**
     * Implicitly set to true when -exp_cl_csv_frags is encountered on command line
     */
    bool should_export_cl_csv_frag_lib;
    /**
     * The path to the input PBD formatted fragment library
     */
    std::string export_cl_csv_frag_lib_in_path;
    /**
     * The path to write the resulting CSV formatted library
     */
    std::string export_cl_csv_frag_lib_out_path;

    //*****************************************************************/
    // Parameters for exporting directories containing simulated loops
    // to CSV format used by clustering tools
    //*****************************************************************/

    /**
     * All these parameters are set together by a single command line switch.
     * This utility will convert a PDB formatted fragment library to a CSV file format
     * usable by the clustering tools.
     * Usage: -exp_cl_csv_pdbs <path_to_input_pdb_dir,|path_next_input_pdb_dir,|...> <path_to_output_csv_file> <-1|max_num_pdbs> <0|1 for backbone only>
     *
     * Example:
     * -exp_cl_csv_pdbs "./2iwv,./2iww" ./clust.csv -1 1 -start 217 -end 234
     *
     * Will merge the loop region(s) 217-234 in PDB folders .2iwv and .2iww into
     * a single CSV file. The resulting CSV file will be written to ./clust.csv.
     * Also, -1 signifies no cap on the number of PDBs contained in the CSV.
     * Finally, 1 signifies to only export the backbone atoms (N,Ca,Cb,C,O).
     *
     * Note: loop regions are specified via 'start', 'end' arguments or via
     * -multiloops <path to text file> arguments (just as in regular simulations)
     */

    /**
     * Implicitly set to true when -exp_cl_csv_pdbs is encountered on command line
     */
    bool should_export_cl_csv_pdbs;
    /**
     * Array of input directories each containing PDB files
     * Specified on command line by quoted csv list:
     * "<path_to_input_pdb_dir,|path_next_input_pdb_dir,|...>"
     * immediately following -exp_cl_csv_pdbs switch.
     * MUST NOT END WITH A PATH SEPARATOR (e.g "/" or "\\")
     */
    std::vector<std::string> export_cl_csv_pdb_in_dirs;
    /**
     * The path to write the resulting combined CSV formatted set
     */
    std::string export_cl_csv_pdb_out_path;
    /**
     * The path to write the list of files representing each (x,y,z) column set in CSV output
     */
    std::string export_cl_csv_pdb_out_ls_path;
    /**
     * The cap on the number of PDBs to export. -1 means no cap.
     */
    int export_cl_csv_pdb_max;
    /**
     * If true, only back bone atoms (N,Ca,Cb,C,O) for the loop
     * region of interest will exported to the merged CSV file.
     */
    bool export_cl_csv_pdb_bbone_only;

    //**********************/
    // Vector parameters
    //**********************/

    /**
     * Vector of start residues for each simulated loop region.
     * Eventually populated via data present within muloop_file
     */
    std::vector<int> starts;
    /**
     * Vector of end residues for each simulated loop region
     * Eventually populated via data present within muloop_file
     */
    std::vector<int> ends;
    /**
     * Vector of start residues for each masked loop region.
     * Eventually populated via data present within mask_file
     */
    std::vector<int> starts_mask;
    /**
     * Vector of end residues for each masked loop region
     * Eventually populated via data present within mask_file
    */
    std::vector<int> ends_mask;

    //**********************/
    // Mutations
    //**********************/

    /**
     * Manages parsing mutation parameters
     */
    MutationManager muts;

    /**
     * Default constructor
     */
    Params();

    /**
     * Initializes parameters from command line
     */
    void parse_command_line(const int argc, const char **argv);

    /**
     * Print parameters
     */
    void print() const;
};

#endif //MDISGRO_PARAMS_H
