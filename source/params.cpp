//
// Created by apr on 10/22/15.
//

#include "params.h"
#include "util.h"

void Usage() {
    cout <<
    "Pretzel: program for simulation of protein structures using multi-loop Distance-guided chain-Growth Monte Carlo method." <<
    endl;
    cout << "options:\n";
    cout << "\t-f:\t input coordinate file of the protein" << endl;
    cout << "\t-start:\t The starting position for growth mode. Default is -1, starting from the beginning." << endl;
    cout << "\t-end:\t The end position for growth mode. Default is -1, ending at the last residue." << endl;
    cout << "\t-n:\t number of conformations" << endl;
    cout << "\t-pdbout: \t the number of the output conformations." << endl;
    cout << "\t-nds:\t number of sampling states" << endl;
    cout << "\t-nscc:\t number of side chain states" << endl;
    cout << "\t-multiloops:\t The file recording the start and end residues of multiple loops." << endl;
    cout << "\t" << endl;
    cout <<
    "Fragments generation example:\n\t./pretzel -f MultiLoop/1UTK_36.pdb -n 10000 -nds 32 -start xx -end xx -eval -confkeep 1000 -pdbout 100" <<
    endl;
    cout <<
    "Multi-loop example:\n\t./pretzel -f MultiLoop/1UTK_36.pdb -n 10000 -nds 32 -multiloops MultiLoop/1UTK_36_mulist -confkeep 500 -nscc 20 -multilooplib pdb_output -pdbout 20" <<
    endl;
    cout << "\t \n";
    cout << "\t \n";
    exit(0);
}

/**
 * Default constructor
 */
Params::Params() :
        // Integers
        num_conf(1000),
        conf_keep(1),
        output_pdb(1),
        start(-1),
        end(-1),
        // @TODO - add 64-bit support through random library
        rand_seed(static_cast<unsigned int>(time(NULL))),
        ang_type(2),
        num_dist_states(32),
        num_sc_states(16),
        max_num_bb_clashes(0),
        max_num_sc_clashes(15),
        max_muloop_restart(25),
        max_frag_checks(50),
        post_sc_refine(-1),
        // Booleans
        clash_check_only(false),
        should_close(true),
        use_multiloop_lib(false),
        score_only(false),
        use_rot_lib(true),
        // Floats
        ofa(0.80),
        ofna(0.90),
        ofmemb(0.85),
        // Strings
        prot_file(""),
        prot_name(""),
        dir("./"),
        param_file("data/parameter.txt"),
        muloop_file(""),
        mulib_frag_dir(""),
        mulib_frag_ini_path(""),
        job_prefix("topconf"),
        pdb_outdir("./pdb_output"),
        stats_outdir("./stats"),
        memb_file(""),
        mask_file(""),
        // Export for frag lib clustering
        should_export_cl_csv_frag_lib(false),
        export_cl_csv_frag_lib_in_path(""),
        export_cl_csv_frag_lib_out_path(""),
        // Export for PDB clustering
        should_export_cl_csv_pdbs(false),
        export_cl_csv_pdb_out_path(""),
        export_cl_csv_pdb_out_ls_path(""),
        export_cl_csv_pdb_max(-1),
        export_cl_csv_pdb_bbone_only(1)
        {}

/**
 * Initializes parameters from command line
 */
void Params::parse_command_line(const int argc, const char *argv[]) {

    // There are no arguments, so something went awry
    if (argc == 1)
        Usage();

    const int MAX_NUM_CONF = 5000000;

    // Get the mode
    int ii = 1;
    while (ii < argc) {

        //**********************/
        // Integer parameters
        //**********************/

        if (strcmp(argv[ii], "-n") == 0) {
            // sample size
            this->num_conf = atoi(argv[ii + 1]);
            if (this->num_conf > MAX_NUM_CONF) {
                cout << "The maximum sample size is " << MAX_NUM_CONF << "." << endl;
                cout << "Please decrease the input sample size (number of confs)"
                << "or modify the program." << endl;
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-confkeep") == 0) {
            this->conf_keep = atoi(argv[ii + 1]);
            if (this->conf_keep <= 0) {
                cout << "Error number of memory conformations to keep -confkeep must be positive.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-pdbout") == 0) {
            // Number of conformations to write to disk
            this->output_pdb = atoi(argv[ii + 1]);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-start") == 0) {
            this->start = atoi(argv[ii + 1]);
            this->starts.push_back(this->start);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-end") == 0) {
            this->end = atoi(argv[ii + 1]);
            this->ends.push_back(this->end);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-rs") == 0) {
            // random seed
            this->rand_seed = strtoul(argv[ii + 1], NULL, 0);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-angt") == 0) {
            // torsion angle representation type
            this->ang_type = atoi(argv[ii + 1]);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-nds") == 0) {
            this->num_dist_states = atoi(argv[ii + 1]);
            // number of backbone distance states
            if (this->num_dist_states <= 0) {
                cout << "Error num distance states -nds must be positive.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-nscc") == 0) {
            // number of side chain states in side chain sampling
            this->num_sc_states = atoi(argv[ii + 1]);
            if (this->num_sc_states <= 0) {
                cout << "Error num side chain states -nscc must be positive.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-MaxBBClashes") == 0) {
            // maximum number of backbone clashes allowed
            this->max_num_bb_clashes = atoi(argv[ii + 1]);
            if (this->max_num_bb_clashes < 0) {
                cout << "Error: MaxBBClashes must be non-negative.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-MaxSCClashes") == 0) {
            // maximum number of side chain clashes allowed
            this->max_num_sc_clashes = atoi(argv[ii + 1]);
            if (this->max_num_sc_clashes < 0) {
                cout << "Error: MaxSCClashes must be non-negative.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-murestart") == 0) {
            // maximum number of multi-loop restarts allowed
            this->max_muloop_restart = atoi(argv[ii + 1]);
            if (this->max_muloop_restart < 0) {
                cout << "Error: murestart must be non-negative.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-MaxFragChecks") == 0) {
            // maximum number of fragments to check when restarting a loop
            // during multi-loop growth
            this->max_frag_checks = atoi(argv[ii + 1]);
            if (this->max_frag_checks < 0) {
                cout << "Error: MaxFragChecks must be non-negative.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-PostScRefine") == 0) {
            // number of side chain refinement steps for post-processing
            // utility only
            this->post_sc_refine = atoi(argv[ii + 1]);
            ii += 2;
        }

        //**********************/
        // Boolean parameters
        //**********************/
        
        else if (strcmp(argv[ii], "-ClashCheckOnly") == 0) {
            // check parameter PDB for clashes then exit
            this->clash_check_only = true;
            ++ii;
        }
        else if (strcmp(argv[ii], "-close") == 0) {
            // use analytic closure in SMC fragment regrowth
            this->should_close = atoi(argv[ii + 1]) != 0;
            ii += 2;
        }
        else if (strcmp(argv[ii], "-score") == 0) {
            // score parameter pdb then exit program
            this->score_only = true;
            ++ii;
        }
        else if (strcmp(argv[ii], "-norot") == 0) {
            // disable use of loop specific rotamers
            // - this may be useful to explore phi,psi combinations
            // that are more clash free
            this->use_rot_lib = false;
            ++ii;
        }

        //**********************/
        // Float parameters
        //**********************/

        else if (strcmp(argv[ii], "-ofa") == 0) {
            // overlap factor for atoms within adjacent residues
            this->ofa = atof(argv[ii + 1]);
            if (this->ofa <= 0.0) {
                cout << "Error: adjacent residue collision overlap factor (read value="
                    << this->ofa << ") must be positive.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-ofna") == 0) {
            // overlap factor for atoms within non-adjacent residues
            this->ofna = atof(argv[ii + 1]);
            if (this->ofna <= 0.0) {
                cout << "Error: non-adjacent residue collision overlap factor (read value="
                    << this->ofna << ") must be positive.\n";
                exit(0);
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-ofmemb") == 0) {
            // overlap factor for membrane atoms
            this->ofmemb = atof(argv[ii + 1]);
            if (this->ofmemb <= 0.0) {
                cout << "Error: membrane bilayer collision overlap factor (read value="
                    << this->ofmemb << ") must be positive.\n";
                exit(0);
            }
            ii += 2;
        }

        //**********************/
        // String parameters
        //**********************/

        else if (strcmp(argv[ii], "-f") == 0) {
            // protein coordinate file
            this->prot_file = argv[ii + 1];
            if (this->prot_name == "")
                this->prot_name = File2ProtName(this->prot_file);
            ii += 2;
        }
        else if (strcmp(argv[ii], "-ProtName") == 0) {
            // overload protein name
            this->prot_name = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-e") == 0) {
            if (strcmp(argv[ii + 1], "loodis") == 0)
                PF::cal[EM_LOODIS] = true;
            else {
                cout << "Unrecognized energy term: " << argv[ii + 1] << " !" << endl;
            }
            ii += 2;
        }
        else if (strcmp(argv[ii], "-dir") == 0) {
            this->dir = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-ParameterFile") == 0 ||
                 strcmp(argv[ii], "-param") == 0) {
            this->param_file = argv[ii + 1];
            cout << "parameter file: " << this->param_file << endl;
            ii += 2;
        }
        else if (strcmp(argv[ii], "-multiloops") == 0) {
            // multiple loops growing
            this->muloop_file = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-multilooplib_dir") == 0) {
            this->use_multiloop_lib = true;
            this->mulib_frag_dir = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-multilooplib_ini") == 0) {
            this->use_multiloop_lib = true;
            this->mulib_frag_ini_path = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-job_prefix") == 0) {
            this->job_prefix = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-pdb_outdir") == 0) {
            this->pdb_outdir = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-stats_outdir") == 0) {
            this->stats_outdir = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-memb") == 0) {
            this->memb_file = argv[ii + 1];
            ii += 2;
        }
        else if (strcmp(argv[ii], "-mask") == 0) {
            this->mask_file = argv[ii + 1];
            ii += 2;
        }

        //*****************************************************************/
        // Parameters for exporting to CSV format used by clustering tools
        //*****************************************************************/

        else if (strcmp(argv[ii], "-exp_cl_csv_frags") == 0) {
            this->should_export_cl_csv_frag_lib = true;
            this->export_cl_csv_frag_lib_in_path = argv[ii + 1];
            this->export_cl_csv_frag_lib_out_path = argv[ii + 2];
            ii += 3;
        }


        //*****************************************************************/
        // Parameters for exporting directories containing simulated loops
        // to CSV format used by clustering tools
        //*****************************************************************/

        else if (strcmp(argv[ii], "-exp_cl_csv_pdbs") == 0) {
            this->should_export_cl_csv_pdbs = true;
            // parse csv list
            split(argv[ii + 1], ',', this->export_cl_csv_pdb_in_dirs);
            // remove leading and trailing whitespace and quotes
            for (auto i = this->export_cl_csv_pdb_in_dirs.begin();
                 i != this->export_cl_csv_pdb_in_dirs.end(); ++i) {
                trim(*i);
                trim(*i, "\"");
            }
            this->export_cl_csv_pdb_out_path = argv[ii + 2];
            this->export_cl_csv_pdb_out_ls_path = argv[ii + 3];
            this->export_cl_csv_pdb_max = atoi(argv[ii + 4]);
            this->export_cl_csv_pdb_bbone_only = (atoi(argv[ii + 5]) == 1);
            ii += 6;
        }

        //**********************/
        // Mutations
        //**********************/

        else if (strcmp(argv[ii], "-muts") == 0) {
            // read mutation file path
            std::string mutfpath = argv[ii + 1];
            // load mutation file
            this->muts.load(mutfpath);
            ii += 2;
        }

        //**********************/
        // Unrecognized tokens
        //**********************/

        else {
            string tmpStr;
            tmpStr = argv[ii];
            if (tmpStr.substr(tmpStr.length() - 4, 4) == ".pdb") {
                this->prot_file = tmpStr;
                ++ii;
            }
            else {
                cout << "Unrecognized command: " << argv[ii] << endl;
                Usage();
            }
        }
    }
}

/**
 * Print parameters
 */
void Params::print() const {
    cout << "Params:\n";
    cout << "\tprotein file: " << this->prot_file << " (" << this->prot_name << ")\n";
    cout << "\tpdb_outdir: " << this->pdb_outdir << endl;
    cout << "\tjob_prefix: " << this->job_prefix << endl;
    cout << "\trand_seed: " << this->rand_seed << endl;
    cout << "\tnumber of distance states: " << this->num_dist_states << endl;
    cout << "\tnumber of side chains states: " << this->num_sc_states << endl;
    cout << "\tnumber of conformations to attempt: " << this->num_conf << endl;
    cout << "\tnumber of conformations to keep: " << this->conf_keep << endl;
    cout << "\tnumber of conformations to write: " << this->output_pdb << endl;
    cout << "\tmaximum number of backbone clashes allowed per sample: " << this->max_num_bb_clashes << endl;
    cout << "\tmaximum number of side chain clashes allowed per sample: " << this->max_num_sc_clashes << endl;
    cout << "\tmaximum number of loop restarts per multi-loop sample: " << this->max_muloop_restart << endl;
    cout << "\tmaximum number of checked fragments during multi-loop restart: " << this->max_frag_checks << endl;
    cout << "\tuse backbone rotamer library: " << this->use_rot_lib << endl;
    cout << "\tstats directory: " << this->stats_outdir << endl;
    cout << "\tscore only: " << this->score_only << endl;
    if (this->should_export_cl_csv_frag_lib) {
        cout << "\tExport Cl CSV Lib In: " << this->export_cl_csv_frag_lib_in_path << endl;
        cout << "\tExport Cl CSV Lib Out: " << this->export_cl_csv_frag_lib_out_path << endl;
    }
    if (this->use_multiloop_lib) {
        cout << "\tMulti-loop lib dir: " << this->mulib_frag_dir << endl;
        cout << "\tMulti-loop ini: " << this->mulib_frag_ini_path << endl;
    }
}
