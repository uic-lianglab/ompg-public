//
// Created by apr on 10/22/15.
//

#include "energy_stats.h"
#include "energy_stats_writer.h"
#include "params.h"
#include "util.h"

#include <fstream>

// Utility for getting an output pdb file name
extern void get_out_pdb_fname(string &out_fname, const Params &params, const int rank);

// Utility for getting pdb lib filename
extern void get_out_pdb_lib_fname(string &out_fname, const Params &params, const Structure &Conf);

//////////////////////////////////////////////
// Base Writer
//////////////////////////////////////////////

// Utility for writing energy stats to file
void EnergyStatsWriterBase::write(const vector<Structure*> &Topconflist, const Params &params) {
    // Determine stats file name
    const string stats_fname = params.stats_outdir + "/" + params.job_prefix + "_stats.csv";
    // Open file handle
    ofstream stats_file;
    const bool b_file_already_exists = FileExists(stats_fname);
    stats_file.open(stats_fname, std::ios::app);
    // Write header row
    if (!b_file_already_exists) {
        stats_file <<
        "pdb_name, e_loodis, energy" <<
        endl;
    }
    string out_pdbname;
    this->init_pdb_name(out_pdbname, Topconflist, params);
    for (int i = 0; i < Topconflist.size(); ++i) {
        const Structure &topconf = *(Topconflist[i]);
        this->update_pdb_name(out_pdbname, Topconflist, params, i);
        stats_file << out_pdbname << ", "
                   << topconf._enStats.loodis_term << ", "
                   << topconf._energy << endl;
    }
    stats_file.close();
}

//////////////////////////////////////////////
// Fragment Library Writer
//////////////////////////////////////////////

// Initializes pdb file name
void EnergyStatsWriterFragLib::init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                             const struct Params &params) {
    assert(!Topconflist.empty());
    assert(NULL != Topconflist[0]);
    get_out_pdb_lib_fname(out_pdbname, params, *(Topconflist[0]));
}

// Called within each loop iteration when writing stats for each top conformation
void EnergyStatsWriterFragLib::update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                               const struct Params &params, const int i) {
    /* do nothing */
}

//////////////////////////////////////////////
// Multi-loop Writer
//////////////////////////////////////////////

// Initializes pdb file name
void EnergyStatsWriterMultiLoops::init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                                const struct Params &params) {
    /* do nothing */
}

// Called within each loop iteration when writing stats for each top conformation
void EnergyStatsWriterMultiLoops::update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                                  const struct Params &params, const int i) {
    get_out_pdb_fname(out_pdbname, params, i);
}

//////////////////////////////////////////////
// Score Only Writer
//////////////////////////////////////////////

// Initializes pdb file name
void EnergyStatsWriterScoreOnly::init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                               const struct Params &params) {
    out_pdbname = params.prot_file;
}

// Called within each loop iteration when writing stats for each top conformation
void EnergyStatsWriterScoreOnly::update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                                 const struct Params &params, const int i) {
    /* do nothing */
}