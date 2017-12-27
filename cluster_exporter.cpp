/**
 * Utilities for exporting fragment libraries and full PDB folders to CSV
 * file usable by clustering tools
 */

#include "build.h"
#include "params.h"
#include "frag_lib_manager.h"
#include "pdb_reader.h"
#include "dirent_xplat.h"

#include <limits.h>
#include <fstream>

namespace {

    /**
     * Writes a collection of point arrays to a single CSV file
     * @return true if export successful, false o/w
     */
    bool write_pdbs_to_csv(const std::vector<SlimStruct> &confs, const std::string &outpath) {

        // Do nothing if nothing read
        const int num_confs = (int)confs.size();
        if (num_confs <= 0) {
            std::cout << "No conformations read. Nothing to export.\n";
            return false;
        }

        const int conf_len = static_cast<int>(confs.front().size());
        if (conf_len <= 0) {
            std::cout << "Conformations are of length 0. Nothing to export.\n";
            return false;
        }

        std::ofstream fout(outpath);
        if (!fout.good()) {
            std::cout << "Trouble opening export path: " << outpath << std::endl;
            std::cout << "No export made.\n";
            return false;
        }

        // Output format is:
        // row_1: conf_1.atom_1.x, conf_1.atom_1.y, conf_1.atom_1.z, .... conf_n.atom_1.x, conf_n.atom_1.y, conf_n.atom_1.z
        // ...
        // row_k: conf_1.atom_k.x, conf_k.atom_1.y, conf_k.atom_1.z, .... conf_n.atom_k.x, conf_n.atom_k.y, conf_n.atom_k.z

        std::cout << "Writing reformatted PDB set to: " << outpath << std::endl;
        for (int i_atom = 0; i_atom < conf_len; ++i_atom) {
            // Write first conformation
            fout << confs.front()[i_atom].x << "," << confs.front()[i_atom].y << "," << confs.front()[i_atom].z;
            // Write remaining fragments
            for (int i_struct = 1; i_struct < num_confs; ++i_struct) {
                const SlimStruct &conf = confs[i_struct];
                if (conf.size() <= i_atom) {
                    std::cout << "Trouble exporting - conformations not of same size.\n";
                    return false;
                }
                const Point &atom = conf[i_atom];
                fout << "," << atom.x << "," << atom.y << "," << atom.z;
            }
            fout << std::endl;
        }

        fout.close();
        return true;
    }

    /**
     * Based on discussion:
     * http://stackoverflow.com/questions/874134/find-if-string-ends-with-another-string-in-c
     * @return 1 iff str ends with suffix, 0 o/w 
     */
    inline bool str_ends_with(const char *str, const size_t str_len, const char *suffix, const size_t suffix_len) {
        assert(str != NULL);
        assert(suffix != NULL);
        assert(str_len == strlen(str));
        assert(suffix_len == strlen(suffix));
        return suffix_len <= str_len ?
            (0 == strncmp(str + str_len - suffix_len, suffix, suffix_len)) : false;
    }

    /**
     * Not the most robust utility. Expects a PDB name to have .pdb extension
     * and also at be longer than its extension string (i.e. - more than 4 chars)
     * @return true if we think file name is a PDB file, false o/w
     */
    inline bool is_pdb_fname(const char *fname) {
        const char *suffix = ".pdb";
        const size_t suffix_len = 4;
        const size_t fname_len = strlen(fname);
        return (fname_len > suffix_len) && str_ends_with(fname, fname_len, suffix, suffix_len);
    }

    /**
     * @return true if line is a PDB atom record for a backbone atom, false o/w
     */
    inline bool is_pdb_bbone_atom_rec(const std::string &line) {
        std::string atom_type;
        Pdb::get_atom_rec_token(atom_type, line, PdbAtomName);
        return (atom_type == "N")
            || (atom_type == "CA")
            || (atom_type == "C")
            || (atom_type == "O")
            || (atom_type == "CB");
    }

    /**
     * Utility parses a PDB file and extracts loop regions
     * @param conf - The empty point array to load PDB loop regions into
     * @param pdb_path - path to PDB file
     * @param starts - array of start loop residue sequence numbers
     * @param ends - array of end loop residue sequence numbers
     * @param bbone_only - if true, only backbone atoms are kept, else side chains are included
     * @return true if PDB successfully loaded and parsed, false o/w
     */
    bool load_conf(SlimStruct &conf,
                   const std::string &pdb_path,
                   const std::vector<int> &starts,
                   const std::vector<int> ends,
                   const bool bbone_only) {
        assert(conf.empty());
        assert(!pdb_path.empty());
        assert(!starts.empty());
        assert(!ends.empty());

        // Early out if PDB path is invalid
        if (!FileExists(pdb_path)) {
            std::cout << "Warning: PDB file not found: " << pdb_path << std::endl;
            return false;
        }

        // Parse PDB file line by line
        Point atom;
        std::ifstream fpdb(pdb_path);
        for (std::string line; std::getline(fpdb, line);) {

            // Skip empty lines
            if (line.empty()) {
                continue;
            }

            // Skip non-atom records
            if (!Pdb::is_atom_rec(line)) {
                continue;
            }

            // Skip if non-loop region
            std::string res_seq_no_str;
            Pdb::get_atom_rec_token(
                res_seq_no_str,
                line,
                PdbAtomResSeqNo);
            const int res_seq_no = atoi(res_seq_no_str.c_str());
            if (!is_in_simulated_region(res_seq_no, starts, ends)) {
                continue;
            }
            
            // Skip if filtering side chains and record is side chain atom
            if (bbone_only && !is_pdb_bbone_atom_rec(line)) {
                continue;
            }

            // Parse atom coordinates
            Pdb::parse_atom_coords(atom, line);

            // Add atom to conformation
            conf.push_back(atom);

        } // end iteration over lines of PDB file

        // Close file handle
        fpdb.close();
        return true;
    }
} // End of helper namespace

/**
 * Exports a fragment library to a CSV format expected by clustering tools
 */
void do_export_cl_csv_frag_lib(const Params &params) {
    std::vector<SlimStruct> Looplib;
    FragLibManager::load(Looplib, params.export_cl_csv_frag_lib_in_path, true /*ignore_H_atoms*/, true /*backbone_only*/);

    if (write_pdbs_to_csv(Looplib, params.export_cl_csv_frag_lib_out_path)) {
        std::cout << "Fragment library export finished.\n";
    }
}

/**
 * Exports a set of directories each containing PDB files to a single
 * CSV file for use with other tool such as clustering.
 * 
 * Note: Any file with .pdb extension in parameter folder is processed.
 * No attempt is made to check if file is actually a folder and not a PDB
 *
 * Note: INPUT DIRECTORY PATHS MUST NOT END WITH A PATH SEPARATOR CHAR
 *
 * Note: starts, ends loop region vectors MUST NOT BE REMAPPED. They must
 * have same residue sequence number as that in PDB
 */
void do_export_cl_csv_pdb(const Params &params) {

    // Store loop regions as point arrays
    std::vector<SlimStruct> confs;

    // File path to each PDB used for clustering
    std::vector<std::string> confs_ls;

    // Determine maximum number of PDBs to load
    const int MAX_CONFS = (params.export_cl_csv_pdb_max > 0)
        ? params.export_cl_csv_pdb_max : INT_MAX;

    // Load loops
    int num_confs = 0;
    for (auto idir = params.export_cl_csv_pdb_in_dirs.cbegin();
         idir != params.export_cl_csv_pdb_in_dirs.cend(); ++idir) {

        const std::string &in_dir_path = *idir;
        std::cout << "Loading PDBs from: " << in_dir_path << std::endl;

        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir(in_dir_path.c_str())) != NULL) {
            // Process all PDB files in this directory
            while ((ent = readdir(dir)) != NULL) {
                
                // Skip files without .pdb extension
                if (!is_pdb_fname(ent->d_name)) {
                    continue;
                }

                // Assume input directory does end in path separator.
                // Determine path to this PDB file
                const std::string pdb_path = in_dir_path + "/" + std::string(ent->d_name);

                // Load 
                confs.push_back(SlimStruct());
                if (load_conf(
                          confs.back()
                        , pdb_path
                        , params.starts
                        , params.ends
                        , params.export_cl_csv_pdb_bbone_only)) {
                    ++num_confs;
                    confs_ls.push_back(pdb_path);
                }
                else {
                    std::cout << "Error: unable to load " << pdb_path << "std::endl";
                    closedir(dir);
                    return;
                }

                if (num_confs == MAX_CONFS) {
                    std::cout << "Warning: Maximum number of conformations reached.\n";
                    break;
                }

                // Heartbeat
                if ((num_confs % DISGRO_SMC_REPORT_INTERVAL) == 0) {
                    std::cout << "Processed file #" << num_confs << ": " << pdb_path << std::endl;
                }
            }
            closedir(dir);
        }
        else {
            std::cout << "Error: could not open PDB directory: " << in_dir_path;
            return;
        }

    } // End iteration over PDB  input directories

    assert(confs.size() == confs_ls.size());

    if (!write_pdbs_to_csv(confs, params.export_cl_csv_pdb_out_path)) {
        // Export failed"
        return;
    }

    // write ls file
    std::ofstream fls(params.export_cl_csv_pdb_out_ls_path);
    if (!fls.good()) {
        std::cout << "Error: failed to write cluster file list at: "
                  << params.export_cl_csv_pdb_out_ls_path << std::endl;
        return;
    }
    for (auto ils = confs_ls.cbegin(); ils != confs_ls.cend(); ++ils) {
        fls << *ils << std::endl;
    }
    fls.close();

    std::cout << "Full PDB to cluster CSV export finished successfully.\n";
}
