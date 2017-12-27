#ifndef EnergyStatsWriter_h
#define EnergyStatsWriter_h

#include "structure.h"
#include <vector>
#include <string>

// Base class for writing energy stats to disk
class EnergyStatsWriterBase {
public:

    // Utility for writing energy stats to file
    void write(const std::vector<Structure*> &Topconflist, const struct Params &params);

protected:

    // Initializes pdb file name
    virtual void init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                               const struct Params &params) = 0;

    // Called within each loop iteration when writing stats for each top conformation
    virtual void update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                 const struct Params &params, const int i) = 0;
};

// Specialized class for writing fragment library stats
class EnergyStatsWriterFragLib : public EnergyStatsWriterBase {
private:
    // Initializes pdb file name
    virtual void init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                               const struct Params &params);

    // Called within each loop iteration when writing stats for each top conformation
    virtual void update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                 const struct Params &params, const int i);
};

// Specialized class for writing multi loop top conformations stats
class EnergyStatsWriterMultiLoops : public EnergyStatsWriterBase {
private:
    // Initializes pdb file name
    virtual void init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                               const struct Params &params);

    // Called within each loop iteration when writing stats for each top conformation
    virtual void update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                 const struct Params &params, const int i);
};

// Specialized class for writing parameter (native) pdb energy stats
class EnergyStatsWriterScoreOnly : public EnergyStatsWriterBase {
private:
    // Initializes pdb file name
    virtual void init_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                               const struct Params &params);

    // Called within each loop iteration when writing stats for each top conformation
    virtual void update_pdb_name(std::string &out_pdbname, const std::vector<Structure*> &Topconflist,
                                 const struct Params &params, const int i);
};

#endif // EnergyStatsWriter_h
