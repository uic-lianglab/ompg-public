#ifndef collision_logger_h
#define collision_logger_h

#include "build.h"

#ifdef DISGRO_BUILD_ENABLE_COLLISION_LOGGER

#include "residue.h"
#include "collision_grid_prot.h"

#include <string>
#include <vector>
#include <fstream>

// Output path for collision logs
#define DISGRO_COLLISION_FILE_PATH "./logs/collisions.csv"

class CollisionLogger
{
private:

    struct atom_detailed_info {
        int res_ix;
        int atom_ix;
        std::string res_name;
        std::string atom_name;
        double r_vdw;
        double x, y, z;
        void init(const Atom &atom,
                  const struct atom_info& a_info) {
            res_ix = a_info.res_ix;
            atom_ix = a_info.atom_ix;
            res_name = Residue::Name3[a_info.res_type];
            atom_name = atom._name;
            r_vdw = collision_grid_prot::elem_radius(atom);
            x = atom.x;
            y = atom.y;
            z = atom.z;
        }
    };

    struct collision_info {
        atom_detailed_info a;
        atom_detailed_info b;
        double dist;
        double min_dist;
        std::string tag;
    };

public:

    void record(const Atom &atom,
                const Atom &other_atom,
                const struct atom_info atom_infos[2],
                const char *tag) {
        collision_infos.push_back(struct collision_info());
        collision_info &c_info = collision_infos.back();
        c_info.a.init(atom, atom_infos[0]);
        c_info.b.init(other_atom, atom_infos[1]);
        c_info.dist = atom.dis(other_atom);
        c_info.min_dist = c_info.a.r_vdw + c_info.b.r_vdw;
        c_info.tag = tag;
    }

    void dump() {
        std::ofstream flog(DISGRO_COLLISION_FILE_PATH, std::ios::out);
        // Header row
        flog << "tag,rixa,aixa,rna,ana,rixb,aixb,rnb,anb,d,md,ofc,rva,rvb,xa,ya,za,xb,yb,zb\n";
        for (size_t i = 0; i < collision_infos.size(); ++i) {
            const collision_info &c_info = collision_infos[i];
            flog << c_info.tag << ","
                 << c_info.a.res_ix << ","
                 << c_info.a.atom_ix << ","
                 << c_info.a.res_name << ","
                 << c_info.a.atom_name << ","
                 << c_info.b.res_ix << ","
                 << c_info.b.atom_ix << ","
                 << c_info.b.res_name << ","
                 << c_info.b.atom_name << ","
                 << c_info.dist << ","
                 << c_info.min_dist << ","
                 << c_info.dist / c_info.min_dist << ","
                 << c_info.a.r_vdw << ","
                 << c_info.b.r_vdw << ","
                 << c_info.a.x << ","
                 << c_info.a.y << ","
                 << c_info.a.z << ","
                 << c_info.b.x << ","
                 << c_info.b.y << ","
                 << c_info.b.z << "\n";
        }
        flog.close();
    }

    // Implicitly dumps log on destruction
    ~CollisionLogger() {
        dump();
    }

private:

    std::vector<collision_info> collision_infos;
};

#define COLLISION_LOGGER_DECLARE(logger) mutable CollisionLogger logger;

#define COLLISION_LOGGER_RECORD(logger, atom_a, atom_b, atominfos, tag) \
    logger.record(atom_a, atom_b, atominfos, tag)

#else

#define COLLISION_LOGGER_DECLARE(logger)

#define COLLISION_LOGGER_RECORD(logger, atom_a, atom_b, atominfos, tag) ((void)0)

#endif // DISGRO_BUILD_ENABLE_COLLISION_LOGGER

#endif // collision_logger_h
