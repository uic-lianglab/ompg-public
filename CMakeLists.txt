cmake_minimum_required(VERSION 2.8)
project(pretzel)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_COMPILER_IS_GNUCXX)
    # Add warning if non-void method is missing return statement
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
endif()

set(SOURCE_FILES
    ./source/atom.cpp
    ./source/cal_energy.cpp
    ./source/cluster_exporter.cpp
    ./source/collision_frontend.cpp
    ./source/collision_grid_prot.cpp
    ./source/energy_stats.cpp
    ./source/frag_lib_manager.cpp
    ./source/grow_one_util.cpp
    ./source/main.cpp
    ./source/matrix.cpp
    ./source/membrane.cpp
    ./source/mutation_manager.cpp
    ./source/params.cpp
    ./source/pdb_reader.cpp
    ./source/potential.cpp
    ./source/reprst.cpp
    ./source/residue.cpp
    ./source/rotamer.cpp
    ./source/sample_states.cpp
    ./source/smc.cpp
    ./source/structure.cpp
    ./source/util.cpp
    ./source/vdw_utils.cpp
    ./source/config/ConfigFile.cpp
)

set(HEADER_FILES
    ./source/atom.h
    ./source/build.h
    ./source/cal_energy.h
    ./source/collision_frontend.h
    ./source/collision_grid.h
    ./source/collision_grid_prot.h
    ./source/collision_logger.h
    ./source/dirent_win.h
    ./source/dirent_xplat.h
    ./source/energy_stats.h
    ./source/energy_stats_writer.h
    ./source/frag_lib_manager.h
    ./source/matrix.h
    ./source/membrane.h
    ./source/misc_structs.h
    ./source/mutation_manager.h
    ./source/params.h
    ./source/pdb_reader.h
    ./source/point.h
    ./source/potential_frontend_loodis.h
    ./source/potential.h
    ./source/reprst.h
    ./source/residue.h
    ./source/rotamer.h
    ./source/sample_states.h
    ./source/smc.h
    ./source/structure.h
    ./source/util.h
    ./source/vdw_utils.h
    ./source/vla.h
)

# Default executable name
set(TARGET_EXE_NAME "pretzel")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Note: for full release build, do not want to link any excluded files
    add_executable(${TARGET_EXE_NAME} ${SOURCE_FILES})
else()
    add_executable(${TARGET_EXE_NAME} ${SOURCE_FILES} ${HEADER_FILES})
endif()

