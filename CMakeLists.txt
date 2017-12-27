cmake_minimum_required(VERSION 2.8)
project(mDisgro)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_COMPILER_IS_GNUCXX)
    # Add warning if non-void method is missing return statement
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
endif()

set(SOURCE_FILES
    atom.cpp
    cal_energy.cpp
    cluster_exporter.cpp
    collision_frontend.cpp
    collision_grid_prot.cpp
    energy_stats.cpp
    frag_lib_manager.cpp
    grow_one_util.cpp
    main.cpp
    matrix.cpp
    membrane.cpp
    mutation_manager.cpp
    params.cpp
    pdb_reader.cpp
    potential.cpp
    reprst.cpp
    residue.cpp
    rotamer.cpp
    sample_states.cpp
    smc.cpp
    structure.cpp
    util.cpp
    vdw_utils.cpp
    ./config/ConfigFile.cpp
)

set(HEADER_FILES
    atom.h
    build.h
    cal_energy.h
    collision_frontend.h
    collision_grid.h
    collision_grid_prot.h
    collision_logger.h
    dirent_win.h
    dirent_xplat.h
    energy_stats.h
    energy_stats_writer.h
    frag_lib_manager.h
    matrix.h
    membrane.h
    misc_structs.h
    mutation_manager.h
    params.h
    pdb_reader.h
    point.h
    potential_frontend_loodis.h
    potential.h
    reprst.h
    residue.h
    rotamer.h
    sample_states.h
    smc.h
    structure.h
    util.h
    vdw_utils.h
    vla.h
)

# Default executable name
set(TARGET_EXE_NAME "mDisgro")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Note: for full release build, do not want to link any excluded files
    add_executable(${TARGET_EXE_NAME} ${SOURCE_FILES})
else()
    add_executable(${TARGET_EXE_NAME} ${SOURCE_FILES} ${HEADER_FILES})
endif()

