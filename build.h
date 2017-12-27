/**
 * Header meant to encapsulate global build flags
 */

#ifndef Build_h
#define Build_h

/**
 * Switch which enables side chains to be placed during grow_one calls
 * (Useful if implementing a physical potential with electrostatics)
 */
#ifndef DISGRO_BUILD_ENABLE_SC_GROW_ONE
    // Uncomment to enable side chain placement during grow_one() calls
    //#define DISGRO_BUILD_ENABLE_SC_GROW_ONE
#endif // DISGRO_BUILD_ENABLE_SC_GROW_ONE

/**
 * If enabled, will perform clash detection against all candidates and
 * filter candidates which are not clash free
 */
#ifndef DISGRO_BUILD_ENABLE_CLASH_FREE_BACKBONE
    // Uncomment to enable clash free backbone checking
    #define DISGRO_BUILD_ENABLE_CLASH_FREE_BACKBONE
#endif

/**
 * How often to report number of sample attempts vs closed conformations
 * @TODO - make this a command line parameter
 */
#define DISGRO_SMC_REPORT_INTERVAL 500

#ifndef DISGRO_BUILD_ENABLE_COLLISION_LOGGER
    // Uncomment to enable collision logging
    //#define DISGRO_BUILD_ENABLE_COLLISION_LOGGER
#endif // DISGRO_BUILD_ENABLE_COLLISION_LOGGER

#endif // Build_h
