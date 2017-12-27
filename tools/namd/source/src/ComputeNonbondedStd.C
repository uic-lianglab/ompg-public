/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

// DMK - CHECK/DEBUG - Atom Separation (water vs. non-water)
#include "common.h"
#include "NamdTypes.h"
// BEGIN INTERACTION RECORDER MOD
#include "InteractionRecorder.h"    
// END INTERACTION RECORDER MOD
#if NAMD_SeparateWaters != 0
  #define DEFINE_CHECK_WATER_SEPARATION
#endif


#include "ComputeNonbondedInl.h"

#define NBTYPE NBPAIR
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBaseStd.h"
#define CALCENERGY
#include "ComputeNonbondedBaseStd.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE


#define INTFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBaseStd.h"
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#undef MERGEELECT
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBaseStd.h"
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#undef MERGEELECT
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY
#undef INTFLAG

// moved to ComputeNonbondedPprof.C
#if 0

#define PPROFFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBaseStd.h"
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#undef MERGEELECT
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBaseStd.h"
#define FULLELECT
#include "ComputeNonbondedBaseStd.h"
#define MERGEELECT
#include "ComputeNonbondedBaseStd.h"
#undef MERGEELECT
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY
#undef PPROFFLAG

#endif

