//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
// 
// Copyright Trustees of Indiana University 2010, all rights reserved
// 
// This software written by Matthew Shepherd, Ryan Mitchell, and 
//                  Hrayr Matevosyan at Indiana University, Bloomington
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// 3. Neither the name of the University nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// 
// Creation of derivative forms of this software for commercial
// utilization may be subject to restriction; written permission may be
// obtained from the Trustees of Indiana University.
// 
// INDIANA UNIVERSITY AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES, 
// EXPRESS OR IMPLIED.  By way of example, but not limitation, INDIANA 
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR 
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR 
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS, 
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be 
// held liable for any liability with respect to any claim by the user or 
// any other party arising from use of the program.
//******************************************************************************

#ifndef __GPU_CUSTOM_TYPES_H__ 
#define __GPU_CUSTOM_TYPES_H__

#include <iostream>

using namespace std;

#define GPU_MAX_AMPS 4095 // SHARED MEMORY ARRAY, MUST ACCOMODATE AiAj* for each event!
#define GPU_MAX_PART 10   // MAXIMUM NUMBER OF PARTICLES

#define GPU_BLOCK_SIZE_X 16
#define GPU_BLOCK_SIZE_Y 16
#define GPU_BLOCK_SIZE_SQ 256

// a few helpers to make operations on GPU a little more user friendly

// standard arguments to kernel launches
#define GPU_AMP_PROTO GDouble* pfDevData, GDouble* pfDevUserVars, \
                      WCUComplex* pcDevAmp, int* piDevPerm, int iNParticles, int iNEvents
#define GPU_AMP_ARGS  pfDevData, pfDevUserVars, pcDevAmp, piDevPerm, iNParticles, iNEvents

// how to get the event we are working with:
#define GPU_THIS_EVENT threadIdx.x + GPU_BLOCK_SIZE_X * threadIdx.y + \
                      ( blockIdx.x + blockIdx.y * gridDim.x ) * GPU_BLOCK_SIZE_SQ

// help for parsing arrays
#define GPU_E  0
#define GPU_PX 1
#define GPU_PY 2
#define GPU_PZ 3

#define GPU_KIN(id,val) pfDevData[ ( piDevPerm[id] * 4 + val ) * iNEvents + iEvent]
#define GPU_P4(id) { GPU_KIN(id,GPU_E)  , GPU_KIN(id,GPU_PX), \
                     GPU_KIN(id,GPU_PY) , GPU_KIN(id,GPU_PZ)}

#define GPU_UVARS(val) pfDevUserVars[val*iNEvents+iEvent]

#define COPY_P4(a1,a2) GDouble a2[4]; for( int zzz = 0; zzz < 4; ++zzz ) a2[zzz] = a1[zzz];


// Use this flag to turn on double precision
// (MPI portions of the code depend on this also)
// 
// NOTE: 64-bit floating point math is only
// supported on hardware with 1.3 or higher
// compute capability
//
// single precision should be used with caution;
// some intensity calculations may be sensitive to
// cancellation in terms which is difficult to 
// achieve numerically with single precision

#define USE_DOUBLE_PRECISION

#ifdef USE_DOUBLE_PRECISION  // double precision:

typedef double GDouble;

//MACROS FOR DIFFERENT PRECISION OPERATIONS
#define G_SIN(a)	sin(a)
#define G_ASIN(a) asin(a)
#define G_COS(a)  cos(a)
#define G_ACOS(a) acos(a)

#define G_ATAN2(a,b)	atan2(a,b)
#define G_ATAN(a) atan(a)
#define G_LOG(a)	log(a)

#define G_FABS(a)	fabs(a)
#define G_POW(a,b)	pow(a,b)
#define G_EXP(a) exp(a)
#define G_SQRT(a)	sqrt(a)

#define PI 3.141592653589793
#define D_EPS		1.e-10

#else // single precision:

typedef float GDouble;

//MACROS FOR DIFFERENT PRECISION OPERATIONS
#define G_SIN(a)	sinf(a)
#define G_ASIN(a) asinf(a)
#define G_COS(a)  cosf(a)
#define G_ACOS(a) acosf(a)

#define G_ATAN2(a,b)	atan2f(a,b)
#define G_ATAN(a) atanf(a)
#define G_LOG(a)	logf(a)

#define G_FABS(a)	fabsf(a)
#define G_POW(a,b)	powf(a,b)
#define G_EXP(a) expf(a)
#define G_SQRT(a)	sqrtf(a)

#define PI 3.141592653589793f
#define D_EPS		1.e-10f

#endif // USE_DOUBLE_PRECISION

//Shortcut macros
#define SQ(a) ((a)*(a)) 

//forward declaration
class dim3;

#endif /*__GPU_CUSTOM_TYPES_H__*/
