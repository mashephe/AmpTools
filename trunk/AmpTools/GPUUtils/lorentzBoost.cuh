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

#ifndef LORENTZBOOST
#define LORENTZBOOST

#include "GPUManager/GPUCustomTypes.h"

// a routine to boost p4Vec 

static __device__ void
boost( GDouble* p4Vec, GDouble* beta ){

  GDouble bx = beta[0];
  GDouble by = beta[1];
  GDouble bz = beta[2];

  GDouble bp2 = bx*bx + by*by + bz*bz;

  GDouble pt = p4Vec[0];
  GDouble px = p4Vec[1];
  GDouble py = p4Vec[2];
  GDouble pz = p4Vec[3];

  GDouble gamma = 1.0 / G_SQRT(1.0 - bp2);
  GDouble bgamma = gamma * gamma / (1.0 + gamma);
  
  // the boost transform
  GDouble mtt, mtx, mty, mtz;
  GDouble mxt, mxx, mxy, mxz;
  GDouble myt, myx, myy, myz;
  GDouble mzt, mzx, mzy, mzz;
  
  mtt = gamma;
  mxx = 1.0 + bgamma * bx * bx;
  myy = 1.0 + bgamma * by * by;
  mzz = 1.0 + bgamma * bz * bz;
  mxy = myx = bgamma * bx * by;
  mxz = mzx = bgamma * bx * bz;
  myz = mzy = bgamma * by * bz;
  mxt = mtx = gamma * bx;
  myt = mty = gamma * by;
  mzt = mtz = gamma * bz;
  
  p4Vec[0] = pt*mtt + px*mtx + py*mty + pz*mtz;
  p4Vec[1] = pt*mxt + px*mxx + py*mxy + pz*mxz;
  p4Vec[2] = pt*myt + px*myx + py*myy + pz*myz;
  p4Vec[3] = pt*mzt + px*mzx + py*mzy + pz*mzz;
}

// assuming that p4Vec1 and p4Vec2 are measured
// in the same frame this routine boosts the 
// Lorentz vector p4Vec1 to the rest frame of p4Vec2

static __device__ void
boostToRest( GDouble* p4Vec1, GDouble* p4Vec2 ){

  // the boost vector
  GDouble boostVec[] =  { -p4Vec2[1] / p4Vec2[0],
                          -p4Vec2[2] / p4Vec2[0],
                          -p4Vec2[3] / p4Vec2[0] };
  
  boost( p4Vec1, boostVec );
}

#endif // LORENTZBOOST

