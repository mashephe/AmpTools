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

#ifndef __GPU_MANAGER_H__
#define __GPU_MANAGER_H__

#include <vector>
#include <complex>
#include <map>
#include <cassert>

#include "GPUManager/GPUCustomTypes.h"

using namespace std;
using namespace __gnu_cxx;

class Kinematics;
class AmplitudeManager;
class AmpVecs;
class Amplitude;

class GPUManager
{
  
public:
  
  GPUManager();
  GPUManager( const AmpVecs& a );
  ~GPUManager();
  
  void clearAll();
  void clearAmpCalc();
  void clearLikeCalc();
  
  void init( const AmpVecs& a, bool ampCalcOnly = false );
  
  // Interface Utils
  // First Amplitude calculation interface
  void copyDataToGPU( const AmpVecs& a );
  
  void calcAmplitudeAll( const Amplitude* amp, GDouble* pcResAmp, 
                         const vector< vector< int > >* pvPermutations );
  
  // Now the intensity calculator
  void copyAmpsToGPU( const AmpVecs& a );
  double calcSumLogIntensity( const vector< complex< double > >& prodCoef,
                              const vector< vector< bool > >& cohMtx );
  
  // General utils:
  static int calcNEventsGPU( int iNEvents ){
    
    //Should be a power of 2 for reduction to work, also multiple of GPU_BLOCK_SIZE_SQ    
    int iPow = 0;
    while( ( 1 << iPow ) < iNEvents ) iPow++;
    return 1 << iPow; 
  }
  
private:
  
  static bool m_cudaDisplay;
  
  bool m_ampCalcOnly;
  
  // array dimensions
  unsigned int m_iNParticles;
  unsigned int m_iNEvents;
  unsigned int m_iNTrueEvents;
  unsigned int m_iNAmps;
  unsigned int m_iNAmpsH;
  
  // array sizes
  unsigned int m_iAmpArrSize;
  unsigned int m_iEventArrSize;
  unsigned int m_iTrueEventArrSize;
  unsigned int m_iVArrSize;
  
  //Host Arrays
  GDouble* m_pcCalcAmp;
  
  GDouble* m_pfAmpRe;
  GDouble* m_pfAmpIm;
  GDouble* m_pfVRe;
  GDouble* m_pfVIm;
  GDouble* m_pfRes;
  
  //Device Arrays 
  GDouble* m_pfDevData;
  GDouble* m_pfDevWeights;
  GDouble* m_pcDevCalcAmp;
  int*     m_piDevPerm;
  
  GDouble* m_pfDevAmpRe;
  GDouble* m_pfDevAmpIm;
  GDouble* m_pfDevVRe;
  GDouble* m_pfDevVIm;
  
  // intensity sums maintained at double precision
  GDouble* m_pfDevRes;
  GDouble* m_pfDevREDUCE;
  
  // CUDA Thread and Grid sizes
  unsigned int m_iDimGridX;
  unsigned int m_iDimGridY;
  unsigned int m_iDimThreadX;
  unsigned int m_iDimThreadY;
  
  unsigned int m_iNBlocks; 
  unsigned int m_iNThreads;
  
  // Internal Utils
  void calcCUDADims();
  
};

#endif //__GPU_MANAGER_H__
