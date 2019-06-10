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

#include <stdlib.h>
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"

#include "cuda_runtime.h"

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr,"GPU ERROR: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

using namespace std;
// using namespace __gnu_cxx;

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
  
  void init( const AmpVecs& a, bool use4Vectors = true );
  
  // Interface Utils
  // First Amplitude calculation interface
  void copyDataToGPU( const AmpVecs& a, bool use4Vectors = true  );
  void copyUserVarsToGPU( const AmpVecs& a );
  
  void calcAmplitudeAll( const Amplitude* amp, unsigned long long offset,
                         const vector< vector< int > >* pvPermutations,
			 unsigned long long userVarsOffset );
  
  void assembleTerms( int iAmpInd, int nFact, int nPerm );
  
  void copyAmpsFromGPU( AmpVecs& a );

  // Now the intensity calculator
  double calcSumLogIntensity( const vector< complex< double > >& prodCoef,
                              const vector< vector< bool > >& cohMtx );

  void calcIntegral( GDouble* result, int iAmp, int jAmp, int iNGenEvents );
  
  // General utils:
  static int calcNEventsGPU( int iNEvents ){
    
    //Should be a power of 2 for reduction to work, also multiple of GPU_BLOCK_SIZE_SQ    
    int iPow = 0;
    while( ( 1 << iPow ) < iNEvents ) iPow++;
    return 1 << iPow; 
  }
  
private:
  
  static bool m_cudaDisplay;
  
  // array dimensions
  unsigned int m_iNParticles;
  unsigned long long m_iNEvents;
  unsigned long long m_iNTrueEvents;
  unsigned int m_iNAmps;
  unsigned int m_iNUserVars;
  
  // array sizes
  unsigned long long m_iEventArrSize;
  unsigned long long m_iTrueEventArrSize;
  unsigned long long m_iAmpArrSize;
  unsigned int m_iVArrSize;
  
  //Host Arrays
  GDouble* m_pcCalcAmp;
  
  GDouble* m_pfVVStar;
  GDouble* m_pfRes;
  
  //Device Arrays 
  GDouble* m_pfDevData;
  GDouble* m_pfDevUserVars;
  GDouble* m_pfDevWeights;
  GDouble* m_pcDevCalcAmp;
  int*     m_piDevPerm;
  
  GDouble* m_pfDevAmps;
  GDouble* m_pfDevVVStar;
  
  GDouble* m_pfDevResRe;
  GDouble* m_pfDevResIm;
  GDouble* m_pfDevREDUCE;
  
  // CUDA Thread and Grid sizes
  unsigned int m_iDimGridX;
  unsigned int m_iDimGridY;
  unsigned int m_iDimThreadX;
  unsigned int m_iDimThreadY;
  
  unsigned int m_iNBlocks; 
  unsigned int m_iNThreads;
  
  // Internal Utils
  
  unsigned int m_devProp_major;

  void calcCUDADims();
  
};

#endif //__GPU_MANAGER_H__
