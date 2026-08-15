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
#include "IUAmpTools/report.h"

#include "cuda_runtime.h"

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

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
  ~GPUManager();
  
  void clearAll();
  void clearData();
  void clearTerms();
  
  void init( const AmpVecs& a, bool use4Vectors = true );
  
  void initData( const AmpVecs& a, bool use4Vectors = true );
  void useDataFrom( const AmpVecs& a );
  void initTerms( const AmpVecs& a, size_t chunkSize = 0 );
  
  void copyDataToGPU( const AmpVecs& a, bool use4Vectors = true  );
  void copyUserVarsToGPU( const AmpVecs& a );
  
  void calcAmplitudeAll( const Amplitude* amp, size_t uAmpFactOffset,
                         const vector< vector< int > >* pvPermutations,
                         size_t userVarsOffset,
                         size_t startEvent = 0, size_t chunkSize = 0 );
  
  void assembleTerms( int iAmpInd, int nFact, int nPerm, size_t nEvents );
  
  void copyAmpsFromGPU( AmpVecs& a );

  // Now the intensity calculator
  double calcSumLogIntensity( const vector< complex< double > >& prodCoef,
                              const vector< vector< bool > >& cohMtx );

  void calcIntegrals( double* result, int nElements, 
                      const vector<int>& iIndex, const vector<int>& jIndex, 
                      size_t startEvent, size_t nEvents );
  
  // General utils:
  static size_t calcNEventsGPU( size_t iNEvents ){
    
    //Should be a power of 2 for reduction to work, also multiple of GPU_BLOCK_SIZE_SQ    
    int iPow = 0;
    while( ( 1L << iPow ) < iNEvents ) iPow++;
    return(  ( 1L << iPow ) < GPU_BLOCK_SIZE_SQ ? GPU_BLOCK_SIZE_SQ : ( 1L << iPow ) );
  }
  
  bool m_ownsData;

private:
  
  static bool m_cudaDisplay;
  
  // array dimensions
  size_t m_iNParticles;
  size_t m_iNEvents;
  size_t m_iNTrueEvents;
  size_t m_iNAmps;
  size_t m_iNUserVars;
  
  // array sizes (data)
  size_t m_iGDoubleDataArrSize;
  
  // array sizes (terms)
  size_t m_iDoubleIntenArrSize;
  size_t m_iDoubleReduceArrSize;
  size_t m_iGDoubleFactArrSize;
  size_t m_iAmpArrSize;
  size_t m_iVArrSize;
  size_t m_iNICalcSize;

  //Host Arrays
  GDouble* m_pfVVStar;
  double* m_pdRes;
  
  //Device Arrays 
  GDouble* m_pfDevData;
  GDouble* m_pfDevUserVars;
  GDouble* m_pfDevWeights;
  GDouble* m_pcDevAmpFact;
  int*     m_piDevPerm;
  
  GDouble* m_pfDevAmps;
  GDouble* m_pfDevVVStar;
  double* m_pdDevNICalc;
  
  double* m_pdDevRes;
  double* m_pdDevREDUCE;
  
  // CUDA Thread and Grid sizes
  size_t m_iDimGridX;
  size_t m_iDimGridY;
  size_t m_iDimGridXAmpFact;
  size_t m_iDimGridYAmpFact;
  size_t m_iDimThreadX;
  size_t m_iDimThreadY;
  
  size_t m_iNBlocks;
  size_t m_iNThreads;
  
  // Internal Utils
  
  unsigned int m_devProp_major;
  unsigned int m_devProp_minor;
  size_t m_maxShared_bytes;

  void calcCUDADims();
  void calcCUDADimsAmpFact( size_t chunkSize = 0 );

  static const char* kModule;
};

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) {

    report( ERROR, "GPUManager" ) << "CUDA ERROR: " << cudaGetErrorString( code ) << " in " 
                                  << file << " at line " << line << endl;

    if (abort) exit(code);
  }
}

#endif //__GPU_MANAGER_H__
