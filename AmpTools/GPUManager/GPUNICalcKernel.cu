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

#include "GPUCustomTypes.h"
#include "stdio.h"

__global__ void
ni_calc_kernel( int nElements, GDouble* pfDevNICalc,
                GDouble* pfDevAmps, GDouble* pfDevWeights,
                int nEvents, int nTrueEvents )
{

  // used shared memory block for amplitude indices and results
  extern __shared__ int s[];
  
  unsigned int resultSize = 2*sizeof(GDouble)*nElements;
  unsigned int indexSize = 2*sizeof(int)*nElements;
  
  // the first thread of the block should copy the indices to
  // shared memory for other blocks to use and zero the result
  if( ( threadIdx.x + threadIdx.y ) == 0 ){

    memcpy( &s[2*nElements], &pfDevNICalc[2*nElements], indexSize );
    memset( s, 0, resultSize );
  }
  
  __syncthreads();
  
  // get addresses of arrays for the two indices in shared memory
  // and also the result
  int* iIndex = (int*)&pfDevNICalc[2*nElements];
  int* jIndex = &(iIndex[nElements]);
  GDouble* result = (GDouble*)s;

  int iEvt = threadIdx.x + GPU_BLOCK_SIZE_X * threadIdx.y +
            ( blockIdx.x + blockIdx.y * gridDim.x ) * GPU_BLOCK_SIZE_SQ;


  if( iEvt < nTrueEvents )
  for( int i = 0; i < nElements; ++i ){
  
    // these are the indices to the relevant amplitudes in the amplitude array
    int aInd = 2*iEvt + 2*nEvents*iIndex[i];
    int bInd = 2*iEvt + 2*nEvents*jIndex[i];

    GDouble thisRe, thisIm = 0;

    thisRe = pfDevWeights[iEvt] * (
                     pfDevAmps[aInd]   * pfDevAmps[bInd]  +
                     pfDevAmps[aInd+1] * pfDevAmps[bInd+1] );

    atomicAdd( &result[2*i], thisRe );

    if( aInd == bInd ) continue; // diagonal elements are real

    thisIm = pfDevWeights[iEvt] * (
                       pfDevAmps[aInd+1] * pfDevAmps[bInd] -
                       pfDevAmps[aInd]   * pfDevAmps[bInd+1] );

    atomicAdd( &result[2*i+1], thisIm );
  }
  
  __syncthreads();
  
  // now the first thread should accumulate global device memory
  // this can be made much faster by having each thread add one
  // element in shared memory to one element in device memory...
  // ... but let's try a simple approach to get it to work
  if( ( threadIdx.x + threadIdx.y ) == 0 ){

    for( int i = 0; i < 2*nElements; ++i ){

      atomicAdd( &pfDevNICalc[i], result[i] );
    }
  }
}


extern "C" void GPU_ExecNICalcKernel( dim3 dimGrid, dim3 dimBlock,
       	   			      unsigned int sharedSize,
                                      int nElements, GDouble* pfDevNICalc,
                                      GDouble* pfDevAmps, GDouble* pfDevWeights,
                                      int nEvents, int nTrueEvents )
{
  ni_calc_kernel<<< dimGrid, dimBlock, sharedSize >>>
     ( nElements, pfDevNICalc, pfDevAmps, pfDevWeights, nEvents, nTrueEvents );
}
