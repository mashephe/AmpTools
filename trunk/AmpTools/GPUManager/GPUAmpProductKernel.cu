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

__global__ void
amp_kernel( GDouble* pfDevAmps, GDouble* pfDevVVStar, GDouble* pfDevWeights, 
            int nAmps, int nEvents, GDouble* pfDevRes )
{
	int i = threadIdx.x + GPU_BLOCK_SIZE_X * threadIdx.y + 
            ( blockIdx.x + blockIdx.y * gridDim.x ) * GPU_BLOCK_SIZE_SQ;

	int iA, iB;
	GDouble fSumRe = 0.;

  // index to amplitude A_alpha for the ith event
  int aIndA = 2*i;

	for( iA = 0; iA < nAmps; ++iA ){

  // index to amplitude A_alpha for the ith event
  // (reduce computations by incrementing at end of loop)
  //	  int aIndA = i + 2*nEvents*iA;

    // index in the array of V_alpha * conj( V_beta )
    int vInd = iA*(iA+1)/2;
    
    // index to amplitude A_beta for the ith event
    int aIndB = 2*i;
	  for( iB = 0; iB <= iA; ++iB ){

      // index in the array of V_alpha * conj( V_beta )
      // (reduce computations by incrementing at end of loop)
      //	   int vInd = iA*(iA+1)/2+iB;

	    // index to amplitude A_beta for the ith event
      // (reduce computations by incrementing at end of loop)
      //     int aIndB = i + 2*nEvents*iB;

      // only compute the real part of the intensity since
      // the imaginary part should sum to zero
	    GDouble term = pfDevVVStar[2*vInd] *
               ( pfDevAmps[aIndA]   * pfDevAmps[aIndB] +
                 pfDevAmps[aIndA+1] * pfDevAmps[aIndB+1] );

      term -= pfDevVVStar[2*vInd+1] *
               ( pfDevAmps[aIndA+1] * pfDevAmps[aIndB] -
                 pfDevAmps[aIndA]   * pfDevAmps[aIndB+1] );

	    // we're only summing over the lower diagonal so we need
      // to double the contribution for off diagonal elements
      if( iA != iB ) term *= 2;

 	    fSumRe += term;

      ++vInd;
      aIndB += 2*nEvents;
	  }
    
    aIndA += 2*nEvents;
	}
	
	pfDevRes[i] = pfDevWeights[i] * G_LOG( fSumRe );
}

extern "C" void GPU_ExecAmpKernel( dim3 dimGrid, dim3 dimBlock, 
     GDouble* pfDevAmps, GDouble* pfDevVVStar, GDouble* pfDevWeights,
     int nAmps, int nEvents, GDouble* pfDevRes )
{
	amp_kernel<<< dimGrid, dimBlock >>>( pfDevAmps, pfDevVVStar, pfDevWeights, 
                                       nAmps, nEvents, pfDevRes );
}

