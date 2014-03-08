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

__constant__ int da_iNAmps;  
__constant__ int da_iNEvents; // number of events padded to the closest 2^n
__constant__ GDouble da_pfDevVRe[GPU_MAX_AMPS]; // fit parameters stored in
__constant__ GDouble da_pfDevVIm[GPU_MAX_AMPS]; // const memory space

extern "C" GDouble* da_pfDevVRe_addr() { return da_pfDevVRe;  }
extern "C" GDouble* da_pfDevVIm_addr() { return da_pfDevVIm;  }
extern "C" int*     da_iNAmps_addr()   { return &da_iNAmps;   }
extern "C" int*     da_iNEvents_addr() { return &da_iNEvents; }



__global__ void
amp_kernel( GDouble* pfDevAmpRe, GDouble* pfDevAmpIm, 
            GDouble* pfDevWeights, GDouble* pfDevRes )
{
	int i = threadIdx.x + GPU_BLOCK_SIZE_X * threadIdx.y + 
            ( blockIdx.x + blockIdx.y * gridDim.x ) * GPU_BLOCK_SIZE_SQ;

	int iA, iB;
	GDouble fSumRe = 0.;

  // index to amplitude A_alpha for the ith event
  int aIndA = i;

	for( iA = 0; iA < da_iNAmps; ++iA ){

  // index to amplitude A_alpha for the ith event
  // (reduce computations by incrementing at end of loop)
  //	  int aIndA = i + da_iNEvents*iA;

    // index in the array of V_alpha * conj( V_beta )
    int vInd = iA*(iA+1)/2;
    
    // index to amplitude A_beta for the ith event
    int aIndB = i;
	  for( iB = 0; iB <= iA; ++iB ){

      // index in the array of V_alpha * conj( V_beta )
      // (reduce computations by incrementing at end of loop)
      //	   int vInd = iA*(iA+1)/2+iB;

	    // index to amplitude A_beta for the ith event
      // (reduce computations by incrementing at end of loop)
      //     int aIndB = i + da_iNEvents*iB;

      // only compute the real part of the intensity since
      // the imaginary part should sum to zero
	    GDouble term = da_pfDevVRe[vInd] * 
               ( pfDevAmpRe[aIndA] * pfDevAmpRe[aIndB] +
                 pfDevAmpIm[aIndA] * pfDevAmpIm[aIndB] );

      term -= da_pfDevVIm[vInd] *
               ( pfDevAmpIm[aIndA] * pfDevAmpRe[aIndB] -
                 pfDevAmpRe[aIndA] * pfDevAmpIm[aIndB] );

	    // we're only summing over the lower diagonal so we need
      // to double the contribution for off diagonal elements
      if( iA != iB ) term *= 2;

 	    fSumRe += term;

      ++vInd;
      aIndB += da_iNEvents;
	  }
    
    aIndA += da_iNEvents;
	}
	
	pfDevRes[i] = pfDevWeights[i] * G_LOG( fSumRe );
}

extern "C" void GPU_ExecAmpKernel( dim3 dimGrid, dim3 dimBlock, 
     GDouble* pfDevAmpRe, GDouble* pfDevAmpIm, GDouble* pfDevWeights, 
     GDouble* pfDevRes )
{
	amp_kernel<<< dimGrid, dimBlock >>>( pfDevAmpRe, pfDevAmpIm, 
                                             pfDevWeights, pfDevRes );
}




__global__ void
int_element_kernel( int iA, int iB, GDouble* pfDevAmpRe, GDouble* pfDevAmpIm,
                    GDouble* pfDevWeights, GDouble* pfDevResRe,
                    GDouble* pfDevResIm )
{
	int i = threadIdx.x + GPU_BLOCK_SIZE_X * threadIdx.y + 
            ( blockIdx.x + blockIdx.y * gridDim.x ) * GPU_BLOCK_SIZE_SQ;

  int aInd = i + da_iNEvents*iA;
  int bInd = i + da_iNEvents*iB;

  pfDevResRe[i] = pfDevAmpRe[aInd] * pfDevAmpRe[bInd]  +
                  pfDevAmpIm[aInd] * pfDevAmpIm[bInd];
  
  pfDevResIm[i] = pfDevAmpIm[aInd] * pfDevAmpRe[bInd] - 
                  pfDevAmpRe[aInd] * pfDevAmpIm[bInd];
}

extern "C" void GPU_ExecIntElementKernel( dim3 dimGrid, dim3 dimBlock,
     int iA, int iB, GDouble* pfDevAmpRe, GDouble* pfDevAmpIm, 
     GDouble* pfDevWeights, GDouble* pfDevResRe, GDouble* pfDevResIm )
{
	int_element_kernel<<< dimGrid, dimBlock >>>( iA, iB, pfDevAmpRe, pfDevAmpIm,
                                               pfDevWeights, pfDevResRe,
                                               pfDevResIm );
}

