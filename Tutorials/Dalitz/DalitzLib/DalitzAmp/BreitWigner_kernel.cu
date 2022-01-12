
#include <stdio.h>

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"


__global__ void
BreitWigner_kernel( GPU_AMP_PROTO, GDouble mass0, GDouble width0, 
                       int daught1, int daught2 ){

  int iEvent = GPU_THIS_EVENT;

  /****
  *  Equivalently one could switch needsUserVarsOnly to false
  *  and uncomment the following lines to get the same
  *  result.  The role of user-data is to optimize this function
  *  call in the instance it is repeated multiple times throughout
  *  a fit.

  GDouble p1[4] = GPU_P4(daught1-1);
  GDouble p2[4] = GPU_P4(daught2-1);

  GDouble m = G_SQRT(SQ(p1[0]+p2[0]) - SQ(p1[1]+p2[1])
                                     - SQ(p1[2]+p2[2])
                                     - SQ(p1[3]+p2[3]));

  WCUComplex bwBot = { SQ( m ) - SQ( mass0 ), mass0 * width0 };
  */


  // here we need to be careful to index the user-defined
  // data with the proper integer corresponding to the
  // enumeration in the C++ header file

  GDouble m2 = GPU_UVARS(0);
  WCUComplex bwBot = { m2 - SQ( mass0 ), mass0 * width0 };

  pcDevAmp[iEvent] = ( 1.0 / bwBot );

}


void
BreitWigner_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO, 
                     GDouble mass, GDouble width,
                     int daught1, int daught2 )
{  

  BreitWigner_kernel<<< dimGrid, dimBlock >>>
    ( GPU_AMP_ARGS, mass, width, daught1, daught2 );

}
