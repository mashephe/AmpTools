/*
* Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
*
* Please refer to the NVIDIA end user license agreement (EULA) associated
* with this source code for terms and conditions that govern your use of
* this software. Any use, reproduction, disclosure, or distribution of
* this software and related documentation outside the terms of the EULA
* is strictly prohibited.
*
* CUDA Data-Parallel Primitives Library (CUDPP) is the proprietary
* property of The Regents of the University of California ("The
* Regents") and NVIDIA Corporation ("NVIDIA").
* 
* Copyright (c) 2007 The Regents of the University of California, Davis
* campus and NVIDIA Corporation. All Rights Reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
* 
* - Redistributions of source code must retain the above copyright
*   notice, this list of conditions and the following disclaimer. 
*
* - Redistributions in binary form must reproduce the above copyright
*   notice, this list of conditions and the following disclaimer in the
*   documentation and/or other materials provided with the distribution. 
* 
* - Neither the name of The Regents, NVIDIA, nor the names of its
*   contributors may be used to endorse or promote products derived from
*   this software without specific prior written permission.
*
* The end-user understands that the program was developed for research
* purposes and is advised not to rely exclusively on the program for any
* reason.
*
* THE SOFTWARE PROVIDED IS ON AN "AS IS" BASIS, AND THE REGENTS, NVIDIA
* AND CONTRIBUTORS HAVE NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
* UPDATES, ENHANCEMENTS, OR MODIFICATIONS. THE REGENTS, NVIDIA AND
* CONTRIBUTORS SPECIFICALLY DISCLAIM ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE REGENTS, NVIDIA OR CONTRIBUTORS BE LIABLE TO ANY
* PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, EXEMPLARY OR
* CONSEQUENTIAL DAMAGES, INCLUDING BUT NOT LIMITED TO PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES, LOSS OF USE, DATA OR PROFITS, OR
* BUSINESS INTERRUPTION, HOWEVER CAUSED AND UNDER ANY THEORY OF
* LIABILITY WHETHER IN CONTRACT, STRICT LIABILITY OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
* 
* If you do not agree to these terms, do not download or use the
* software.  This license may be modified only in a writing signed by
* authorized signatory of all parties.  For The Regents contact
* copyright@ucdavis.edu.
*
* Relating to funding received by the Regents-
* 
* Acknowledgment: This material is based upon work supported by
* the Department of Energy under Award Numbers DE-FG02-04ER25609 and
* DE-FC02-06ER25777.
*
* Disclaimer: This report was prepared as an account of work sponsored
* by an agency of the United States Government. Neither the United
* States Government nor any agency thereof, nor any of their employees,
* makes any warranty, express or implied, or assumes any legal liability
* or responsibility for the accuracy, completeness, or usefulness of any
* information, apparatus, product, or process disclosed, or represents
* that its use would not infringe privately owned rights. Reference
* herein to any specific commercial product, process, or service by
* trade name, trademark, manufacturer, or otherwise does not necessarily
* constitute or imply its endorsement, recommendation, or favoring by
* the United States Government or any agency thereof. The views and
* opinions of authors expressed herein do not necessarily state or
* reflect those of the United States Government or any agency hereof.
*/

/*
    Parallel reduction kernels
*/

#ifndef _REDUCE_KERNEL_H_
#define _REDUCE_KERNEL_H_

#include <stdio.h>

#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T*()
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

    __device__ inline operator const T*() const
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }
};

// specialize for double to avoid unaligned memory 
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double*()
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }

    __device__ inline operator const double*() const
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }
};

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
*/

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)

    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory. 
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.  
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <class T, unsigned int blockSize, bool nIsPow2>
__global__ void
reduce6(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    
    T mySum = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {         
        mySum += g_idata[i];
        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n) 
            mySum += g_idata[i+blockSize];  
        i += gridSize;
    } 

    // each thread puts its local sum into shared memory 
    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] = mySum = mySum + sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >= 512 ) { if (tid < 256) { sdata[tid] = mySum = mySum + sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256 ) { if (tid < 128) { sdata[tid] = mySum = mySum + sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128 ) { if (tid <  64) { sdata[tid] = mySum = mySum + sdata[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        // now that we are using warp-synchronous programming (below)
        // we need to declare our shared memory volatile so that the compiler
        // doesn't reorder stores to it and induce incorrect behavior.
        volatile T* smem = sdata;
        if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; EMUSYNC; }
        if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; EMUSYNC; }
        if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; EMUSYNC; }
        if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; EMUSYNC; }
        if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; EMUSYNC; }
        if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; EMUSYNC; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) 
        g_odata[blockIdx.x] = sdata[0];
}


////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void 
reduce(int size, int threads, int blocks, 
       T *d_idata, T *d_odata)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

    switch (threads)
    {
      case 1024:
          reduce6<T, 1024, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 512:
          reduce6<T, 512, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 256:
           reduce6<T, 256, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 128:
           reduce6<T, 128, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 64:
          reduce6<T,  64, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 32:
          reduce6<T,  32, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case 16:
          reduce6<T,  16, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case  8:
          reduce6<T,   8, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case  4:
          reduce6<T,   4, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case  2:
          reduce6<T,   2, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
      case  1:
          reduce6<T,   1, true><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;

    }
}

// Instantiate the reduction function for 3 types
template void 
reduce<int>(int size, int threads, int blocks, 
            int *d_idata, int *d_odata);

template void 
reduce<float>(int size, int threads, int blocks, 
              float *d_idata, float *d_odata);

template void 
reduce<double>(int size, int threads, int blocks, 
               double *d_idata, double *d_odata);

#endif // #ifndef _REDUCE_KERNEL_H_
