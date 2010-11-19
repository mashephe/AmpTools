/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.  This source code is a "commercial item" as 
 * that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer software" and "commercial computer software 
 * documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein.
 */

/*
    Parallel reduction kernels
*/


#ifndef __GPU_REDUCE_KERNEL_H__
#define __GPU_REDUCE_KERNEL_H__


#define SMVERSION sm13 //HHM ADDITION!!!

#include <stdio.h>
#include "GPUSharedMem.cuh"


#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

// Macros to append an SM version identifier to a function name
// This allows us to compile a file multiple times for different architecture
// versions
// The second macro is necessary to evaluate the value of the SMVERSION macro
// rather than appending "SMVERSION" itself


#define FUNCVERSION(x, y) x ## _ ## y
#define XFUNCVERSION(x, y) FUNCVERSION(x, y)
#define FUNC(NAME) XFUNCVERSION(NAME, SMVERSION) 

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
*/

/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved 
   inactivity means that no whole warps are active, which is also very 
   inefficient */
template <class T, unsigned int blockSize>
__global__ void
FUNC(reduce5)(T *g_idata, T *g_odata)
{
    SharedMemory<T> smem;
    T *sdata = smem.getPointer();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    sdata[tid] = g_idata[i] + g_idata[i+blockSize];
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; EMUSYNC; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; EMUSYNC; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; EMUSYNC; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; EMUSYNC; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; EMUSYNC; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; EMUSYNC; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)
*/
template <class T, unsigned int blockSize>
__global__ void
FUNC(reduce6)(T *g_idata, T *g_odata, unsigned int n)
{
    SharedMemory<T> smem;
    T *sdata = smem.getPointer();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        sdata[tid] += g_idata[i] + g_idata[i+blockSize];  
        i += gridSize;
    } 
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; EMUSYNC; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; EMUSYNC; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; EMUSYNC; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; EMUSYNC; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; EMUSYNC; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; EMUSYNC; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void 
FUNC(reduce)(int size, int threads, int blocks, 
             int whichKernel, T *d_idata, T *d_odata)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = threads * sizeof(T);

    // choose which of the optimized versions of reduction to launch
    switch (whichKernel)
    {
        case 5:
        switch (threads)
        {
        case 512:
            FUNC(reduce5)<T, 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case 256:
            FUNC(reduce5)<T, 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case 128:
            FUNC(reduce5)<T, 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case 64:
            FUNC(reduce5)<T,  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case 32:
            FUNC(reduce5)<T,  32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case 16:
            FUNC(reduce5)<T,  16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case  8:
            FUNC(reduce5)<T,   8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case  4:
            FUNC(reduce5)<T,   4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case  2:
            FUNC(reduce5)<T,   2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        case  1:
            FUNC(reduce5)<T,   1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata); break;
        }
        break;       
    case 6:
    default:
        switch (threads)
        {
        case 512:
            FUNC(reduce6)<T, 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case 256:
            FUNC(reduce6)<T, 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case 128:
            FUNC(reduce6)<T, 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case 64:
            FUNC(reduce6)<T,  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case 32:
            FUNC(reduce6)<T,  32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case 16:
            FUNC(reduce6)<T,  16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case  8:
            FUNC(reduce6)<T,   8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case  4:
            FUNC(reduce6)<T,   4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case  2:
            FUNC(reduce6)<T,   2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        case  1:
            FUNC(reduce6)<T,   1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
        }
        break;       
    }
}

extern "C"
void FUNC(reduceInt)(int size, int threads, int blocks, 
                     int whichKernel, int *d_idata, int *d_odata)
{
    FUNC(reduce)<int>(size, threads, blocks, whichKernel, d_idata, d_odata);
}

extern "C"
void FUNC(reduceFloat)(int size, int threads, int blocks, 
                       int whichKernel, float *d_idata, float *d_odata)
{
    FUNC(reduce)<float>(size, threads, blocks, whichKernel, d_idata, d_odata);
}

extern "C"
void FUNC(reduceDouble)(int size, int threads, int blocks, 
                        int whichKernel, double *d_idata, double *d_odata)
{
    FUNC(reduce)<double>(size, threads, blocks, whichKernel, d_idata, d_odata);
}

#endif // __GPU_REDUCE_KERNEL_H__
