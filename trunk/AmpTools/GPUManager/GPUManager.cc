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

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "cuda_runtime.h"

#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/AmpVecs.h"

#include "GPUManager/GPUKernel.h"
#include "GPUManager/GPUManager.h"

bool GPUManager::m_cudaDisplay = false;

template <class T>
void reduce(int size, int threads, int blocks,
                 int whichKernel, T *d_idata, T *d_odata);

GPUManager::GPUManager() : 
  m_ampCalcOnly( false )
{
	m_iNEvents=0;
	m_iNTrueEvents=0;

	m_iNAmps=0;
	m_iNAmpsH=0;
	m_iAmpArrSize=0;

	m_iEventArrSize=0;
	m_iVArrSize=0;
	
	//Host Arrays
	m_pfAmpRe=0;
	m_pfAmpIm=0;
	m_pfVRe=0;
	m_pfVIm=0;
	m_pfRes=0;
	
	//Device Arrays	
	m_pfDevData=0;
	m_pcDevCalcAmp=0;
  m_piDevPerm=0;
	m_pfDevAmpRe=0;
	m_pfDevAmpIm=0;
  m_pfDevWeights=0;
	m_pfDevVRe=0;
	m_pfDevVIm=0;
	m_pfDevRes=0;
	m_pfDevREDUCE=0;
	
	//CUDA Thread and Grid sizes
	m_iDimGridX=0;
	m_iDimGridY=0;
	m_iDimThreadX=0;
	m_iDimThreadY=0;
  
  int thisDevice = 0;

  if( !m_cudaDisplay )
    cout<<"\n################### CUDA DEVICE ##################\n";    
  
#ifdef USE_MPI
  
  // Note that a better algorithm would be to utilize the jobs "local
  // rank" on the machine instead of global rank -- this needs development.
  // The obvious problem wiht the technique below is that, e.g., in a 
  // two-GPU per node cluster, all even rank jobs will use device zero.
  // If two even rank jobs land on the same node device zero will be
  // overloaded and device 1 will be unused.
  
  int rank, devs;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  cudaGetDeviceCount(&devs);
  thisDevice = rank % devs;
  
  if( !m_cudaDisplay ) {
    cout << "\nParallel GPU configuration requested." << endl;
    cout << "\nNumber of CUDA devices available on this node:  " << devs << endl;
    cout << "\nMPI process " << rank << " is using device " << thisDevice << endl;
  }
#endif
  
	///////CUDA INITIALIZATION
  cudaSetDevice( thisDevice );
  cudaDeviceProp devProp;
  cudaGetDeviceProperties( &devProp, thisDevice );
  
	if( ! m_cudaDisplay ){
    
    cout<<"Current GPU Properites:\n";
    cout<<"\t Name: "<<devProp.name<<endl; 
    cout<<"\t Total global memory: "<<devProp.totalGlobalMem/((float)1024*1024)<<" MB"<<endl; 
    cout<<"\t Number of cores: " << 8*devProp.multiProcessorCount << endl;
    cout<<"\t Rev.: "<<devProp.major<<"."<<devProp.minor<<endl;
    cout<<"\t Precision (size of GDouble): " << sizeof(GDouble) << " bytes" << endl;	
    cout<<"##################################################\n\n";
    ///////END OF CUDA INITIALIZATION
    m_cudaDisplay = true;
  }
  
  if( ( devProp.major == 1 ) && devProp.minor < 3 ){
    
    // double precision operations need 1.3 hardware or higher
    assert( sizeof( GDouble ) <= 4 );
  }
  
	cudaError_t cerrKernel=cudaGetLastError();
	if( cerrKernel != cudaSuccess  ){
    
		cout<<"\n\nDEVICE INIT ERROR: "<< cudaGetErrorString(cerrKernel) << endl;  
    assert( false );
  }
}

GPUManager::GPUManager( const AmpVecs& a )
{
	GPUManager();
	init( a );
}

GPUManager::~GPUManager()
{
	clearAll();
}


// Initialization routines:

void 
GPUManager::init( const AmpVecs& a, bool ampCalcOnly )
{
  clearAll();
  
  m_ampCalcOnly = ampCalcOnly;
  
  // copy over some info from the AmpVecs object for array dimensions
	m_iNTrueEvents = a.m_iNTrueEvents;
	m_iNEvents     = a.m_iNEvents;
  m_iNParticles  = a.m_iNParticles;
	m_iNAmps       = a.m_iNAmps;
  
  // the rest of the data are derived:
	m_iEventArrSize     = sizeof(GDouble) * m_iNEvents;
  m_iTrueEventArrSize = sizeof(GDouble) * m_iNTrueEvents;

	// size of upper half of AiAj* matrix
  m_iNAmpsH = m_iNAmps * ( m_iNAmps + 1 ) / 2;
  
  // size needed to store amplitudes for each event
	m_iAmpArrSize     = sizeof(GDouble) * m_iNEvents * m_iNAmpsH;

	// size of upper half of ViVj* matrix
  m_iVArrSize = sizeof(GDouble) * m_iNAmpsH;
  
  // save memory by not allocating it if we are only using the GPU
  // to do amplitude calcuations (as in normalization integrals)
  if( !m_ampCalcOnly ){
  
    // host memory needed for intensity calculation
    cudaMallocHost( (void**)&m_pfAmpRe , m_iAmpArrSize   );
    cudaMallocHost( (void**)&m_pfAmpIm , m_iAmpArrSize   );
    cudaMallocHost( (void**)&m_pfVRe   , m_iVArrSize     );
    cudaMallocHost( (void**)&m_pfVIm   , m_iVArrSize     );
    cudaMallocHost( (void**)&m_pfRes   , m_iEventArrSize );

    // device memory needed for intensity calculation and sum
    cudaMalloc( (void**)&m_pfDevAmpRe   , m_iAmpArrSize                       );
    cudaMalloc( (void**)&m_pfDevAmpIm   , m_iAmpArrSize                       );
    cudaMalloc( (void**)&m_pfDevVRe     , m_iVArrSize                         );
    cudaMalloc( (void**)&m_pfDevVIm     , m_iVArrSize                         );
    cudaMalloc( (void**)&m_pfDevWeights , m_iEventArrSize                     );
    cudaMalloc( (void**)&m_pfDevRes     , m_iEventArrSize                     );
    cudaMalloc( (void**)&m_pfDevREDUCE  , m_iEventArrSize                     );
  }
  
  
  // allocate device memory needed for amplitude calculations
  cudaMalloc( (void**)&m_pfDevData    , 4 * m_iNParticles * m_iEventArrSize );
  cudaMalloc( (void**)&m_pcDevCalcAmp , 2 * m_iEventArrSize                 );
  cudaMalloc( (void**)&m_piDevPerm    , m_iNParticles * sizeof( int )       );
	
	cudaError_t cerrKernel = cudaGetLastError();
	if( cerrKernel!= cudaSuccess  ){
        
    // an error here probably means we don't have enough
    // RAM in the GPU for the number of events and amplitudes
    // we plan to work with
    
		cout << "\n\nMEMORY ALLOCATION ERROR: " 
         << cudaGetErrorString( cerrKernel ) << endl;
    
    assert( false );
  }
    
  cout << "GPU memory allocated for " << m_iNAmps << " amplitudes and "
       << m_iNEvents << " events (" << m_iNTrueEvents << " actual events)"
       << endl;
  
  // check for errors
  cerrKernel = cudaGetLastError();
	if( cerrKernel!= cudaSuccess  ){
    
		cout << "\nError initializing constant memory: " 
         << cudaGetErrorString( cerrKernel ) 
         << endl;
    assert( false );
  }
  
	//CUDA Dims
	calcCUDADims();
}


void
GPUManager::copyDataToGPU( const AmpVecs& a )
{		
  
  // make sure AmpVecs has been loaded with data
  assert( a.m_pdData );
    
  // copy the data into the device
  cudaMemcpy( m_pfDevData, a.m_pdData, 4 * m_iNParticles * m_iEventArrSize, 
              cudaMemcpyHostToDevice );
   
  // check for errors
 	cudaError_t cerrKernel = cudaGetLastError();
	if( cerrKernel!= cudaSuccess  ){
    
		cout << "\nError copying data to device: " << cudaGetErrorString( cerrKernel ) 
         << endl;
    assert( false );
  }
}


void
GPUManager::copyAmpsToGPU( const AmpVecs& a )
{
	if(!m_pfAmpRe) {
    
		cout << "GPUManager::InitAmps is called without initalization or this\n" 
         << "instance of GPUManager is for amplitude calculation only." << endl;
    assert( false );
	}
	   
	unsigned int i,j,iEvent;
	for( iEvent = 0; iEvent < m_iNTrueEvents; iEvent++ )
	{
		//Saving only the upper half of the AiAj*
		for( i = 0; i < m_iNAmps; i++ )
			for( j = 0; j <= i; j++ )
			{
				m_pfAmpRe[iEvent+m_iNEvents*(i*(i+1)/2+j)] = 
          a.m_pdAmps[2*m_iNEvents*i+2*iEvent]   * a.m_pdAmps[2*m_iNEvents*j+2*iEvent] + 
          a.m_pdAmps[2*m_iNEvents*i+2*iEvent+1] * a.m_pdAmps[2*m_iNEvents*j+2*iEvent+1];
				
				m_pfAmpIm[iEvent+m_iNEvents*(i*(i+1)/2+j)] = 
          -a.m_pdAmps[2*m_iNEvents*i+2*iEvent]   * a.m_pdAmps[2*m_iNEvents*j+2*iEvent+1] + 
           a.m_pdAmps[2*m_iNEvents*i+2*iEvent+1] * a.m_pdAmps[2*m_iNEvents*j+2*iEvent];
				
				//Doubling the off-diagonal elements to sum over only upper triangle
				if(j !=i )
				{
					m_pfAmpRe[iEvent+m_iNEvents*(i*(i+1)/2+j)]*=2.;
					m_pfAmpIm[iEvent+m_iNEvents*(i*(i+1)/2+j)]*=2.;
				}	
			}
	}
  
	//Now padding the upper half to make sure there are no nans in log
	for( ; iEvent < m_iNEvents; iEvent++ ) {
		for( i = 0; i < m_iNAmps; i++ ) {
			for( j = 0; j <= i; j++ ) {
				m_pfAmpRe[iEvent+m_iNEvents*(i*(i+1)/2+j)] = 1.;
				m_pfAmpIm[iEvent+m_iNEvents*(i*(i+1)/2+j)] = 0.;
      }
    }
  }

  
  /* // useful block for debugging:
  for( iEvent = 0; iEvent < m_iNEvents; iEvent++ ){
    
    cout << "Event " << iEvent << endl;
    for( i = 0; i < m_iNAmps; i++ ){
      
      cout << "Amp " << i << ":\t";
			for( j = 0; j <= i; j++ ){
			
        cout << "(" <<
        m_pfAmpRe[iEvent+m_iNEvents*(i*(i+1)/2+j)]
        << ", " <<
				m_pfAmpIm[iEvent+m_iNEvents*(i*(i+1)/2+j)]
        << ")\t";
      }
      cout << endl;
    }
  }
   */
   
	cudaMemcpy( m_pfDevAmpRe, m_pfAmpRe, m_iAmpArrSize, cudaMemcpyHostToDevice );
	cudaMemcpy( m_pfDevAmpIm, m_pfAmpIm, m_iAmpArrSize, cudaMemcpyHostToDevice );

  // copy the weights to the GPU
  cudaMemcpy( m_pfDevWeights, a.m_pdWeights, m_iEventArrSize, cudaMemcpyHostToDevice );
  
  cudaError_t cerrKernel = cudaGetLastError();
	if( cerrKernel!= cudaSuccess  ){
    
		cout << "\nError copying amplitudes or weights to GPU: " 
    << cudaGetErrorString(cerrKernel) << endl;
    assert( false );
  }
  
}

void
GPUManager::setParamPtrs(  const vector< const complex< double >* >& pvpFitPars )
{
	m_vpFitPars = pvpFitPars;
}

void
GPUManager::setCoherenceMatrix( const vector< vector < bool > >& cohMtx ){
  
  m_vbSumCoherently = cohMtx;
}

void 
GPUManager::calcAmplitudeAll( const Amplitude* amp, GDouble* pcResAmp, 
                              const vector< vector< int > >* pvPermutations )
{
  
	dim3 dimBlock( m_iDimThreadX, m_iDimThreadY );
	dim3 dimGrid( m_iDimGridX, m_iDimGridY );

  // do the computation for all events for each permutation in the
  // vector of permunations
  vector< vector< int > >::const_iterator permItr = pvPermutations->begin();

  // if this is not true, AmplitudeManager hasn't been setup properly
  assert( permItr->size() == m_iNParticles );
  
  int permOffset = 0;
  for( ; permItr != pvPermutations->end(); ++permItr ){
    
    // copy the permutation to global memory
    cudaMemcpy( m_piDevPerm, &((*permItr)[0]), m_iNParticles * sizeof( int ),
                cudaMemcpyHostToDevice );

    // check for errors
    cudaError_t cerrPreKernel = cudaGetLastError();
    if( cerrPreKernel != cudaSuccess  ){
      cout << "\nError copying permutation to GPU: " 
           << cudaGetErrorString( cerrPreKernel ) << endl;
      assert( false );
    }
    
    // calculate amplitude factor for all events --
    // casting amp array to WCUComplex for 8 or 16 bit write 
    // operation of both real and complex parts at once
    
    amp->calcAmplitudeGPU( dimGrid, dimBlock, m_pfDevData, 
                          (WCUComplex*)m_pcDevCalcAmp, 
                          m_piDevPerm, m_iNParticles, m_iNEvents,
                          *permItr );

    cudaThreadSynchronize();    

    // check to be sure kernel execution was OK
    cudaError_t cerrKernel=cudaGetLastError();
    if( cerrKernel!= cudaSuccess  ){
      
      cout << "\nKERNEL LAUNCH ERROR [" << amp->name() << "]: " 
           << cudaGetErrorString( cerrKernel ) << endl;
      assert( false );
    }
    
    // now copy the result out of the GPU into the correct place in the
    // pcResAmp array for this particular permutation
    cudaMemcpy( &(pcResAmp[permOffset]), m_pcDevCalcAmp, 2 * m_iEventArrSize, 
               cudaMemcpyDeviceToHost);
    
    cerrKernel = cudaGetLastError();
    if( cerrKernel!= cudaSuccess  ){
      
      cout << "\nError copying amplitudes from GPU [" << amp->name() << "]: " 
      << cudaGetErrorString(cerrKernel) << endl;
      assert( false );
    }
    
    // increment the offset so that we place the computation for the
    // next permutation after the previous in pcResAmp
    permOffset += 2 * m_iNEvents;
  }    
}

double GPUManager::calcSumLogIntensity()
{
  
  // be sure memory has been allocated for intensity computation
  assert( !m_ampCalcOnly );
  
	unsigned int i,j;
  
  // precompute the real and imaginary parts of ViVj* and copy to 
  // GPU global memory
	complex< double > cdFij;
	for( i = 0; i< m_iNAmps; i++) {
		for( j = 0; j <= i; j++ ) {
      
			cdFij = (*(m_vpFitPars[i])) * conj(*(m_vpFitPars[j]));
      
      // here is the transition from double -> GDouble
			m_pfVRe[i*(i+1)/2+j] = 
        ( m_vbSumCoherently[i][j] ? static_cast< GDouble >( cdFij.real() ) : 0 );
			m_pfVIm[i*(i+1)/2+j] = 
        ( m_vbSumCoherently[i][j] ? static_cast< GDouble >( cdFij.imag() ) : 0 );
    }
  }
	  
  //Init global memory on GPU
	cudaMemcpyToSymbol( "da_pfDevVRe" , m_pfVRe     , m_iVArrSize );
	cudaMemcpyToSymbol( "da_pfDevVIm" , m_pfVIm     , m_iVArrSize );
	cudaMemcpyToSymbol( "da_iNAmpsH"  , &m_iNAmpsH  , sizeof(int) );
	cudaMemcpyToSymbol( "da_iNEvents" , &m_iNEvents , sizeof(int) );

  cudaError_t cerrPreKernel=cudaGetLastError();
  if( cerrPreKernel != cudaSuccess  ){
    cout << "\nError copying weights or parameters to GPU: " 
    << cudaGetErrorString( cerrPreKernel ) << endl;
    assert( false );
  }
  
  // compute the intensities
	dim3 dimBlock( m_iDimThreadX, m_iDimThreadY );
	dim3 dimGrid( m_iDimGridX, m_iDimGridY );
	GPU_ExecAmpKernel( dimGrid, dimBlock, m_pfDevAmpRe, m_pfDevAmpIm, 
                     m_pfDevWeights, m_pfDevRes );  
  
	// Now the summation of the results	-- do this on the CPU for small
  // numbers of events or cases where double precision GPU is not enabled
	double dGPUResult = 0;
	if( m_iNTrueEvents <= m_iNBlocks || sizeof( GDouble ) <= 4 )
	{
		cudaMemcpy(m_pfRes,m_pfDevRes,m_iTrueEventArrSize,cudaMemcpyDeviceToHost);
		for(i=0; i<m_iNTrueEvents; i++)
			dGPUResult += m_pfRes[i];
	}
	else
	{
		cudaThreadSynchronize();		
		//Zeroing out the padding as not to alter the results after the reduction
		cudaMemset(m_pfDevRes+m_iNTrueEvents,0,sizeof(GDouble)*(m_iNEvents-m_iNTrueEvents));
		int whichKernel = 6;   		
		// execute the kernel to sum partial sums from each block on CPU  
		reduce<GDouble>(m_iNEvents, m_iNThreads, m_iNBlocks, whichKernel, m_pfDevRes, m_pfDevREDUCE);
		// copy result from device to host
		cudaMemcpy( m_pfRes, m_pfDevREDUCE, m_iNBlocks*sizeof(GDouble), cudaMemcpyDeviceToHost);
		for(i=0; i<m_iNBlocks; i++) 
			dGPUResult += m_pfRes[i];
	}
	return dGPUResult;
}

// Methods to clear memory:

void GPUManager::clearAll()
{
	clearAmpCalc();
  
	if( !m_ampCalcOnly ) clearLikeCalc();
}

void GPUManager::clearAmpCalc()
{
	m_iNParticles=0;
	m_iNEvents=0;
	m_iNTrueEvents=0;
	m_iEventArrSize=0;
	m_iTrueEventArrSize=0;
		
	//Device Memory	
	if(m_pfDevData)
		cudaFree(m_pfDevData);
	m_pfDevData=0;
	
	if(m_pcDevCalcAmp)
		cudaFree(m_pcDevCalcAmp);
	m_pcDevCalcAmp=0;
  
  if(m_piDevPerm)
    cudaFree(m_piDevPerm);
  m_piDevPerm=0;
	
}

void GPUManager::clearLikeCalc()
{
	m_iNEvents=0;
	m_iNTrueEvents=0;
	m_iNAmps=0;
	m_iNAmpsH=0;

	m_iEventArrSize=0;
  m_iTrueEventArrSize=0;
	m_iAmpArrSize=0;
	m_iVArrSize=0;
	
	//Host Memory
	
	//Allocated pointers
	if(m_pfAmpRe)
		cudaFreeHost(m_pfAmpRe);
	m_pfAmpRe=0;
	
	if(m_pfAmpIm)
		cudaFreeHost(m_pfAmpIm);
	m_pfAmpIm=0;
	
	if(m_pfVRe)
		cudaFreeHost(m_pfVRe);
	m_pfVRe=0;
	
	if(m_pfVIm)
		cudaFreeHost(m_pfVIm);
	m_pfVIm=0;
	  
	if(m_pfRes)
		cudaFreeHost(m_pfRes);
	m_pfRes=0;
	
	//Device Memory	
		
	if(m_pfDevAmpRe)
		cudaFree(m_pfDevAmpRe);
	m_pfDevAmpRe=0;
	
	if(m_pfDevAmpIm)
		cudaFree(m_pfDevAmpIm);
	m_pfDevAmpIm=0;
	
	if(m_pfDevVRe)
		cudaFree(m_pfDevVRe);
	m_pfDevVRe=0;
	
	if(m_pfDevVIm)
		cudaFree(m_pfDevVIm);
	m_pfDevVIm=0;
	
  if(m_pfDevWeights)
    cudaFree(m_pfDevWeights);
  m_pfDevWeights=0;
  
	if(m_pfDevRes)
		cudaFree(m_pfDevRes);
	m_pfDevRes=0;
	
	if(m_pfDevREDUCE)
		cudaFree(m_pfDevREDUCE);
	m_pfDevREDUCE=0;
	
	//CUDA Thread and Grid sizes
	m_iDimGridX=0;
	m_iDimGridY=0;
	m_iDimThreadX=0;
	m_iDimThreadY=0;
}

// Internal utilities:
void GPUManager::calcCUDADims()
{
	if(m_iNEvents<1)
		return;
	
	m_iDimThreadX=GPU_BLOCK_SIZE_X;
	m_iDimThreadY=GPU_BLOCK_SIZE_Y;	
	
	unsigned int iBlockSizeSq=GPU_BLOCK_SIZE_SQ;
	unsigned int iNBlocks=m_iNEvents/iBlockSizeSq;
	if(iNBlocks<=1)
	{
		m_iDimGridX=1;
		m_iDimGridY=1;
	}
	else
	{
		unsigned int iDivLo=1,iDivHi=iNBlocks;
		for(iDivLo=static_cast<int>(sqrt(iNBlocks));iDivLo>=1;iDivLo--)
		{
			iDivHi=iNBlocks/iDivLo;
			if(iDivLo*iDivHi==iNBlocks)
				break;
		}
		m_iDimGridX=iDivLo;
		m_iDimGridY=iDivHi;
	}

	// cout<<"\tThread dimensions:  ("<<m_iDimThreadX<<","<<m_iDimThreadY<<")\n";
	// cout<<"\tGrid dimensions:  ("<<m_iDimGridX<<","<<m_iDimGridY<<")\n";
	
	//Reduction Parameters
	unsigned int maxThreads = 256;  // number of threads per block
	unsigned int maxBlocks = 256;		
	
  if (m_iNEvents == 1) 
		m_iNThreads = 1;
	else
		m_iNThreads = (m_iNEvents < maxThreads*2) ? m_iNEvents / 2 : maxThreads;
	
	m_iNBlocks = m_iNEvents / (m_iNThreads * 2);	
	m_iNBlocks = min(maxBlocks, m_iNBlocks);
	
	// cout<<"Reduction:\n";
	// cout<<"\tNumber of threads:  "<<m_iNThreads<<"\n";
  // cout<<"\tNumber of blocks:   "<<m_iNBlocks<<"\n\n\n"<<flush;	
}

