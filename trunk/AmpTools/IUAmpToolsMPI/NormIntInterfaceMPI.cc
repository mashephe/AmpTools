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

#include <iostream>

#include <mpi.h>
#include <cstring>

#include "IUAmpTools/IntensityManager.h"

#include "IUAmpToolsMPI/NormIntInterfaceMPI.h"
#include "IUAmpToolsMPI/MPITag.h"

using namespace std;

NormIntInterfaceMPI::NormIntInterfaceMPI( DataReader* genMCData, 
                                          DataReader* accMCData, 
                                          const IntensityManager& intenManager ):
NormIntInterface( genMCData, accMCData, intenManager )
{
  setupMPI();  
}

NormIntInterfaceMPI::NormIntInterfaceMPI( const string& normIntFile ) :
NormIntInterface( normIntFile )
{}

NormIntInterfaceMPI::~NormIntInterfaceMPI() {
  
}

complex< double > 
NormIntInterfaceMPI::normInt( string amp, string conjAmp, bool forceUseCache ) const {
  
  // in the case that we have a free parameter, recompute the NI's in parallel
  // note that evaluation order is important; if the NI interface comes from a file
  // the the ampManager pointer will be NULL. We rely on hasAccessToMC to short
  // circuit the evaluation to avoid a segmentation fault.
  if( hasAccessToMC() && intenManager()->hasTermWithFreeParam() && !forceUseCache )
    forceCacheUpdate( true );
  
  // then we can use the parent class to return the value from the 
  // updated cache
  return NormIntInterface::normInt( amp, conjAmp, true );
}

void NormIntInterfaceMPI::forceCacheUpdate( bool normIntOnly ) const {
  
  if( !m_isMaster ) NormIntInterface::forceCacheUpdate( normIntOnly );
  
  if( !normIntOnly ) sumIntegrals( kAmpInt );
  sumIntegrals( kNormInt );
}

void
NormIntInterfaceMPI::setupMPI()
{
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  
  MPI_Status status;
  
  m_isMaster = ( m_rank == 0 );
  
  int totalGenEvents = 0;
  int totalAccEvents = 0;
    
  if( m_isMaster ){
    
    for( int i = 1; i < m_numProc; ++i ){
      
      int thisEvents;

      // trigger sending of events from workers -- data is irrelevant
      MPI_Send( &thisEvents, 1, MPI_INT, i, MPITag::kAcknowledge,
               MPI_COMM_WORLD );
      
      // now receive actual data
      MPI_Recv( &thisEvents, 1, MPI_INT, i, MPITag::kIntSend,
               MPI_COMM_WORLD, &status );
      totalGenEvents += thisEvents;
      
      MPI_Recv( &thisEvents, 1, MPI_INT, i, MPITag::kIntSend,
               MPI_COMM_WORLD, &status );
      totalAccEvents += thisEvents;
      
      // send acknowledgment 
      MPI_Send( &thisEvents, 1, MPI_INT, i, MPITag::kAcknowledge,
               MPI_COMM_WORLD );
    }
    
    setGenEvents( totalGenEvents );
    setAccEvents( totalAccEvents );
  }
  else{
    
    int thisEvents;

    // if we are not the master, send generated and accepted events
    // to master and wait for acknowledge

    // workers can beat the master to this point in the code if the master is
    // still distributing data -- need to pause here and wait for master
    // to signal that it is ready to accept numbers of events

    // data is irrelevant for this receive
    MPI_Recv( &thisEvents, 1, MPI_INT, 0, MPITag::kAcknowledge, MPI_COMM_WORLD,
              &status );

    thisEvents = numGenEvents();
    MPI_Send( &thisEvents, 1, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD );
    
    thisEvents = numAccEvents();
    MPI_Send( &thisEvents, 1, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD );
    
    MPI_Recv( &thisEvents, 1, MPI_INT, 0, MPITag::kAcknowledge, MPI_COMM_WORLD,
             &status );
  }
}

void
NormIntInterfaceMPI::sumIntegrals( IntType type ) const
{  
  // first collect integrals from all nodes -- const cast is bad
  // but we're going to reset the memory in the NormIntInterface
  // at the end of this routine.  This avoids copying.
  
  GDouble* integrals =
    const_cast< GDouble* >( type == kNormInt ? normIntMatrix() : ampIntMatrix() );
  
  // scale up the integrals
  for( int i = 0; i < cacheSize(); ++i ) integrals[i] *= numGenEvents();
  
  GDouble* result = new GDouble[cacheSize()];
  
  // zero out the result array
  if( m_isMaster ) memset( integrals, 0, cacheSize() * sizeof( GDouble ) );
  
  // at this point, the master should be holding an array full of zeroes
  // and other nodes will hold integrals for their subsets of data
  
  // sum over all nodes and distribute results to each
  MPI_Datatype t =
    ( sizeof( GDouble ) == sizeof( double ) ? MPI_DOUBLE : MPI_FLOAT );  
  MPI_Allreduce( integrals, result, cacheSize(), t, MPI_SUM, MPI_COMM_WORLD );

  // now broadcast the total number of events from the master to the
  // workers so that they may renormalize the sum properly
  int totalEvents = numGenEvents();
  MPI_Bcast( &totalEvents, 1, MPI_INT, 0, MPI_COMM_WORLD );
  
  // and renormalize the sum
  for( int i = 0; i < cacheSize(); ++i ) result[i] /= totalEvents;
  
  // reset the memory in the parent class
  if( type == kNormInt ){
    
    setNormIntMatrix( result );
  }
  else{
    
    setAmpIntMatrix( result );
  }
  
  delete[] result;
}
