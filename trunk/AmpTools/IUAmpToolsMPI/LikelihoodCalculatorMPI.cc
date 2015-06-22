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

#include <mpi.h>
#include <pthread.h>

#include "IUAmpToolsMPI/LikelihoodCalculatorMPI.h"
#include "IUAmpToolsMPI/LikelihoodManagerMPI.h"

LikelihoodCalculatorMPI::
LikelihoodCalculatorMPI( const IntensityManager& intenManager,
                        const NormIntInterface& normInt,
                        DataReader* dataReader,
                        DataReader* bkgReader,
                        ParameterManagerMPI& parManager ) :
LikelihoodCalculator( intenManager, normInt, dataReader, bkgReader, parManager ),
m_intenManager( intenManager ),
m_parManager( parManager ),
m_thisId( m_idCounter++ ),
m_firstPass( true )
{
  setupMPI();
  
  LikelihoodManagerMPI::registerCalculator( m_thisId, this );
  
  if( !m_isMaster ){
    
    // check back in with the master after registration
    MPI_Send( &m_thisId, 1, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD );
  }
  else{
    
    // wait for the workers to check in before proceeding -- if this is 
    // not done, the immediate calls to operator()() could behave
    // unexpectedly
    
    int id;
    MPI_Status status;
    
    for( int i = 1; i < m_numProc; ++i ){
      
      MPI_Recv( &id, 1, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD, &status );
      
      // ids should match
      assert( m_thisId == id );
    }
  }  
}

LikelihoodCalculatorMPI::~LikelihoodCalculatorMPI(){
  
  // have the master break all of the workers out of their
  // deliverLikelihood loops
  if( m_isMaster && m_thisId == kFirstId ){
    
    // a two element array to hold commands that go to the workers
    // the first element is the id of the likelihood calculator
    // the second element is the flag of the command
    int cmnd[2];
    cmnd[0] = m_thisId;
    
    // break the likelihood manager out of its loop on the workers
    cmnd[1] = LikelihoodManagerMPI::kExit;
    for( int i = 1; i < m_numProc; ++i ){
      
      MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
    }
  }
}

double
LikelihoodCalculatorMPI::operator()()
{
  assert( m_isMaster );
  
  MPI_Status status;
  
  // a two element array to hold commands that go to the workers
  // the first element is the id of the likelihood calculator
  // the second element is the flag of the command
  int cmnd[2];
  cmnd[0] = m_thisId;
  
  // tell all of the workers to update parameters
  // note: this is a little inefficient since all instances of the
  // likelihood calculator share the same parameter manager -- this
  // cause a little extra overhead in multiple final-state fits
  cmnd[1] = LikelihoodManagerMPI::kUpdateParameters;
  for( int i = 1; i < m_numProc; ++i ){
    
    MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
  }
  
  // tell the master to do parameter update
  m_parManager.updateParameters();
  
  // tell all of the workers to send the partial sums 
  cmnd[1] = LikelihoodManagerMPI::kComputeLikelihood;
  for( int i = 1; i < m_numProc; ++i ){
    
    MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
  }
  
  double lnL = 0;
  double sumBkgWeights = 0;
  double numBkgEvents  = 0;
  double numDataEvents = 0;
  double data[4];
  
  // collect the sums
  for( int i = 1; i < m_numProc; ++i ){
    
    MPI_Recv( &data, 4, MPI_DOUBLE, i, MPITag::kDoubleSend,
             MPI_COMM_WORLD, &status );
    
    lnL           += data[0];
    sumBkgWeights += data[1];
    numBkgEvents  += data[2];
    numDataEvents += data[3];
  }
  
  setSumBkgWeights( sumBkgWeights );
  setNumBkgEvents ( numBkgEvents  );
  setNumDataEvents( numDataEvents );
  
  // if we have an amplitude with a free parameter, the call to normIntTerm()
  // will trigger recomputation of NI's -- we need to put the workers in the
  // loop to send the recomputed NI's back to the master
  if( m_intenManager.hasTermWithFreeParam() || m_firstPass ){
    
    cmnd[1] = LikelihoodManagerMPI::kComputeIntegrals;
    for( int i = 1; i < m_numProc; ++i ){
      
      MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
    }
  }
  
  // this call will utilize the NormIntInterface on the master which
  // ultimately gets handled by the instance of NormIntInterfaceMPI
  // that is passed into the constructor of this class
  lnL -= normIntTerm();
  
  m_firstPass = false;
  
  return -2 * lnL;
}

void
LikelihoodCalculatorMPI::updateParameters()
{
  assert( !m_isMaster );
  
  // do the update on the worker nodes
  m_parManager.updateParameters();
}

void
LikelihoodCalculatorMPI::updateAmpParameter()
{
  assert( !m_isMaster );
  
  // do the update on the worker nodes
  m_parManager.updateAmpParameter();
}

void
LikelihoodCalculatorMPI::computeLikelihood()
{
  assert( !m_isMaster );
  
  double data[4];
  
  data[0] = dataTerm();
  data[1] = sumBkgWeights();
  data[2] = numBkgEvents();
  data[3] = numDataEvents();
  
  MPI_Send( &data, 4, MPI_DOUBLE, 0, MPITag::kDoubleSend, MPI_COMM_WORLD );
}

void
LikelihoodCalculatorMPI::setupMPI()
{
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  
  m_isMaster = ( m_rank == 0 );
}

int LikelihoodCalculatorMPI::m_idCounter = kFirstId;
