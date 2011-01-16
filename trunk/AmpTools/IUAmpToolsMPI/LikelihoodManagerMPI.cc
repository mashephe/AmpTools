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
#include <map>

#include <mpi.h>

#include "IUAmpToolsMPI/LikelihoodManagerMPI.h"
#include "IUAmpToolsMPI/LikelihoodCalculatorMPI.h"
#include "IUAmpToolsMPI/MPITag.h"

void
LikelihoodManagerMPI::registerCalculator( int id, 
                                         LikelihoodCalculatorMPI* calc )
{
  m_calcMap[id] = calc;
}

void
LikelihoodManagerMPI::deliverLikelihood()
{
  if( !m_mpiSetup ) setupMPI();
  
  if( m_isMaster ){
    
    cerr << "ERROR:  deliverLikelihood() should not run on master node" 
    << endl;
    assert( false );
  }
  
  MPI_Status status;
  
  int cmnd[2];
  int* id = &(cmnd[0]);
  int* fitFlag = &(cmnd[1]);
  
  LikelihoodCalculatorMPI* likCalc;
  map< int, LikelihoodCalculatorMPI* >::iterator mapItr;
  
  MPI_Recv( cmnd, 2, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD, &status );
  
  while( *fitFlag != kExit ){
    
    mapItr = m_calcMap.find( *id );
    // we must have a matching likelihood calculator to proceed
    assert( mapItr != m_calcMap.end() );
    likCalc = mapItr->second;
    
    switch( *fitFlag ){
        
      case kUpdateParameters:
        
        likCalc->updateParameters();
        break;
  
      case kUpdateAmpParameter:
        
        likCalc->updateAmpParameter();
        break;
        
      case kComputeIntegrals:
        
        // this actually has a return value, but serves the purpose
        // of triggering a NI recalculation, if necessary, on the workers
        likCalc->normIntTerm();
        break;
        
      case kComputeLikelihood:
        
        likCalc->computeLikelihood();
        break;
        
      default:
        
        cerr << "Unknown command flag!" << endl;
        assert( false );
    }
    
    MPI_Recv( &cmnd, 2, MPI_INT, 0, MPITag::kIntSend, 
             MPI_COMM_WORLD, &status );
  }
  
  assert( *fitFlag == kExit );
}

void
LikelihoodManagerMPI::broadcastToFirst( FitCommand command ){
  
  if( !m_mpiSetup ) setupMPI();
  
  // this broadcasts a particular command to the first registered
  // calculator -- note that there can be multiple likelihood
  // calculators per node as there is one likelihood calculator
  // for every reaction
  
  // should only be called on the master:
  assert( m_isMaster );

  int cmnd[2];
  cmnd[1] = command;

  map< int, LikelihoodCalculatorMPI* >::iterator mapItr = m_calcMap.begin();
  
  // if this is false then there are no registered calculators
  assert( mapItr != m_calcMap.end() );
  
  cmnd[0] = mapItr->first;    
   
  // this will send a command to just the first registered calculator
  for( int i = 1; i < m_numProc; ++i ){
      
    MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
  }
}

void
LikelihoodManagerMPI::setupMPI()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  m_isMaster = ( rank == 0 );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  m_mpiSetup = true;
}

bool LikelihoodManagerMPI::m_mpiSetup = false;
bool LikelihoodManagerMPI::m_isMaster = false;
int LikelihoodManagerMPI::m_numProc = 0;
map< int, LikelihoodCalculatorMPI* > LikelihoodManagerMPI::m_calcMap;
