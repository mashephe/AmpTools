
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
LikelihoodManagerMPI::setupMPI()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  m_isMaster = ( rank == 0 );
  m_mpiSetup = true;
}

bool LikelihoodManagerMPI::m_mpiSetup = false;
bool LikelihoodManagerMPI::m_isMaster = false;
map< int, LikelihoodCalculatorMPI* > LikelihoodManagerMPI::m_calcMap;
