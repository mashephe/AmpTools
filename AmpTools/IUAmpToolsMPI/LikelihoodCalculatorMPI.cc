
#include <mpi.h>
#include <pthread.h>

#include "IUAmpTools/LikelihoodCalculatorMPI.h"
#include "IUAmpTools/LikelihoodManagerMPI.h"

LikelihoodCalculatorMPI::
LikelihoodCalculatorMPI( const AmplitudeManager& ampManager,
			 const NormIntInterface& normInt,
			 DataReader& dataReader,
			 ParameterManagerMPI& parManager ) :
LikelihoodCalculator( ampManager, normInt, dataReader, parManager ),
m_parManager( parManager ),
m_functionEvaluated( false ),
m_thisId( m_idCounter++ )
{
  setupMPI();

  if( !m_isMaster ){

    LikelihoodManagerMPI::registerCalculator( m_thisId, this );
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

void
LikelihoodCalculatorMPI::update( const MISubject* subject ){

  assert( m_isMaster );

  // parameters have changed -- tell workers to update parameters
  // this will trigger expensive calculations in workers
  //
  // note: this is a little inefficient since all instances of the
  // likelihood calculator share the same parameter manager; however,
  // need to be sure parameters are updated before beginning expensive
  // calculations in the update routine of the workers

  // a two element array to hold commands that go to the workers
  // the first element is the id of the likelihood calculator
  // the second element is the flag of the command
  int cmnd[2];
  cmnd[0] = m_thisId;
  
  // tell all of the workers to update parameters
  cmnd[1] = LikelihoodManagerMPI::kUpdateParameters;
  for( int i = 1; i < m_numProc; ++i ){

    MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
  }

  // tell the master to do parameter update
  m_parManager.updateParameters();
}

double
LikelihoodCalculatorMPI::operator()()
{
  assert( m_isMaster );

  MISubject* dummy( NULL );
  if( !m_functionEvaluated ) update( dummy );

  MPI_Status status;

  // a two element array to hold commands that go to the workers
  // the first element is the id of the likelihood calculator
  // the second element is the flag of the command
  int cmnd[2];
  cmnd[0] = m_thisId;
  
  // tell all of the workers to send the partial sums 
  cmnd[1] = LikelihoodManagerMPI::kComputeLikelihood;
  for( int i = 1; i < m_numProc; ++i ){

    MPI_Send( cmnd, 2, MPI_INT, i, MPITag::kIntSend, MPI_COMM_WORLD );
  }

  double lnL = 0;
  double partialSum;

  // collect the sums
  for( int i = 1; i < m_numProc; ++i ){

    MPI_Recv( &partialSum, 1, MPI_DOUBLE, i, MPITag::kDoubleSend,
	      MPI_COMM_WORLD, &status );

    lnL += partialSum;
  }
  
  lnL -= normIntTerm();

  m_functionEvaluated = true;

  return -2 * lnL;
}

void
LikelihoodCalculatorMPI::updateParameters()
{
  assert( !m_isMaster );

  // do the update on the worker nodes
  m_parManager.updateParameters();

  // now use parent class function to do expensive calculations
  MISubject* dummy( NULL );
  LikelihoodCalculator::update( dummy );
}

void
LikelihoodCalculatorMPI::computeLikelihood()
{
  assert( !m_isMaster );

  double lnL = dataTerm();
  MPI_Send( &lnL, 1, MPI_DOUBLE, 0, MPITag::kDoubleSend, MPI_COMM_WORLD );
}

void
LikelihoodCalculatorMPI::setupMPI()
{
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  
  m_isMaster = ( m_rank == 0 );
}

int LikelihoodCalculatorMPI::m_idCounter = kFirstId;
