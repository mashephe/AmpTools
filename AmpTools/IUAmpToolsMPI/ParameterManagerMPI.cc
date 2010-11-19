
#include <mpi.h>
#include <sstream>

#include "IUAmpTools/ParameterManagerMPI.h"
#include "IUAmpTools/MPITag.h"

ParameterManagerMPI::
ParameterManagerMPI( MinuitMinimizationManager& minuitManager,
                    AmplitudeManager* ampManager ) :
ParameterManager( minuitManager, ampManager )
{
  setupMPI();
  
  if( !m_isMaster ){
    
    cerr << "Instance of MinuitMinimizationManager exists on worker node"
    << endl;
    
    assert( false );
  }
}

ParameterManagerMPI::
ParameterManagerMPI( MinuitMinimizationManager& minuitManager,
                    const vector<AmplitudeManager*>& ampManagers ) :
ParameterManager( minuitManager, ampManagers )
{
  setupMPI();
  
  if( !m_isMaster ){
    
    cerr << "Instance of MinuitMinimizationManager exists on worker node"
    << endl;
    
    assert( false );
  }
}

ParameterManagerMPI::
ParameterManagerMPI( AmplitudeManager* ampManager ) :
// feed the ParameterManager constructor a bogus reference,
// this is OK as long as we override calls that use this reference
ParameterManager( *( static_cast< MinuitMinimizationManager* >( NULL ) ),
                 ampManager ),
m_ampManagers( 0 )
{
  setupMPI();
  
  if( m_isMaster ){
    
    cerr << "Master ParameterManager has no MinuitMinimizationManager"
    << endl;
    
    assert( false );
  }
  
  m_ampManagers.push_back( ampManager );
}
ParameterManagerMPI::
ParameterManagerMPI( const vector< AmplitudeManager* >& ampManagers ) :
// feed the ParameterManager constructor a bogus reference,
// this is OK as long as we override calls that use this reference
ParameterManager( *( static_cast< MinuitMinimizationManager* >( NULL ) ),
                 ampManagers ),
m_ampManagers( ampManagers )
{
  setupMPI();
  
  if( m_isMaster ){
    
    cerr << "Master ParameterManager has no MinuitMinimizationManager"
    << endl;
    
    assert( false );
  }
}

ParameterManagerMPI::~ParameterManagerMPI()
{
  
}

void
ParameterManagerMPI::setupMPI()
{
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );
  
  m_isMaster = ( m_rank == 0 );
}

void ParameterManagerMPI::addProductionParameter( const string& ampName ){
  
  // utilize the base class functionality for the master node
  ParameterManager::addProductionParameter( ampName );
  
  // then use member data to hang onto a pointer to the complex number
  m_parMap[ ampName ] = getParPtr( ampName );
}

void ParameterManagerMPI::updateParameters()
{
  
  int nameLength;
  char parName[kMaxNameLength];
  double value[2];
  
  MPI_Status status;
  
  if( m_isMaster ){
    
    for( int proc = 1; proc < m_numProc; ++proc ){
      
      for( map< string, complex< double >* >::iterator 
          parItr = m_parMap.begin();
          parItr != m_parMap.end();
          ++parItr ){
        
        nameLength = parItr->first.length();
        assert( nameLength <= kMaxNameLength );
        strcpy( parName, parItr->first.c_str() );
        
        MPI_Send( &nameLength, 1, MPI_INT, proc,
                 MPITag::kIntSend, MPI_COMM_WORLD );
        MPI_Send( parName, nameLength, MPI_CHAR, proc,
                 MPITag::kCharSend, MPI_COMM_WORLD );
        
        value[0] = parItr->second->real();
        value[1] = parItr->second->imag();
        MPI_Send( value, 2, MPI_DOUBLE, proc,
                 MPITag::kDoubleSend, MPI_COMM_WORLD );
      }
      
      // send zero namelength to indicate to the workers that the
      // parameter send is complete
      nameLength = 0;
      MPI_Send( &nameLength, 1, MPI_INT, proc,
               MPITag::kIntSend, MPI_COMM_WORLD );
    }
  }
  else{
    
    // the workers collect the parameters watching for zero namelength
    // which indicates the send is over
    
    MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
             MPI_COMM_WORLD, &status );
    
    while( nameLength != 0 ){
      
      MPI_Recv( parName, nameLength, MPI_CHAR, 0, MPITag::kCharSend,
               MPI_COMM_WORLD, &status );
      MPI_Recv( value, 2, MPI_DOUBLE, 0, MPITag::kDoubleSend,
               MPI_COMM_WORLD, &status );
      
      ostringstream name;
      for( int i = 0; i < nameLength; ++i ) name << parName[i];
      
      complex< double > parVal( value[0], value[1] );
      
      // check that the name is sane and the reset the value
      map< string, complex< double >* >::iterator parItr = 
      m_parMap.find( name.str() );
      assert( parItr != m_parMap.end() );
      (*(parItr->second)) = parVal;
      
      // fetch the length of the next parameter name which will
      // be zero if we have reached the end of the chain of parameters
      MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
               MPI_COMM_WORLD, &status );
    }
  }
}

