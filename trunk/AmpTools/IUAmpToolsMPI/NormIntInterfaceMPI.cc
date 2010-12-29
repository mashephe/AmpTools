
#include <iostream>

#include <mpi.h>

#include "IUAmpToolsMPI/NormIntInterfaceMPI.h"
#include "IUAmpToolsMPI/MPITag.h"

using namespace std;

NormIntInterfaceMPI::NormIntInterfaceMPI( DataReader* genMCData, 
                                         DataReader* accMCData, 
                                         const AmplitudeManager& ampManager ):
NormIntInterface( genMCData, accMCData, ampManager )
{
  
  setupMPI();
  
  sumIntegrals( kAmpInt );
  sumIntegrals( kNormInt );
}

NormIntInterfaceMPI::NormIntInterfaceMPI( const string& normIntFile ) :
NormIntInterface( normIntFile )
{}

NormIntInterfaceMPI::~NormIntInterfaceMPI() {
  
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
    
    // if we are not the master, send generated and accepted events
    // to master and wait for acknowledge
    
    int thisEvents = numGenEvents();
    MPI_Send( &thisEvents, 1, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD );
    
    thisEvents = numAccEvents();
    MPI_Send( &thisEvents, 1, MPI_INT, 0, MPITag::kIntSend, MPI_COMM_WORLD );
    
    MPI_Recv( &thisEvents, 1, MPI_INT, 0, MPITag::kAcknowledge, MPI_COMM_WORLD,
             &status );
  }
}

void
NormIntInterfaceMPI::sumIntegrals( IntType type )
{  
  // first collect integrals from all nodes
  map< string, map< string, complex< double > > >
  integrals = ( type == kNormInt ? getNormIntegrals() : getAmpIntegrals() );
  
  // at this point, the master should be holding an array full of zeroes
  // and other nodes will hold integrals for their subsets of data
  
  MPI_Status status;
  
  // now loop over the pairs of amplitudes
  for( map< string, map< string, complex< double > > >::const_iterator
      ampItr = integrals.begin();
      ampItr != integrals.end();
      ++ampItr ){
    
    string ampName = ampItr->first;
    
    for( map< string, complex< double > >::const_iterator cnjItr =
        ampItr->second.begin();
        cnjItr != ampItr->second.end();
        ++cnjItr ){
      
      string cnjName = cnjItr->first;
      
      // if we are the master node, then need to add up contributions from 
      // other nodes -- otherwise we need to send integrals to the master node
      if( m_isMaster ){
        
        for( int i = 1; i < m_numProc; ++i ){
          
          // cout << "Waiting for " << ampName << "*" << cnjName 
          //     << " from process " << i << endl;
          
          double thisIntegral[2];
          MPI_Recv( thisIntegral, 2, MPI_DOUBLE, i, MPITag::kDoubleSend,
                   MPI_COMM_WORLD, &status );
          integrals[ampName][cnjName] += complex< double >( thisIntegral[0],
                                                            thisIntegral[1] );
        }
        
        // renormalize integrals
        integrals[ampName][cnjName] /= numGenEvents();
        
        // send the renormalized integral to all processes
        for( int i = 1; i < m_numProc; ++i ){
          
          double thisIntegral[] = 
          { numGenEvents() * integrals[ampName][cnjName].real(),
            numGenEvents() * integrals[ampName][cnjName].imag() };
          
          MPI_Send( thisIntegral, 2, MPI_DOUBLE, i, MPITag::kDoubleSend,
                   MPI_COMM_WORLD );
        }
      }
      else{
        
        // scale integral up by generated events
        // master will divide by total number generated in the end
        
        double thisIntegral[] = 
        { numGenEvents() * integrals[ampName][cnjName].real(),
          numGenEvents() * integrals[ampName][cnjName].imag() };
        
        // cout << "Process " << m_rank << " sending "
        //     << ampName << "*" << cnjName << endl;
        
        MPI_Send( thisIntegral, 2, MPI_DOUBLE, 0, MPITag::kDoubleSend,
                 MPI_COMM_WORLD );
        
        MPI_Recv( thisIntegral, 2, MPI_DOUBLE, 0, MPITag::kDoubleSend,
                 MPI_COMM_WORLD, &status );
        
        integrals[ampName][cnjName] = complex< double >( thisIntegral[0],
                                                        thisIntegral[1] );
      }
      
      // regardless of node, set local copy to be the summed version
      if( type == kNormInt ){
        
        setNormIntegral( ampName, cnjName, integrals[ampName][cnjName] );
      }
      else{
        
        setAmpIntegral( ampName, cnjName, integrals[ampName][cnjName] );
      }
    }
  }
  
  // cout << "Parallel NI calculation complete!" << endl;
}
