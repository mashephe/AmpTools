
#include <mpi.h>
#include <sstream>

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

#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpToolsMPI/MPITag.h"

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
ParameterManager( ampManager ),
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
ParameterManager( ampManagers ),
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
 
  // we need to deallocate the memory on the worker nodes that
  // was used to hold the parameter values
  
  if( !m_isMaster ){
  
    for( map< string, complex< double >* >::iterator mapItr = m_prodParMap.begin();
        mapItr != m_prodParMap.end();
        ++mapItr )
      delete mapItr->second;
  
    for( map< string, double* >::iterator mapItr = m_ampParMap.begin();
        mapItr != m_ampParMap.end();
        ++mapItr )
      delete mapItr->second;
  }
}

void
ParameterManagerMPI::setupMPI()
{
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );

  m_isMaster = ( m_rank == 0 );
}

void ParameterManagerMPI::addProductionParameter( const string& ampName, 
                                                  bool real )
{
  
  if( m_isMaster ){
    
    // utilize the base class functionality for the master node
    ParameterManager::addProductionParameter( ampName, real );
  
    // then use member data to hang onto a pointer to the complex number
    m_prodParMap[ ampName ] = getProdParPtr( ampName );
  }
  else {

    // find the Amplitude Manager that has this amplitude
    
    vector< AmplitudeManager* >::iterator ampManPtr = m_ampManagers.begin();
    for( ; ampManPtr != m_ampManagers.end(); ++ampManPtr ){
      if( (*ampManPtr)->hasProductionAmp( ampName ) ) break;
    }
    
    if( ampManPtr == m_ampManagers.end() ){

      cout << "ParameterManager ERROR: Could not find production amplitude for " 
           << ampName << endl;
      assert( false );
    }
    
    // now allocate memory to hold the parameter values on the worker nodes
    // initialize with the inital value from the amplitude manager.
    
    m_prodParMap[ampName] = 
      new complex< double >( (**ampManPtr).productionAmp( ampName ) );
    
    // and tell the amplitude manager to look to this memory for updates
    (**ampManPtr).setExternalProductionAmplitude( ampName,
                                                  m_prodParMap[ampName] );
  }
}

void ParameterManagerMPI::addAmplitudeParameter( const string& ampName, 
                                                 const ParameterInfo* parInfo  )
{

  const string& parName = parInfo->parName();

  if( m_isMaster ){
    
    // utilize the base class functionality for the master node
    ParameterManager::addAmplitudeParameter( ampName, parInfo );
  
    // then use member data to hang onto a pointer to the double
    m_ampParMap[parName] = getAmpParPtr( parName );
  }
  else{
    
    // if this is a new parameter, we need to allocate memory for it    
    if( m_ampParMap.find( parName ) == m_ampParMap.end() ){
      
      // need to allocate memory for this parameter
      m_ampParMap[parName] = new double( parInfo->value() );
    }

    // now tell amplitudes to look to this memory for the value
    bool foundOne = false;
    vector< AmplitudeManager* >::iterator ampManPtr = m_ampManagers.begin();
    for( ; ampManPtr != m_ampManagers.end(); ++ampManPtr ){
      
      if( !(*ampManPtr)->hasProductionAmp( ampName ) ) continue;

      foundOne = true;
      
      if( parInfo->fixed() ){
        
        // if it is fixed just go ahead and set the parameter by value
        // this prevents Amplitude class from thinking that is has
        // a free parameter
        
        (**ampManPtr).setAmpParValue( ampName, parName, parInfo->value() );
      }
      else{
        
        (**ampManPtr).setAmpParPtr( ampName, parName, m_ampParMap[parName] );
      }      
    }
    
    if( !foundOne ){
      
      cout << "WARNING:  could not find amplitude named " << ampName 
           << "          while trying to set parameter " << parName << endl;
    }
  }
}


void ParameterManagerMPI::updateParameters()
{
  
  int nameLength;
  char parName[kMaxNameLength];
  double cplxValue[2], value;
  
  MPI_Status status;
  
  if( m_isMaster ){
    
    for( int proc = 1; proc < m_numProc; ++proc ){
      
      for( map< string, complex< double >* >::iterator 
          parItr = m_prodParMap.begin();
          parItr != m_prodParMap.end();
          ++parItr ){
        
        nameLength = parItr->first.length();
        assert( nameLength <= kMaxNameLength );
        strcpy( parName, parItr->first.c_str() );
        
        MPI_Send( &nameLength, 1, MPI_INT, proc,
                 MPITag::kIntSend, MPI_COMM_WORLD );
        MPI_Send( parName, nameLength, MPI_CHAR, proc,
                 MPITag::kCharSend, MPI_COMM_WORLD );
        
        cplxValue[0] = parItr->second->real();
        cplxValue[1] = parItr->second->imag();
        MPI_Send( cplxValue, 2, MPI_DOUBLE, proc,
                 MPITag::kDoubleSend, MPI_COMM_WORLD );
      }

      // send zero namelength to indicate to the workers that the
      // production parameter send is complete
      nameLength = 0;
      MPI_Send( &nameLength, 1, MPI_INT, proc,
               MPITag::kIntSend, MPI_COMM_WORLD );
      
      for( map< string, double* >::iterator parItr = m_ampParMap.begin();
           parItr != m_ampParMap.end();
           ++parItr ){
        
        nameLength = parItr->first.length();
        assert( nameLength <= kMaxNameLength );
        strcpy( parName, parItr->first.c_str() );
        
        MPI_Send( &nameLength, 1, MPI_INT, proc,
                 MPITag::kIntSend, MPI_COMM_WORLD );
        MPI_Send( parName, nameLength, MPI_CHAR, proc,
                 MPITag::kCharSend, MPI_COMM_WORLD );
        
        value = *(parItr->second);
        
        MPI_Send( &value, 1, MPI_DOUBLE, proc,
                 MPITag::kDoubleSend, MPI_COMM_WORLD );
      }
      
      // send zero namelength to indicate to the workers that the
      // amplitude parameter send is complete
      nameLength = 0;
      MPI_Send( &nameLength, 1, MPI_INT, proc,
               MPITag::kIntSend, MPI_COMM_WORLD );
    }
  }
  else{
    
    // the workers collect the production parameters watching for zero namelength
    // which indicates the send is over
    
    MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
             MPI_COMM_WORLD, &status );
    
    while( nameLength != 0 ){
      
      MPI_Recv( parName, nameLength, MPI_CHAR, 0, MPITag::kCharSend,
               MPI_COMM_WORLD, &status );
      MPI_Recv( cplxValue, 2, MPI_DOUBLE, 0, MPITag::kDoubleSend,
               MPI_COMM_WORLD, &status );
      
      ostringstream name;
      for( int i = 0; i < nameLength; ++i ) name << parName[i];
      
      complex< double > parVal( cplxValue[0], cplxValue[1] );
      
      // check that the name is sane and the reset the value
      map< string, complex< double >* >::iterator parItr = 
      m_prodParMap.find( name.str() );
      assert( parItr != m_prodParMap.end() );
      (*(parItr->second)) = parVal;
      
      // fetch the length of the next parameter name which will
      // be zero if we have reached the end of the chain of parameters
      MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
               MPI_COMM_WORLD, &status );
    }
    
    // the workers collect the amplitude parameters watching for zero namelength
    // which indicates the send is over -- note the difference in data size
    // require two such loops
    
    MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
             MPI_COMM_WORLD, &status );
    
    while( nameLength != 0 ){
      
      MPI_Recv( parName, nameLength, MPI_CHAR, 0, MPITag::kCharSend,
               MPI_COMM_WORLD, &status );
      MPI_Recv( &value, 1, MPI_DOUBLE, 0, MPITag::kDoubleSend,
                MPI_COMM_WORLD, &status );
      
      ostringstream name;
      for( int i = 0; i < nameLength; ++i ) name << parName[i];
            
      // check that the name is sane and the reset the value
      map< string, double* >::iterator parItr = m_ampParMap.find( name.str() );
      assert( parItr != m_ampParMap.end() );
      (*(parItr->second)) = value;
      
      // fetch the length of the next parameter name which will
      // be zero if we have reached the end of the chain of parameters
      MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
               MPI_COMM_WORLD, &status );
    }
  }
}


void ParameterManagerMPI::updateAmpParameter( const string& parName )
{
  
  int nameLength;
  char parNameArr[kMaxNameLength];
  double value;
  
  MPI_Status status;
  
  if( m_isMaster ){
    
    for( int proc = 1; proc < m_numProc; ++proc ){
      
      map< string, double* >::iterator parItr =
      m_ampParMap.find( parName );
      
      // we should *always* be able to find the parameter with the 
      // the appropriate name
      assert( parItr != m_ampParMap.end() );
              
      nameLength = parItr->first.length();
      assert( nameLength <= kMaxNameLength );
      strcpy( parNameArr, parItr->first.c_str() );
        
      MPI_Send( &nameLength, 1, MPI_INT, proc,
               MPITag::kIntSend, MPI_COMM_WORLD );
      MPI_Send( parNameArr, nameLength, MPI_CHAR, proc,
               MPITag::kCharSend, MPI_COMM_WORLD );
      
      value = *(parItr->second);
        
      MPI_Send( &value, 1, MPI_DOUBLE, proc,
               MPITag::kDoubleSend, MPI_COMM_WORLD );
    }
  }
  else{
    
    MPI_Recv( &nameLength, 1, MPI_INT, 0, MPITag::kIntSend,
             MPI_COMM_WORLD, &status );
    MPI_Recv( parNameArr, nameLength, MPI_CHAR, 0, MPITag::kCharSend,
             MPI_COMM_WORLD, &status );
    
    MPI_Recv( &value, 1, MPI_DOUBLE, 0, MPITag::kDoubleSend,
             MPI_COMM_WORLD, &status );
      
    ostringstream name;
    for( int i = 0; i < nameLength; ++i ) name << parNameArr[i];
      
    // check that the name is sane and the reset the value
    map< string, double* >::iterator parItr = m_ampParMap.find( name.str() );
    assert( parItr != m_ampParMap.end() );
    (*(parItr->second)) = value;
    
    // now we need to trigger the update method using the base class
    // this will call the update routines in the individual amplitudes
        
    ParameterManager::update( name.str() );
  }
}

void
ParameterManagerMPI::update( const string& parName ){
  
  // should only be called directly on the master via
  // the call to virtual ParameterManager::update 
  
  assert( m_isMaster );
 
  // this puts workers into the updateAmpParameter routine above 
  LikelihoodManagerMPI::broadcastToFirst( LikelihoodManagerMPI::kUpdateAmpParameter );
  
  // now put the master there
  updateAmpParameter( parName );
}

