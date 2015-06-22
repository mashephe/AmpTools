
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
ParameterManagerMPI( MinuitMinimizationManager* minuitManager,
                    IntensityManager* intenManager ) :
ParameterManager( minuitManager, intenManager )
{
  setupMPI();
  
  if( !m_isMaster ){
    
    cerr << "Instance of MinuitMinimizationManager exists on worker node"
    << endl;
    
    assert( false );
  }
}

ParameterManagerMPI::
ParameterManagerMPI( MinuitMinimizationManager* minuitManager,
                    const vector<IntensityManager*>& intenManagers ) :
ParameterManager( minuitManager, intenManagers )
{
  setupMPI();
  
  if( !m_isMaster ){
    
    cerr << "Instance of MinuitMinimizationManager exists on worker node"
    << endl;
    
    assert( false );
  }
}

ParameterManagerMPI::
ParameterManagerMPI( IntensityManager* intenManager ) :
ParameterManager( intenManager ),
m_intenManagers( 0 )
{
  setupMPI();
  
  if( m_isMaster ){
    
    cerr << "Master ParameterManager has no MinuitMinimizationManager"
    << endl;
    
    assert( false );
  }
  
  m_intenManagers.push_back( intenManager );
}

ParameterManagerMPI::
ParameterManagerMPI( const vector< IntensityManager* >& intenManagers ) :
ParameterManager( intenManagers ),
m_intenManagers( intenManagers )
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

void ParameterManagerMPI::addProductionParameter( const string& termName,
                                                  bool real, bool fixed )
{
  
  if( m_isMaster ){
    
    // utilize the base class functionality for the master node
    ParameterManager::addProductionParameter( termName, real, fixed );
  
    // then use member data to hang onto a pointer to the complex number
    m_prodParMap[ termName ] = getProdParPtr( termName );
  }
  else {

    // find the Amplitude Manager that has this amplitude
    
    vector< IntensityManager* >::iterator intenManPtr = m_intenManagers.begin();
    for( ; intenManPtr != m_intenManagers.end(); ++intenManPtr ){
      if( (*intenManPtr)->hasTerm( termName ) ) break;
    }
    
    if( intenManPtr == m_intenManagers.end() ){

      cout << "ParameterManager ERROR: Could not find production amplitude for " 
           << termName << endl;
      assert( false );
    }
    
    // now allocate memory to hold the parameter values on the worker nodes
    // initialize with the inital value from the amplitude manager.
    
    m_prodParMap[termName] = 
      new complex< double >( (**intenManPtr).productionFactor( termName ) );
    
    // and tell the amplitude manager to look to this memory for updates
    (**intenManPtr).setExternalProductionFactor( termName,
						 m_prodParMap[termName] );
  }
}

void ParameterManagerMPI::addAmplitudeParameter( const string& termName, 
                                                 const ParameterInfo* parInfo  )
{

  const string& parName = parInfo->parName();

  if( m_isMaster ){
    
    // utilize the base class functionality for the master node
    ParameterManager::addAmplitudeParameter( termName, parInfo );
  
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
    vector< IntensityManager* >::iterator intenManPtr = m_intenManagers.begin();
    for( ; intenManPtr != m_intenManagers.end(); ++intenManPtr ){
      
      if( !(*intenManPtr)->hasTerm( termName ) ) continue;

      foundOne = true;
      
      if( parInfo->fixed() ){
        
        // if it is fixed just go ahead and set the parameter by value
        // this prevents Amplitude class from thinking that is has
        // a free parameter
        
        (**intenManPtr).setParValue( termName, parName, parInfo->value() );
      }
      else{
        
        (**intenManPtr).setParPtr( termName, parName, m_ampParMap[parName] );
      }      
    }
    
    if( !foundOne ){
      
      cout << "WARNING:  could not find term named " << termName 
           << "          while trying to set parameter " << parName << endl;
    }
  }
}


void ParameterManagerMPI::updateParameters()
{
  // pointers to parameters are stored in maps on both the
  // master and worker nodes -- these maps contain the
  // same keys and should be sorted in the same way
  //
  // iterate over the map and pack an array and then broadcast
  // this from the master to the workers, then copy back
  // out of the array into the map on the workers
  
  double* parData = new double[2*m_prodParMap.size()+m_ampParMap.size()];
  
  int i = 0;
  for( map< string, complex< double >* >::iterator
      parItr = m_prodParMap.begin();
      parItr != m_prodParMap.end();
      ++parItr ) {
    
    parData[i++] = real( *(parItr->second) );
    parData[i++] = imag( *(parItr->second) );
  }
  
  for( map< string, double* >::iterator
      parItr = m_ampParMap.begin();
      parItr != m_ampParMap.end();
      ++parItr ) {

    parData[i++] = *(parItr->second);
  }

  MPI_Bcast( parData, i, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  
  if( !m_isMaster ){
    
    i = 0;
    for( map< string, complex< double >* >::iterator
        parItr = m_prodParMap.begin();
        parItr != m_prodParMap.end();
        ++parItr ) {
      
      (*(parItr->second)) = complex< double >( parData[i],
                                               parData[i+1] );
      i += 2;
    }
    
    for( map< string, double* >::iterator
        parItr = m_ampParMap.begin();
        parItr != m_ampParMap.end();
        ++parItr ) {
      
      (*(parItr->second)) = parData[i++];
    }
  }
  
  delete[] parData;
}


void ParameterManagerMPI::updateAmpParameter( const string& parName )
{
  
  int nameLength;
  char parNameArr[kMaxNameLength];
  double value;
  
  MPI_Status status;
  
  if( m_isMaster ){
    
    map< string, double* >::iterator parItr =
    m_ampParMap.find( parName );
    
    // we should *always* be able to find the parameter with the
    // the appropriate name
    assert( parItr != m_ampParMap.end() );
    
    nameLength = parItr->first.length();
    assert( nameLength <= kMaxNameLength );
    strcpy( parNameArr, parItr->first.c_str() );
    
    MPI_Bcast( &nameLength, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( parNameArr, nameLength, MPI_CHAR, 0, MPI_COMM_WORLD );
    
    value = *(parItr->second);
    
    MPI_Bcast( &value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
  }
  else{
    
    MPI_Bcast( &nameLength, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast( parNameArr, nameLength, MPI_CHAR, 0, MPI_COMM_WORLD );
    MPI_Bcast( &value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
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

