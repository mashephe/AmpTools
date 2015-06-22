
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
#include <fstream>
#include <cassert>
#include <string>

#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/ComplexParameter.h"
#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MinuitParameterManager.h"
#include "GPUManager/GPUCustomTypes.h"

ParameterManager::ParameterManager( MinuitMinimizationManager* minuitManager,
                                    IntensityManager* intenManager ) :
MIObserver(),
m_minuitManager( minuitManager ),
m_intenManagers( 0 )
{ 
  m_minuitManager->attach( this );
  m_intenManagers.push_back(intenManager);
//  cout << "Parameter manager initialized." << endl;
}

ParameterManager::
ParameterManager( MinuitMinimizationManager* minuitManager,
                 const vector<IntensityManager*>& intenManagers ) :
MIObserver(),
m_minuitManager( minuitManager ),
m_intenManagers( intenManagers )
{ 
  m_minuitManager->attach( this );
//  cout << "Parameter manager initialized." << endl;
}

// protected constructors: these are used in MPI implementations where
// the ParameterManager is created on a worker node that does not have
// a MinuitMinimizationManager.  The reference must be initialized, so
// it is initialized to NULL but never used.  The class using these
// constructors must override any method that uses the
// MinuitMinimizationManager in order to avoid dereferencing a null
// pointer.

ParameterManager::ParameterManager( IntensityManager* ampManager ) :
  m_minuitManager( NULL ),
  m_intenManagers( 0 )
{ 
  m_intenManagers.push_back(ampManager);
//  cout << "Parameter manager initialized." << endl;
}

ParameterManager::
ParameterManager( const vector<IntensityManager*>& intenManagers ) :
  m_minuitManager( NULL ),
  m_intenManagers( intenManagers )
{ 
//  cout << "Parameter manager initialized." << endl;
}

ParameterManager::~ParameterManager()
{
  for( vector< ComplexParameter* >::iterator parItr = m_prodPtrCache.begin();
      parItr != m_prodPtrCache.end(); 
      ++parItr ){
    
    delete *parItr;
  }
  
  for( vector< MinuitParameter* >::iterator parItr = m_ampPtrCache.begin();
      parItr != m_ampPtrCache.end(); 
      ++parItr ){
    
    delete *parItr;
  }  
  
  for( vector< GaussianBound* >::iterator boundItr = m_boundPtrCache.begin();
      boundItr != m_boundPtrCache.end();
      ++boundItr ){
    
    delete *boundItr;
  }
}

void
ParameterManager::setupFromConfigurationInfo( ConfigurationInfo* cfgInfo ){
  
  vector< AmplitudeInfo* > amps = cfgInfo->amplitudeList();
  
  m_constraintMap = cfgInfo->constraintMap();
  
  // separate creation of amplitude parameters from production parameters
  // in order to make the error matrix more easily separable
  // practically, this just means two loops below
  
  for( vector< AmplitudeInfo* >::iterator ampItr = amps.begin();
      ampItr != amps.end();
      ++ampItr ){
    
    addProductionParameter( (**ampItr).fullName(), (**ampItr).real(), (**ampItr).fixed() );
  }
  
  for( vector< AmplitudeInfo* >::iterator ampItr = amps.begin();
      ampItr != amps.end();
      ++ampItr ){
    
    vector< ParameterInfo* > pars = (**ampItr).parameters();
    
    for( vector< ParameterInfo* >::const_iterator parItr = pars.begin();
        parItr != pars.end();
        ++parItr ){
      
      addAmplitudeParameter( (**ampItr).fullName(), *parItr );
    }
  }
}


void
ParameterManager::addAmplitudeParameter( const string& termName, const ParameterInfo* parInfo ){
  
  const string& parName = parInfo->parName();
  
  // see if this is a parameter that we already know about
  
  map< string, MinuitParameter* >::iterator mapItr = m_ampParams.find( parName );
  MinuitParameter* parPtr;
  
  if( mapItr == m_ampParams.end() ){
    
//    cout << "Creating new amplitude parameter:  " << parInfo->parName() << endl;
    
    parPtr = new MinuitParameter( parName, m_minuitManager->parameterManager(),
                                 parInfo->value());
    
    // attach to allow the parameter to call back this class when it is updated
    parPtr->attach( this );
    
    if( parInfo->fixed() ){
      
      parPtr->fix();
    }
    
    if( parInfo->bounded() ){
      
      parPtr->bound( parInfo->lowerBound(), parInfo->upperBound() );
    }
    
    if( parInfo->gaussianBounded() ){
      
      GaussianBound* boundPtr = 
      new GaussianBound( m_minuitManager, parPtr, parInfo->centralValue(),
                        parInfo->gaussianError() );
      
      m_boundPtrCache.push_back( boundPtr );
    }
    
    // keep track of new objects that are being allocated
    m_ampPtrCache.push_back( parPtr );
    m_ampParams[parName] = parPtr;
  }
  else{
    
    parPtr = mapItr->second;
  }
  
  // find the Amplitude Manager that has the relevant amplitude
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
      
      (**intenManPtr).setParPtr( termName, parName, parPtr->constValuePtr() );
    }
  }
  
  if( !foundOne ){
    
    cout << "WARNING:  could not find amplitude named " << termName 
         << " while trying to set parameter " << parName << endl;
  }
}

void 
ParameterManager::addProductionParameter( const string& termName, bool real, bool fixed )
{
  
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
  
  // get the parameter's initial value from the amplitude manager
  // (The productionAmp method will return the scaled production amplitude
  // we need to divide out the scale to get the production parameter initial
  // value that was specified in the configuration file.)
  complex< double > initialValue = (**intenManPtr).productionFactor( termName ) /
        (double)(**intenManPtr).getScale( termName );
  
  // find the ComplexParameter for this amplitude or an amplitude constrained to
  //   be the same as this amplitude
  
  ComplexParameter* par = findParameter(termName);
  
  // create ComplexParameter from scratch if it doesn't already exist
  
  if (!par){
//    cout << "ParameterManager:  Creating new complex production amplitude parameter for "
//         << termName << endl;
    par = new ComplexParameter( termName, *m_minuitManager, initialValue, real );
    m_prodPtrCache.push_back( par );
  }
  
  if( fixed ) par->fix();
  
  // update the amplitude manager
  
  (**intenManPtr).setExternalProductionFactor( termName,
                                               par->constValuePtr() );
  
  // record this parameter
  
  m_prodParams[termName] = par;
  
}

complex< double >* 
ParameterManager::getProdParPtr( const string& termName ){
  
  map< string, ComplexParameter* >::iterator mapItr 
  = m_prodParams.find( termName );
  
  // make sure we found one
  assert( mapItr != m_prodParams.end() );
  
  return mapItr->second->valuePtr();
}

double* 
ParameterManager::getAmpParPtr( const string& parName ){
  
  map< string, MinuitParameter* >::iterator mapItr = m_ampParams.find( parName );
  
  // make sure we found one
  assert( mapItr != m_ampParams.end() );
  
  return mapItr->second->valuePtr();
}

bool
ParameterManager::hasConstraints( const string& termName ) const{
  map<string, vector<string> >::const_iterator
  mapItr = m_constraintMap.find(termName);
  return (mapItr != m_constraintMap.end()) ? true : false;
}

bool
ParameterManager::hasParameter( const string& termName ) const{
  map<string, ComplexParameter* >::const_iterator
  mapItr = m_prodParams.find(termName);
  return (mapItr != m_prodParams.end()) ? true : false;
}

ComplexParameter*
ParameterManager::findParameter( const string& termName) const{
  
  // return the parameter associated with this amplitude if it is already defined
  
  map<string, ComplexParameter*>::const_iterator pItr = m_prodParams.find(termName);
  if (pItr != m_prodParams.end()) return pItr->second;
  
  // otherwise look for a parameter associated with an amplitude that is
  //   constrained to be the same as this amplitude
  
  map<string, vector<string> >::const_iterator cItr = m_constraintMap.find(termName);
  if (cItr == m_constraintMap.end()) return NULL;
  
  vector<string> constraints = cItr->second; 
  for (unsigned int i = 0; i < constraints.size(); i++){
    pItr = m_prodParams.find(constraints[i]);
    if (pItr != m_prodParams.end()) return pItr->second;
  }
  
  return NULL;
}

void
ParameterManager::update( const MISubject* parPtr ){
  
  // this method is called whenever any parameter changes
  // if it is an amplitude parameter, we want to notify the
  // amplitude of the change
  
  // first loop over the map containing the parameter pointers and
  // try to find the one that matches parPtr
  
  for( map< string, MinuitParameter* >::const_iterator mapItr = m_ampParams.begin();
      mapItr != m_ampParams.end();
      ++mapItr ){
    
    if( mapItr->second == parPtr ){
      
      // we found the relevant param -- now notify all amplitude managers that
      // the parameter has changed
      
      update( mapItr->first );
    }
  }
  
  updateParCovariance();
}

void
ParameterManager::update( const string& parName ){
  
  // useful to have this method available to update by name
  
  for( vector< IntensityManager* >::const_iterator intenMan = m_intenManagers.begin();
      intenMan != m_intenManagers.end();
      ++intenMan ){
    
    (**intenMan).updatePar( parName );
  }  
}  

void
ParameterManager::updateParCovariance(){
  
  // build a vector that provides the MINUIT parameter index i for
  // the real parts of the production parameters, if there is an imaginary
  // part it will have an index of i + 1
  
  vector< int > prodParMinuitIndex( 0 );
  int numMinuitPars = 0;
  for( vector< ComplexParameter* >::const_iterator par = m_prodPtrCache.begin();
      par != m_prodPtrCache.end();
      ++par ){
    
    if( (**par).isFixed() ){

      // if the production parameter is fixed, then it won't
      // have a MINUIT index and were going to set this to kFixedIndex
      // later in building minuitParIndex, but we need an entry in this
      // vector for that ComplexParameter
      
      prodParMinuitIndex.push_back( kFixedIndex );
    }
    else{
  
      prodParMinuitIndex.push_back( numMinuitPars );
      numMinuitPars += ( (**par).isPurelyReal() ? 1 : 2 );
    }
  }
  
  vector< int > minuitParIndex( 0 );
  
  // build the list of parameters by looping over all amplitude managers and looping
  // over all amplitudes
  
  m_parList.clear();
  m_parValues.clear();
  m_parIndex.clear();
  
  int index = 0;
  for( vector< IntensityManager* >::const_iterator intenMan = m_intenManagers.begin();
      intenMan != m_intenManagers.end();
      ++intenMan ){
    
    const vector< string >& termNames = (**intenMan).getTermNames();
    for( vector< string >::const_iterator name = termNames.begin();
        name != termNames.end();
        ++name ){
      
      // this will return the complex parameter associated with
      // this amplitude or the complex parameter to which this
      // amplitude is constrained
      const ComplexParameter* prodPar = findParameter( *name );
      
      // now determine the index of this parameter in the cache
      vector< ComplexParameter* >::const_iterator parItr = 
      find( m_prodPtrCache.begin(), m_prodPtrCache.end(), prodPar );
      assert( parItr != m_prodPtrCache.end() );
      int cacheIndex = parItr - m_prodPtrCache.begin();
      
      // if the parameter is fixed, the real part won't
      // have a MINUIT index, save kFixedIndex (negative) and watch
      // out for this when building the error matrix
      // (this is a redundant check since the check has already been
      //  done above in building prodParMinuitIndex, but it helps
      //  to make the code a little clearer)
      minuitParIndex.push_back( prodPar->isFixed() ?
                                kFixedIndex : prodParMinuitIndex[cacheIndex] );
      
      // record the other values
      m_parList.push_back( (*name) + "_re" );
      m_parValues.push_back( real( prodPar->value() ) );
      m_parIndex[m_parList.back()] = index++;
      
      // if the parameter is purely real or fixed, the imaginary part won't
      // have a MINUIT index, save kFixedIndex (negative) and watch 
      // out for this when building the error matrix
      minuitParIndex.push_back( prodPar->isPurelyReal() || prodPar->isFixed() ?
                                kFixedIndex : prodParMinuitIndex[cacheIndex] + 1 );
      
      // record the imaginary parameter info
      m_parList.push_back( (*name) + "_im" );
      m_parValues.push_back( imag( prodPar->value() ) );
      m_parIndex[m_parList.back()] = index++;
    }
  }
  
  // finish this list by looping over all the amplitude parameters
  
  for( vector< MinuitParameter* >::const_iterator ampPar = m_ampPtrCache.begin();
      ampPar != m_ampPtrCache.end();
      ++ampPar ){
    
    // if the parameter is not floating, it won't have a row in the covariance matrix
    
    if( (**ampPar).floating() ){
      
      // the MINUIT parameter indices number 
      // sequentially after the production parameters
      
      minuitParIndex.push_back( numMinuitPars );
      ++numMinuitPars;
    }
    else{
      
      minuitParIndex.push_back( kFixedIndex );
    }
    
    m_parList.push_back( (**ampPar).name() );
    m_parValues.push_back( (**ampPar).value() );
    m_parIndex[m_parList.back()] = index++;
  }
  
  // now we have a flat list of all the (real valued) parameters
  // we need to fill out the covariance matrix
  
  int nPar = m_parList.size();
  
  const vector< vector< double > >& minCovMtx = 
  m_minuitManager->parameterManager().covarianceMatrix();
  
  m_covMatrix.clear();
  for( int i = 0; i < nPar; ++i ){
    
    // create a new covariance matrix row
    m_covMatrix.push_back( vector< double >( nPar ) );
    
    for( int j = 0; j < nPar; ++j ){
      
      int iIndex = minuitParIndex[i];
      int jIndex = minuitParIndex[j];
      
      // if either i or j are a fixed parameter we set the
      // index 
      
      if( iIndex == kFixedIndex || jIndex == kFixedIndex ){
        
        m_covMatrix[i][j] = 0;
        continue;
      }
      
      // this shouldn't happen -- if it does, something has
      // gone wrong with the parameter numbering or the 
      // construction of the covariance matrix and we should
      // crash
      
      assert( iIndex < minCovMtx.size() );
      assert( jIndex < minCovMtx.size() );
      
      m_covMatrix[i][j] = minCovMtx[iIndex][jIndex];
    }
  }
}
