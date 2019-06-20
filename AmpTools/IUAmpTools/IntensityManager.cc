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
#include <sstream>
#include <fstream>

#include <sys/time.h>

#include <vector>
#include <string>
#include <map>
#include <set>
#include <list>
#include <complex>
#include <string.h>
#include <algorithm>
#include <cassert>

#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/Term.h"

IntensityManager::IntensityManager( const vector< string >& reaction,
                                   const string& reactionName) :
m_reactionName(reactionName),
m_renormalizeTerms( false ),
m_normInt( NULL ),
m_needsUserVarsOnly( true ),
m_optimizeParIteration( false ),
m_flushFourVecsIfPossible( false ),
m_forceUserVarRecalculation( false )
{

}

double
IntensityManager::calcIntensity( const Kinematics* kinematics ) const
{
  
  // Create and fill an AmpVecs object with a single event
  AmpVecs aVecs;
  aVecs.loadEvent(kinematics);
  aVecs.allocateTerms(*this,true);
  
  // Calculate the intensity based on this one event
  calcIntensities(aVecs);
  GDouble intensity = aVecs.m_pdIntensity[0];
  
  // Deallocate memory and return
  aVecs.deallocAmpVecs();
  
  return intensity;
  
}

int
IntensityManager::addTerm( const string& name,
                           const string& scale ){
  
  // track the name and sum based on index
  m_termNames.push_back( name );
  
  // record the index
  m_termIndex[name] = m_termNames.size() - 1;
    
  m_prodFactorVec.push_back( static_cast< complex< double >* >( 0 ) );
  m_termScaleVec.push_back( AmpParameter( scale ) );
  
  setDefaultProductionFactor( name, complex< double >( 1, 0 ) );

  return m_termIndex[name];
}


const vector< string >&
IntensityManager::getTermNames() const {
  
  return m_termNames;
}

const AmpParameter&
IntensityManager::getScale( const string& name ) const {
  
  return m_termScaleVec[termIndex(name)];
}

int
IntensityManager::termIndex( const string& termName ) const {
  
  map< string, int >::const_iterator mapItr = m_termIndex.find( termName );
  assert( mapItr != m_termIndex.end() );
  
  return mapItr->second;
}

bool
IntensityManager::hasTerm(const string& name) const {
  
  map< string, const complex< double >* >::const_iterator prodItr = m_prodFactor.find( name );
  return (prodItr != m_prodFactor.end()) ? true : false;
}


complex< double >
IntensityManager::productionFactor( const string& name ) const {
  
  if( m_prodFactor.find( name ) == m_prodFactor.end() ){
    
    cout << "ERROR: cannot provide production amplitude for " << name << endl;
    assert( false );
  }
  
  return productionFactor( m_termIndex.at( name ) );
}

complex< double >
IntensityManager::productionFactor( int ampIndex ) const {
  
  return *m_prodFactorVec.at( ampIndex ) *
  static_cast< double >( m_termScaleVec.at( ampIndex ) );
}

void
IntensityManager::prodFactorArray( double* array ) const {
  
  complex< double > value;
  
  for( int i = 0; i < m_prodFactorVec.size(); ++i ){
    
    value = *m_prodFactorVec[i] *
      static_cast< double >( m_termScaleVec[i] );
  
    array[2*i] = real( value );
    array[2*i+1] = imag( value );
  }
}

unsigned int
IntensityManager::userVarsPerEvent() const {
  
  set< string > countedStaticTerms;
  set< string > countedUniqueTerms;
  
  vector< string > termNames = getTermNames();
  
  unsigned int userStorage = 0;
  
  for( int i = 0; i < getTermNames().size(); i++ ) {
    
    unsigned int iNPermutations = getPermutations( termNames[i] ).size();
    const vector< const Term* >& factorVec = getFactors( termNames[i] );
    
    for( int j = 0; j < factorVec.size(); ++j ){
      
      if( factorVec[j]->areUserVarsStatic() ){
        
        // for factors that have static data, we only want to
        // count the allocation once -- if we have counted
        // it already then skip to the next factor
        if( countedStaticTerms.find( factorVec[j]->name() ) ==
           countedStaticTerms.end() )
          countedStaticTerms.insert( factorVec[j]->name() );
        else continue;
      }
      else{
        // user data is not static, so
        // check to see if we have seen an instance of this
        // same amplitude that would behave in the same way
        // (i.e., has the same arguments) and if it exists, we
        // will use user data block from it instead
        if( countedUniqueTerms.find( factorVec[j]->identifier() ) ==
           countedUniqueTerms.end() )
          countedUniqueTerms.insert( factorVec[j]->identifier() );
        else continue;
      }
      
      userStorage += iNPermutations * factorVec[j]->numUserVars();
    }
  }
  
  return userStorage;
}

void
IntensityManager::calcUserVars( AmpVecs& a ) const
{
  
#ifdef VTRACE
  VT_TRACER( "IntensityManager::calcUserVars" );
#endif
 
  const vector< string >& termNames = getTermNames();
  
  int iNTerms = termNames.size();
  
  assert( iNTerms && a.m_iNEvents && a.m_iNTrueEvents && a.m_pdUserVars);
  
  int iTermIndex;
  unsigned long long iUserVarsOffset = 0;
  for( iTermIndex = 0; iTermIndex < iNTerms; iTermIndex++ )
  {
    const vector< vector< int > >& vvPermuations = getPermutations(termNames[iTermIndex]);
    int iNPerms = vvPermuations.size();
    
    const vector< const Term* >& vTerms = getFactors( termNames[iTermIndex] );
    int iFactor, iNFactors = vTerms.size();
    
    const Term* pCurrTerm = 0;
    for( iFactor=0; iFactor < iNFactors; iFactor++ ){
      
      pCurrTerm = vTerms.at( iFactor );
      
      if( pCurrTerm->numUserVars() == 0 ) continue;
      
      // this is the number of variables for the data set
      int iNVars = pCurrTerm->numUserVars();
      int iNData = iNVars * a.m_iNEvents * iNPerms;
      
      // we will set this based on the algorithm below
      unsigned long long thisOffset = 0;
      
      if( pCurrTerm->areUserVarsStatic() ){
        
        // the user variables are static, so let's look at the
        // list of data pointers stored in the relevant AmpVecs
        // object and see if there is one associated with this
        // amplitude name
        
        map< string, unsigned long long >::const_iterator offsetItr =
        a.m_userVarsOffset.find( pCurrTerm->name() );
        
        if( offsetItr == a.m_userVarsOffset.end() ){
          
          // set the offset to where the calculation will end up
          a.m_userVarsOffset[pCurrTerm->name()] = iUserVarsOffset;
          
          // record it to use in the lines below
          thisOffset = iUserVarsOffset;
          
          // increment
          iUserVarsOffset += iNData;
        }
        else // we have done the calculation already
          if( m_forceUserVarRecalculation ){ // but should redo it
            
            thisOffset = offsetItr->second;
          }
          else{ // and don't want to repeat it
            
            continue;
          }
      }
      else{
        // the variables are not static, repeat the algorithm
        // above but search based on identifier of the amplitude
        
        map< string, unsigned long long >::const_iterator offsetItr =
        a.m_userVarsOffset.find( pCurrTerm->identifier() );
        
        if( offsetItr == a.m_userVarsOffset.end() ){
          
          // set the offset to where the calculation will end up
          a.m_userVarsOffset[pCurrTerm->identifier()] = iUserVarsOffset;
          
          // record it to use in the lines below
          thisOffset = iUserVarsOffset;
          
          // increment -- can only happen at most once either above or here
          iUserVarsOffset += iNData;
        }
        else // we have done the calculation already
          if( m_forceUserVarRecalculation ){ // but should redo it
            
            thisOffset = offsetItr->second;
          }
          else{ // and don't want to repeat it
            
            continue;
          }
      }
      
      // calculation of user-defined kinematics data
      // is something that should only be done once
      // per fit, so do it on the CPU no matter what
      
      pCurrTerm->
      calcUserVarsAll( a.m_pdData,
                      a.m_pdUserVars + thisOffset,
                      a.m_iNEvents, &vvPermuations );
      
#ifdef GPU_ACCELERATION
      
      // we want to reorder the userVars if we are working on the
      // GPU so that the same variable for neighboring events is
      // next to each other, this will enhance block read and
      // caching ability
      
      GDouble* tmpVarStorage = new GDouble[iNData];
      
      for( int iPerm = 0; iPerm < iNPerms; ++iPerm ){
        for( int iEvt = 0; iEvt < a.m_iNEvents; ++iEvt ){
          for( int iVar = 0; iVar < iNVars; ++iVar ){
            
            unsigned long long cpuIndex =
            thisOffset + iPerm*a.m_iNEvents*iNVars + iEvt*iNVars + iVar;
            unsigned long long gpuIndex =
            iPerm*a.m_iNEvents*iNVars + iVar*a.m_iNEvents + iEvt;
            
            tmpVarStorage[gpuIndex] = a.m_pdUserVars[cpuIndex];
          }
        }
      }
      
      memcpy( a.m_pdUserVars + thisOffset, tmpVarStorage,
             iNData*sizeof(GDouble) );
      
      delete[] tmpVarStorage;
#endif //GPU_ACCELERATION
      
    }
  }
  
  // and if we are doing a GPU accelerated fit
  // copy the entire user data block to the GPU so we have it there
  
#ifdef GPU_ACCELERATION
  a.m_gpuMan.copyUserVarsToGPU( a );
#endif //GPU_ACCELERATION
  
  if( needsUserVarsOnly() && m_flushFourVecsIfPossible ) a.clearFourVecs();
  
  return;
}


void
IntensityManager::setDefaultProductionFactor( const string& name,
                                              complex< double > prodFactor )
{
  m_defaultProdFactor[name] = prodFactor;
  m_prodFactor[name] = &(m_defaultProdFactor[name]);
  m_prodFactorVec[m_termIndex[name]] = &(m_defaultProdFactor[name]);
}

void
IntensityManager::setExternalProductionFactor( const string& name,
                                               const complex< double >* prodAmpPtr )
{
  map< string, const complex< double >* >::iterator prodItr = m_prodFactor.find( name );
  
  if( prodItr == m_prodFactor.end() ){
    
    cout << "ERROR:  amplitude " << name << " has no factors!" << endl;
    assert( false );
  }
  
  m_prodFactor[name] = prodAmpPtr;
  m_prodFactorVec[m_termIndex[name]] = prodAmpPtr;
}

void
IntensityManager::resetProductionFactors()
{
  for( map< string, complex< double > >::iterator
      prodItr = m_defaultProdFactor.begin();
      prodItr != m_defaultProdFactor.end();  ++prodItr ){
    
    m_prodFactor[prodItr->first] = &(m_defaultProdFactor[prodItr->first]);
    m_prodFactorVec[m_termIndex[prodItr->first]] = &(m_defaultProdFactor[prodItr->first]);
  }
}

void
IntensityManager::setParPtr( const string& name, const string& parName,
                             const double* ampParPtr ){
   
  if( m_termScaleVec[m_termIndex[name]].name() == parName ){
    
    m_termScaleVec[m_termIndex[name]].setExternalValue( ampParPtr );
  }
}

void
IntensityManager::setParValue( const string& name, const string& parName,
                                 double val ){
  
  if( m_termScaleVec[m_termIndex[name]].name() == parName ){
    
    m_termScaleVec[m_termIndex[name]].setValue( val );
  }
}

/*
 * disable these for now until a fix is made to the NormIntInterface to prevent
 * problems with free parameters in amplitudes
 
 void
 IntensityManager::renormalizeAmps( const NormIntInterface* normInt ){
 
 m_normInt = normInt;
 m_renormalizeAmps = true;
 }
 
 void
 IntensityManager::disableRenormalization(){
 
 m_renormalizeAmps = false;
 }
 */

