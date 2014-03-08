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
#include <list>
#include <complex>
#include <string.h>
#include <algorithm>
#include <cassert>

#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/Kinematics.h"

IntensityManager::IntensityManager( const vector< string >& reaction,
                                   const string& reactionName) :
m_reactionName(reactionName),
m_renormalizeTerms( false ),
m_normInt( NULL )
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

