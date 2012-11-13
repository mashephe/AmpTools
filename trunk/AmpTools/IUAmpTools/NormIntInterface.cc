
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
#include <vector>
#include <string>
#include <complex>
#include <cassert>

#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/DataReader.h"

NormIntInterface::NormIntInterface() :
m_pAmpManager( NULL ),
m_accMCReader( NULL ),
m_genMCReader( NULL ),
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{}

NormIntInterface::NormIntInterface( const string& normIntFile ) :
m_pAmpManager( NULL ),
m_accMCReader( NULL ),
m_genMCReader( NULL ),
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{
  
  cout << "Reading cached normalization integral calculation from: "
  << normIntFile << endl;
  
  ifstream inFile( normIntFile.c_str() );
  
  if( !inFile ){
    cout << "NormIntInterface WARNING:  Could not find file "
    << normIntFile << endl;
    assert( false );
  }
  
  inFile >> (*this);
}


NormIntInterface::NormIntInterface( DataReader* genMCData, 
                                   DataReader* accMCData, 
                                   const AmplitudeManager& ampManager ) :
m_pAmpManager( &ampManager ),
m_accMCReader( accMCData ),
m_genMCReader( genMCData ),
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{
  assert( ( m_accMCReader != NULL ) && ( m_genMCReader != NULL ) );
  
  m_nGenEvents = m_genMCReader->numEvents();
  m_nAccEvents = m_accMCReader->numEvents(); 
  
  // start by initalizing the caches to zero -- MPI implementation will
  // iteratate over and update the cache so setting up the keys in the
  // map is essential
  
  vector< string > amps = ampManager.getAmpNames();
  for( vector< string >::iterator amp = amps.begin();
      amp != amps.end();
      ++amp ) {
    for( vector< string >::iterator conjAmp = amps.begin();
        conjAmp != amps.end();
        ++conjAmp ) {
      
      m_normIntCache[*amp][*conjAmp] = complex< double >( 0, 0 );
      m_ampIntCache[*amp][*conjAmp] = complex< double >( 0, 0 );      
    }
  }
  
  // avoid computation of normalization integrals until they are 
  // actually requested -- this avoids unnecessary expensive computations
}


istream&
NormIntInterface::loadNormIntCache( istream& input )
{
  input >> m_nGenEvents >> m_nAccEvents;
  
  int numAmps;
  input >> numAmps;
  
  vector<string> ampNames;
  
  for( int i = 0; i < numAmps; ++i ){
    
    string name;
    input >> name;
    ampNames.push_back( name );
  }
  
  for( int i = 0; i < numAmps; ++i ){
    for( int j = 0; j < numAmps; ++j ){
      
      complex< double > integral;
      input >> integral;
      m_ampIntCache[ampNames[i]][ampNames[j]] = integral;
    }
  }
  
  for( int i = 0; i < numAmps; ++i ){
    for( int j = 0; j < numAmps; ++j ){
      
      complex< double > integral;
      input >> integral;
      m_normIntCache[ampNames[i]][ampNames[j]] = integral;
    }
  }
  
  m_emptyNormIntCache = false;
  m_emptyAmpIntCache = false;
}

bool
NormIntInterface::hasAccessToMC() const
{
  
  return( ( m_accMCReader != NULL ) &&
         ( m_genMCReader != NULL ) &&
         ( m_pAmpManager != NULL ) );
}

bool
NormIntInterface::hasNormInt( string amp, string conjAmp ) const
{
  map< string, map< string, complex < double > > >::const_iterator 
  conjMap = m_normIntCache.find( amp );
  if( conjMap == m_normIntCache.end() ) return false;
  
  map< string, complex< double > >::const_iterator integral = conjMap->second.find( conjAmp );
  if( integral == conjMap->second.end() ) return false;
  
  return true;
}


complex< double > 
NormIntInterface::normInt( string amp, string conjAmp, bool forceUseCache ) const
{
  
  // pass in true flag to delay computation of amplitude integrals as
  // long as possible -- in the case of a fit with free parameters, this
  // should allow just one computation of the amplitude integrals at the
  // end of the fit
  
  if( m_emptyNormIntCache ) forceCacheUpdate( true );
  
  if( !hasAccessToMC() && ( m_pAmpManager != NULL ) && m_pAmpManager->hasAmpWithFreeParam() ){
    
    cout << "ERROR: the AmplitudeManager has amplitudes that contain free\n"
    << "       parameters, but no MC has been provided to recalculate\n"
    << "       the normalization integrals.  Check that config file\n"
    << "       lists appropriate data and MC sources." << endl;
    
    assert( false );
  }
  
  if( !forceUseCache && hasAccessToMC() && 
     ( m_pAmpManager != NULL ) && m_pAmpManager->hasAmpWithFreeParam() ){
    
    m_normIntCache = 
    m_pAmpManager->calcIntegrals( m_mcVecs, m_nGenEvents, false );
  }
  
  map< string, map< string, complex < double > > >::const_iterator 
  conjMap = m_normIntCache.find( amp );
  
  if( conjMap == m_normIntCache.end() ){
    
    cout << "ERROR: normalization integral does not exist for " << amp << endl;
    assert( false );
  }
  
  map< string, complex< double > >::const_iterator integral = conjMap->second.find( conjAmp );
  if( integral == conjMap->second.end() ){
    
    cout << "ERROR: normalization integral does not exist for " << conjAmp << "*" << endl;
    assert( false );
  }
  
  return integral->second;
}


bool
NormIntInterface::hasAmpInt( string amp, string conjAmp ) const
{
  
  map< string, map< string, complex < double > > >::const_iterator 
  conjMap = m_ampIntCache.find( amp );
  if( conjMap == m_ampIntCache.end() ) return false;
  
  map< string, complex< double > >::const_iterator integral = conjMap->second.find( conjAmp );
  if( integral == conjMap->second.end() ) return false;
  
  return true;
}


complex< double > 
NormIntInterface::ampInt( string amp, string conjAmp, bool forceUseCache ) const
{
  
  if( m_emptyAmpIntCache ) forceCacheUpdate();
  
  if( !forceUseCache && ( m_pAmpManager != NULL ) && m_pAmpManager->hasAmpWithFreeParam() ){
    
    cout << "WARNING:  request of for numerical integral of amplitude\n"
    << "    that has floating parameters using *generated* MC.\n"
    << "    (This is not typically used in a fit -- check for bugs!)\n" 
    << "    Providing original cached value which may be incorrect if\n"
    << "    the parameter has changed since initialization."
    << endl;
    
    // problem: this *is* used in a fit if the renormalizeAmps is turned on
    // and there is a floating parameter in one of the amplitudes
    // add code here to load up generated MC and hang onto it
    // print a notice that we're hogging memory in this configuration
    // for now renormalize is disabled in AmplitudeManager to avoid
    // confusing the user
  }
  
  map< string, map< string, complex < double > > >::const_iterator 
  conjMap = m_ampIntCache.find( amp );
  
  if( conjMap == m_ampIntCache.end() ){
    
    cout << "ERROR: amplitude integral does not exist for " << amp << endl;
    assert( false );
  }
  
  map< string, complex< double > >::const_iterator integral = conjMap->second.find( conjAmp );
  if( integral == conjMap->second.end() ){
    
    cout << "ERROR: amplitude integral does not exist for " << conjAmp << "*" << endl;
    assert( false );
  }
  
  return integral->second;
}

void
NormIntInterface::forceCacheUpdate( bool normIntOnly ) const
{
  
  // if the accepted MC is not available, then the data have likely
  // not been loaded into AmpVecs and we can't reclaculate the integrals
  // below
  assert( m_accMCReader != NULL );
  
  if( !m_emptyNormIntCache && normIntOnly ){
    
    // we can assume that m_mcVecs contains the accepted MC since the
    // generated MC is never left in m_mcVecs
    
    m_normIntCache = m_pAmpManager->calcIntegrals( m_mcVecs, m_nGenEvents, false );
    
    // stop here if we just want the norm ints -- this will likely be the normal
    // mode of operation for NI recalculations during a fit
    return;
  }
  
  // now we need to have the generated MC in addition to the accepted
  // MC in order to be able to continue
  assert( m_genMCReader != NULL );
  
  // flush this if anything is loaded
  m_mcVecs.deallocAmpVecs();
  
  cout << "Loading generated Monte Carlo..." << endl;
  m_mcVecs.loadData( m_genMCReader );
  m_mcVecs.allocateAmps( *m_pAmpManager );
  cout << "\tDone.\n" << flush;
  
  cout << "Calculating integrals..." << endl;
  m_ampIntCache = 
  m_pAmpManager->calcIntegrals( m_mcVecs, m_nGenEvents );
  cout << "\tDone." << endl;
  
  m_emptyAmpIntCache = false;
  
  if( m_accMCReader == m_genMCReader )
  {    
    // optimization for perfect acceptance  
    cout << "Perfect acceptance -- using integrals from generated MC" << endl;
    
    m_normIntCache = m_ampIntCache;
  }
  else
  {  
    // load the accepted MC -- always leave accepted MC in memory
    m_mcVecs.deallocAmpVecs();
    
    cout << "Loading acccepted Monte Carlo..." << endl;
    m_mcVecs.loadData( m_accMCReader );
    m_mcVecs.allocateAmps( *m_pAmpManager );
    cout << "\tDone." << endl;
    
    // since we have flushed the amplitudes from AmpVecs we need
    // to recalcualte them or else subsequent calls to calcIntegrals
    // with firstPass set to false will fail
    cout << "Calculating integrals..." << endl;
    m_normIntCache = 
    m_pAmpManager->calcIntegrals( m_mcVecs, m_nGenEvents );
    
    cout << "\tDone." << endl;
  }
  
  m_emptyNormIntCache = false;
  
  // this method always leaves class in state with m_mcVecs holding accepted MC
}  

void
NormIntInterface::exportNormIntCache( const string& fileName, bool renormalize) const
{
  ofstream out( fileName.c_str() );
  exportNormIntCache( out, renormalize );
}

void
NormIntInterface::exportNormIntCache( ostream& out, bool renormalize) const
{
  if( m_emptyNormIntCache || m_emptyAmpIntCache ) forceCacheUpdate();
  
  out << m_nGenEvents << "\t" << m_nAccEvents << endl;
  
  out << m_normIntCache.size() << endl;
  
  for( map< string, map< string, complex< double > > >::const_iterator 
      ampItr = m_normIntCache.begin();
      ampItr != m_normIntCache.end(); 
      ++ampItr ){
    
    out << ampItr->first << endl;
  }
  
  // write out generated matrix first
  for( map< string, map< string, complex< double > > >::const_iterator 
      ampItr = m_ampIntCache.begin();
      ampItr != m_ampIntCache.end(); 
      ++ampItr ){
    
    for( map< string, complex< double > >::const_iterator 
        conjItr = ampItr->second.begin();
        conjItr != ampItr->second.end();
        ++conjItr ){
      
      complex< double > value = conjItr->second;
      
      if( renormalize ){
        
        value /= sqrt( ampInt( ampItr->first, ampItr->first ) *
                      ampInt( conjItr->first, conjItr->first ) );
      }
      
      out << value << "\t";
    }
    
    out << endl;
  }
  
  // then accepted matrix
  for( map< string, map< string, complex< double > > >::const_iterator 
      ampItr = m_normIntCache.begin();
      ampItr != m_normIntCache.end(); 
      ++ampItr ){
    
    for( map< string, complex< double > >::const_iterator 
        conjItr = ampItr->second.begin();
        conjItr != ampItr->second.end();
        ++conjItr ){
      
      complex< double > value = conjItr->second;
      
      if( renormalize ){
        
        value /= sqrt( ampInt( ampItr->first, ampItr->first ) *
                      ampInt( conjItr->first, conjItr->first ) );
      }
      
      out << value << "\t";
    }
    
    out << endl;
  }
}

map< string, map< string, complex< double > > >
NormIntInterface::getAmpIntegrals() const { 
  
  return m_ampIntCache; 
  
}

map< string, map< string, complex< double > > >
NormIntInterface::getNormIntegrals() const { 
  
  return m_normIntCache; 
}


void 
NormIntInterface::setAmpIntegral( string ampName, string cnjName,
                                 complex< double > val ) const { 
  
  m_ampIntCache[ampName][cnjName] = val; 
  m_emptyAmpIntCache = false;
}

void 
NormIntInterface::setNormIntegral( string ampName, string cnjName,
                                  complex< double > val ) const { 
  
  m_normIntCache[ampName][cnjName] = val; 
  m_emptyNormIntCache = false;
}

