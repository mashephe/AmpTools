
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
#include <cstring>
#include <iomanip>

#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/DataReader.h"

NormIntInterface::NormIntInterface() :
m_pIntenManager( NULL ),
m_accMCReader( NULL ),
m_genMCReader( NULL ),
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{}

NormIntInterface::NormIntInterface( const string& normIntFile ) :
m_pIntenManager( NULL ),
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
                                   const IntensityManager& intenManager ) :
m_pIntenManager( &intenManager ),
m_accMCReader( accMCData ),
m_genMCReader( genMCData ),
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true ),
m_termNames( intenManager.getTermNames() )
{
  assert( ( m_accMCReader != NULL ) && ( m_genMCReader != NULL ) );
  
  m_termNames = intenManager.getTermNames();
  
  m_nGenEvents = m_genMCReader->numEvents();
  m_nAccEvents = m_accMCReader->numEvents(); 
  
  initializeCache();
}


istream&
NormIntInterface::loadNormIntCache( istream& input )
{
  input >> m_nGenEvents >> m_nAccEvents;
  
  int numTerms;
  input >> numTerms;
  
  m_termNames.clear();
  m_termIndex.clear();
  
  for( int i = 0; i < numTerms; ++i ){
    
    string name;
    input >> name;
    m_termNames.push_back( name );
    m_termIndex[name] = i;
  }
  
  initializeCache();
  
  for( int i = 0; i < numTerms; ++i ){
    for( int j = 0; j < numTerms; ++j ){

      complex< double > integral;
      input >> integral;
      
      m_ampIntCache[2*i*numTerms+2*j]   = real( integral );
      m_ampIntCache[2*i*numTerms+2*j+1] = imag( integral );
    }
  }
  
  for( int i = 0; i < numTerms; ++i ){
    for( int j = 0; j < numTerms; ++j ){
      
      complex< double > integral;
      input >> integral;
      
      m_normIntCache[2*i*numTerms+2*j]   = real( integral );
      m_normIntCache[2*i*numTerms+2*j+1] = imag( integral );
    }
  }
  
  m_emptyNormIntCache = false;
  m_emptyAmpIntCache = false;
  
  return input;
}

void
NormIntInterface::operator+=( const NormIntInterface& nii )
{
  
  double nAccEvts = nii.numAccEvents();
  double nGenEvts = nii.numGenEvents();
  
  double totalAccEvts = nAccEvts + m_nAccEvents;
  double totalGenEvts = nGenEvts + m_nGenEvents;
  
  string ampName, conjAmpName;
  
  int n = m_termNames.size();
  
  for( int i = 0; i < n; ++i ){
    for( int j = 0; j < n; ++j ){
    
      m_ampIntCache[2*i*n+j]   *= m_nAccEvents;
      m_ampIntCache[2*i*n+j+1] *= m_nAccEvents;

      complex< double > ai = nii.ampInt( m_termNames[i], m_termNames[j]  );
      
      m_ampIntCache[2*i*n+j]   += nAccEvts * real( ai );
      m_ampIntCache[2*i*n+j+1] += nAccEvts * imag( ai );

      m_ampIntCache[2*i*n+j]   /= totalAccEvts;
      m_ampIntCache[2*i*n+j+1] /= totalAccEvts;

      
      m_ampIntCache[2*i*n+j]   *= m_nGenEvents;
      m_ampIntCache[2*i*n+j+1] *= m_nGenEvents;
      
      complex< double > ni = nii.normInt( m_termNames[i], m_termNames[j]  );
      
      m_ampIntCache[2*i*n+j]   += nGenEvts * real( ni );
      m_ampIntCache[2*i*n+j+1] += nGenEvts * imag( ni );
      
      m_ampIntCache[2*i*n+j]   /= totalGenEvts;
      m_ampIntCache[2*i*n+j+1] /= totalGenEvts;
      
    }
  }
  
  
  m_nAccEvents = totalAccEvts;
  m_nGenEvents = totalGenEvts;
  
  m_emptyNormIntCache = false;
  m_emptyAmpIntCache = false;
}

bool
NormIntInterface::hasAccessToMC() const
{
  
  return( ( m_accMCReader != NULL ) &&
          ( m_genMCReader != NULL ) &&
          ( m_pIntenManager != NULL ) );
}

bool
NormIntInterface::hasNormInt( string amp, string conjAmp ) const
{
  
  return( ( m_termIndex.find( amp )     != m_termIndex.end() ) &&
          ( m_termIndex.find( conjAmp ) != m_termIndex.end() ) );
}


complex< double > 
NormIntInterface::normInt( string amp, string conjAmp, bool forceUseCache ) const
{
  
  if( !hasNormInt( amp, conjAmp ) ){
    
    cout << "ERROR: normalization integral does not exist for "
    << amp << ", " << conjAmp << "*" << endl;
    assert( false );
  }

  // pass in true flag to delay computation of amplitude integrals as
  // long as possible
  
  if( m_emptyNormIntCache ) forceCacheUpdate( true );
  
  if( !hasAccessToMC() && ( m_pIntenManager != NULL ) &&
       m_pIntenManager->hasTermWithFreeParam() ){
    
    cout << "ERROR: the AmplitudeManager has amplitudes that contain free\n"
    << "       parameters, but no MC has been provided to recalculate\n"
    << "       the normalization integrals.  Check that config file\n"
    << "       lists appropriate data and MC sources." << endl;
    
    assert( false );
  }
  
  if( !forceUseCache && hasAccessToMC() && 
     ( m_pIntenManager != NULL ) && m_pIntenManager->hasTermWithFreeParam() ){
    
    m_pIntenManager->calcIntegrals( m_mcVecs, m_nGenEvents );
    setNormIntMatrix( m_mcVecs.m_pdIntegralMatrix );
  }
  
  int n = m_termNames.size();
  int i = m_termIndex.find( amp )->second;
  int j = m_termIndex.find( conjAmp )->second;
  
  return complex< double >( m_normIntCache[2*i*n+2*j],
                            m_normIntCache[2*i*n+2*j+1] );
}


bool
NormIntInterface::hasAmpInt( string amp, string conjAmp ) const
{
  
  // it is not possible to have one and not the other anymore
  
  return hasNormInt( amp, conjAmp );
}


complex< double > 
NormIntInterface::ampInt( string amp, string conjAmp, bool forceUseCache ) const
{
  
  if( m_emptyAmpIntCache ) forceCacheUpdate();
  
  if( !forceUseCache && ( m_pIntenManager != NULL ) &&
       m_pIntenManager->hasTermWithFreeParam() ){
    
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
  
  
  if( !hasAmpInt( amp, conjAmp ) ){
    
    cout << "ERROR: amplitude integral does not exist for "
    << amp << ", " << conjAmp << "*" << endl;
    assert( false );
  }
  
  int n = m_termNames.size();
  int i = m_termIndex.find( amp )->second;
  int j = m_termIndex.find( conjAmp )->second;
  
  return complex< double >( m_ampIntCache[2*i*n+2*j],
                            m_ampIntCache[2*i*n+2*j+1] );
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
    
    m_pIntenManager->calcIntegrals( m_mcVecs, m_nGenEvents );
    setNormIntMatrix( m_mcVecs.m_pdIntegralMatrix );
    
    m_emptyNormIntCache = false;
    
    // stop here if we just want the norm ints -- this will likely be the normal
    // mode of operation for NI recalculations during a fit
    return;
  }
  
  if( !normIntOnly ){
  
    // now we need to have the generated MC in addition to the accepted
    // MC in order to be able to continue
    assert( m_genMCReader != NULL );
  
    // flush this if anything is loaded
    m_mcVecs.deallocAmpVecs();
  
    cout << "Loading generated Monte Carlo..." << endl;
    m_mcVecs.loadData( m_genMCReader );
    m_mcVecs.allocateTerms( *m_pIntenManager );
    cout << "\tDone.\n" << flush;
  
    cout << "Calculating integrals..." << endl;
    m_pIntenManager->calcIntegrals( m_mcVecs, m_nGenEvents );
    setAmpIntMatrix( m_mcVecs.m_pdIntegralMatrix );
    cout << "\tDone." << endl;
  
    m_emptyAmpIntCache = false;
  }
  
  if( ( m_accMCReader == m_genMCReader ) && !m_emptyAmpIntCache )
  {    
    // optimization for perfect acceptance  
    cout << "Perfect acceptance -- using integrals from generated MC" << endl;
    
    setNormIntMatrix( m_ampIntCache );
  }
  else
  {  
    // load the accepted MC -- always leave accepted MC in memory
    m_mcVecs.deallocAmpVecs();
    
    cout << "Loading acccepted Monte Carlo..." << endl;
    m_mcVecs.loadData( m_accMCReader );
    m_mcVecs.allocateTerms( *m_pIntenManager );
    cout << "\tDone." << endl;
    
    // since we have flushed the amplitudes from AmpVecs we need
    // to recalcualte them or else subsequent calls to calcIntegrals
    // with firstPass set to false will fail
    cout << "Calculating integrals..." << endl;
    m_pIntenManager->calcIntegrals( m_mcVecs, m_nGenEvents );
    setNormIntMatrix( m_mcVecs.m_pdIntegralMatrix );
    
    cout << "\tDone." << endl;
  }
  
  m_emptyNormIntCache = false;
  
  // this method always leaves class in state with m_mcVecs holding accepted MC
}  

void
NormIntInterface::exportNormIntCache( const string& fileName, bool renormalize) const
{
  ofstream out( fileName.c_str() );
  out.precision( 15 );
  exportNormIntCache( out, renormalize );
}

void
NormIntInterface::exportNormIntCache( ostream& out, bool renormalize) const
{
  if( m_emptyNormIntCache || m_emptyAmpIntCache ) forceCacheUpdate();
  
  out << m_nGenEvents << "\t" << m_nAccEvents << endl;
  
  out << m_termNames.size() << endl;
  
  for( vector< string >::const_iterator name = m_termNames.begin();
       name != m_termNames.end();
       ++name ){
    
    out << *name << endl;
  }
  
  int n = m_termNames.size();

  // write out generated matrix first
  for( int i = 0; i < n ; ++i ){
    for( int j = 0; j < n; ++j ) {
    
      complex< double > value( m_ampIntCache[2*i*n+2*j],
                               m_ampIntCache[2*i*n+2*j+1] );
      
      if( renormalize ){
 
        // diagonal elements are real so we don't need to deal
        // with the imaginary part in the renormalization
        
        value /= sqrt( m_ampIntCache[2*i*n+2*i] *
                       m_ampIntCache[2*j*n+2*j] );
        
      }
      
      out << value << "\t";
    }

    out << endl;
  }
  
  // then accepted matrix
  for( int i = 0; i < n ; ++i ){
    for( int j = 0; j < n; ++j ) {
      
      complex< double > value( m_normIntCache[2*i*n+2*j],
                               m_normIntCache[2*i*n+2*j+1] );
      
      if( renormalize ){
        
        // diagonal elements are real so we don't need to deal
        // with the imaginary part in the renormalization
        // renormalize by generated so one has efficiency on
        // the diagonal for the accepted matrix
        
        value /= sqrt( m_ampIntCache[2*i*n+2*i] *
                       m_ampIntCache[2*j*n+2*j] );
        
      }
      
      out << value << "\t";
    }
    
    out << endl;
  }
}

void
NormIntInterface::initializeCache() {
  
  int n = m_termNames.size();
  
  for( int i = 0; i < n; ++i )
    m_termIndex[m_termNames[i]] = i;
  
  // 2 for real and imaginary
  // allocate for the whole matrix even though it is square under
  // complex conjugation so that memory lookups are easier
  m_cacheSize = 2*n*n;
  
  m_normIntCache = new GDouble[m_cacheSize];
  m_ampIntCache = new GDouble[m_cacheSize];
  
  memset( m_normIntCache, 0, m_cacheSize * sizeof( GDouble ) );
  memset( m_ampIntCache, 0, m_cacheSize * sizeof( GDouble ) );
}

void
NormIntInterface::setAmpIntMatrix( const GDouble* input ) const {
  
  memcpy( m_ampIntCache, input, m_cacheSize * sizeof( GDouble ) );
  m_emptyAmpIntCache = false;
}

void
NormIntInterface::setNormIntMatrix( const double* input ) const {
  
  memcpy( m_normIntCache, input, m_cacheSize * sizeof( GDouble ) );
  m_emptyNormIntCache = false;
}
