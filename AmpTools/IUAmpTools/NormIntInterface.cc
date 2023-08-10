
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
#include "IUAmpTools/report.h"
const char* NormIntInterface::kModule = "NormIntInterface";

#ifndef __ACLIC__
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/DataReader.h"
#endif

NormIntInterface::NormIntInterface() :
#ifndef __ACLIC__
m_pIntenManager( NULL ),
m_accMCReader( NULL ),
m_genMCReader( NULL ),
#endif
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{}

NormIntInterface::NormIntInterface( const string& normIntFile ) :
#ifndef __ACLIC__
m_pIntenManager( NULL ),
m_accMCReader( NULL ),
m_genMCReader( NULL ),
#endif
m_emptyNormIntCache( true ),
m_emptyAmpIntCache( true )
{
  
  report( INFO, kModule ) << "Reading cached normalization integral calculation from: "
       << normIntFile << endl;
  
  ifstream inFile( normIntFile.c_str() );
  
  if( !inFile ){
    report( ERROR, kModule ) << "Could not find file " << normIntFile << endl;
    assert( false );
  }
  
  inFile >> (*this);
}

#ifndef __ACLIC__
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
  
  std::map<DataReader*,AmpVecs*>::iterator genVecs = m_uniqueDataSets.find( m_genMCReader );
  if( genVecs == m_uniqueDataSets.end() ){
  
    report( INFO, kModule ) << "Loading generated Monte Carlo from file..." << endl;
    m_genMCVecs.loadData( m_genMCReader );
    
    m_uniqueDataSets[m_genMCReader] = &m_genMCVecs;
  }
  else{
    
    report( NOTICE, kModule ) << "Duplicated Monte Carlo set detected, "
         << "using previously loaded version" << endl;
    
    genVecs->second->shareDataWith( &m_genMCVecs );
  }
  m_nGenEvents = m_genMCVecs.m_iNTrueEvents;

  std::map<DataReader*,AmpVecs*>::iterator accVecs = m_uniqueDataSets.find( m_accMCReader );
  if( accVecs == m_uniqueDataSets.end() ){
  
    report( INFO, kModule ) << "Loading accepted Monte Carlo from file..." << endl;
    m_accMCVecs.loadData( m_accMCReader );
    
    m_uniqueDataSets[m_accMCReader] = &m_accMCVecs;
  }
  else{
    
    report( NOTICE, kModule ) << "Duplicated Monte Carlo set detected, "
         << "using previously loaded version" << endl;
    
    accVecs->second->shareDataWith( &m_accMCVecs );
  }
  m_sumAccWeights = m_accMCVecs.m_dSumWeights;
  
  initializeCache();
}

NormIntInterface::~NormIntInterface(){

  // Find if this interface's AmpVecs are being used
  // in the unique data set map and delete them from the
  // map if they are.  When the objects go out of scope
  // the destructor of AmpVecs will properly transfer
  // ownership of the data to one of the shared objects.
  
  std::map<DataReader*,AmpVecs*>::iterator genVecs = m_uniqueDataSets.find( m_genMCReader );
  if( genVecs != m_uniqueDataSets.end() &&
     genVecs->second == &m_genMCVecs ){
    
    m_uniqueDataSets.erase( genVecs );
  }

  std::map<DataReader*,AmpVecs*>::iterator accVecs = m_uniqueDataSets.find( m_accMCReader );
  if( accVecs != m_uniqueDataSets.end() &&
     accVecs->second == &m_accMCVecs ){
    
    m_uniqueDataSets.erase( accVecs );
  }
}

#endif

istream&
NormIntInterface::loadNormIntCache( istream& input )
{
  input >> m_nGenEvents >> m_sumAccWeights;
  
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
  
  double totalAccEvts = nAccEvts + m_sumAccWeights;
  double totalGenEvts = nGenEvts + m_nGenEvents;
  
  string ampName, conjAmpName;
  
  int n = m_termNames.size();
  
  for( int i = 0; i < n; ++i ){
    for( int j = 0; j < n; ++j ){
    
      m_ampIntCache[2*i*n+j]   *= m_sumAccWeights;
      m_ampIntCache[2*i*n+j+1] *= m_sumAccWeights;

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
  
  
  m_sumAccWeights = totalAccEvts;
  m_nGenEvents = totalGenEvts;
  
  m_emptyNormIntCache = false;
  m_emptyAmpIntCache = false;
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
    
    report( ERROR, kModule ) << "ERROR: normalization integral does not exist for "
    << amp << ", " << conjAmp << "*" << endl;
    assert( false );
  }

#ifdef __ACLIC__
  // dummy line to avoid compiler warning in ROOT applications
  if( forceUseCache ) report( DEBUG, kModule ) << endl;
#else
  // pass in true flag to delay computation of amplitude integrals as
  // long as possible
  
  if( m_emptyNormIntCache ) forceCacheUpdate( true );
  
  if( !hasAccessToMC() && ( m_pIntenManager != NULL ) &&
       m_pIntenManager->hasTermWithFreeParam() ){
    
    report( ERROR, kModule ) << "The AmplitudeManager has amplitudes that contain free\n"
    << "\tparameters, but no MC has been provided to recalculate\n"
    << "\tthe normalization integrals.  Check that config file\n"
    << "\tlists appropriate data and MC sources." << endl;
    
    assert( false );
  }
  
  if( !forceUseCache && hasAccessToMC() && 
     ( m_pIntenManager != NULL ) && m_pIntenManager->hasTermWithFreeParam() ){
    
    m_pIntenManager->calcIntegrals( m_accMCVecs, m_nGenEvents );
    setNormIntMatrix( m_accMCVecs.m_pdIntegralMatrix );
  }
#endif
  
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
  
#ifdef __ACLIC__
  // dummy line to avoid compiler warning in ROOT applications
  if( forceUseCache ) report( DEBUG, kModule ) << endl;
#else
  if( m_emptyAmpIntCache ) forceCacheUpdate();
  
  if( !forceUseCache && ( m_pIntenManager != NULL ) &&
       m_pIntenManager->hasTermWithFreeParam() ){
    
    report( WARNING, kModule ) << "Request of for numerical integral of amplitude\n"
    << "\tthat has floating parameters using *generated* MC.\n"
    << "\t(This is not typically used in a fit -- check for bugs!)\n"
    << "\tProviding original cached value which may be incorrect if\n"
    << "\tthe parameter has changed since initialization."
    << endl;
  }
#endif

  if( !hasAmpInt( amp, conjAmp ) ){
    
    report( ERROR, kModule ) << "amplitude integral does not exist for "
    << amp << ", " << conjAmp << "*" << endl;
    assert( false );
  }
  
  int n = m_termNames.size();
  int i = m_termIndex.find( amp )->second;
  int j = m_termIndex.find( conjAmp )->second;
  
  return complex< double >( m_ampIntCache[2*i*n+2*j],
                            m_ampIntCache[2*i*n+2*j+1] );
}

#ifndef __ACLIC__
bool
NormIntInterface::hasAccessToMC() const
{
  
  return( ( m_accMCReader != NULL ) &&
          ( m_genMCReader != NULL ) &&
          ( m_pIntenManager != NULL ) );
}

void
NormIntInterface::forceCacheUpdate( bool normIntOnly ) const
{
  
  // if the accepted MC is not available, then the data have likely
  // not been loaded into AmpVecs and we can't reclaculate the integrals
  // below
  assert( m_accMCVecs.m_dataLoaded );
  
  // do "lazy" allocation of memory here -- this is important for MPI jobs
  // where forceCacheUpdate is only called on follower nodes, as it
  // avoids big memory allocations on the lead nodes
  if( m_accMCVecs.m_iNTerms == 0 ) m_accMCVecs.allocateTerms( *m_pIntenManager );
  
  // we won't enter here if the cache is empty, which can happen
  // on the first pass through the data -- for subsequent passes
  // (during fitting) the loop below will execute and return
  if( !m_emptyNormIntCache && normIntOnly ){
        
    m_pIntenManager->calcIntegrals( m_accMCVecs, m_nGenEvents );
    setNormIntMatrix( m_accMCVecs.m_pdIntegralMatrix );
    
    m_emptyNormIntCache = false;
    
    // stop here if we just want the norm ints -- this will likely be the normal
    // mode of operation for NI recalculations during a fit
    return;
  }
  
  // if we make it this far it is the first trip through or a forced
  // update to the cache...
  // we only need to recalculate the generated MC integrals if they are
  // empty or there are floating paramters in the amplitude
  if( !normIntOnly &&
      ( m_emptyAmpIntCache || m_pIntenManager->hasTermWithFreeParam() ) ){
      
    // now we need to have the generated MC in addition to the accepted
    // MC in order to be able to continue
    assert( m_genMCVecs.m_dataLoaded );

    // do "lazy" allocation of memory here -- this is important for MPI jobs
    // where forceCacheUpdate is only called on follower nodes, as it
    // avoids big memory allocations on the lead nodes
    if( m_genMCVecs.m_iNTerms == 0 ) m_genMCVecs.allocateTerms( *m_pIntenManager );
    
    m_pIntenManager->calcIntegrals( m_genMCVecs, m_nGenEvents );
    setAmpIntMatrix( m_genMCVecs.m_pdIntegralMatrix );
  
    m_emptyAmpIntCache = false;
  }
  
  // first trip through or forced update...
  // we either copy the generated MC integrals to accepted or
  // compute thm directly:
  if( ( m_accMCReader == m_genMCReader ) && !m_emptyAmpIntCache ) {
    
    // optimization for perfect acceptance  
    report( NOTICE, kModule ) << "Perfect acceptance -- generated and accepted MC are the same" << endl;
    
    setNormIntMatrix( m_ampIntCache );
  }
  else {
    
    m_pIntenManager->calcIntegrals( m_accMCVecs, m_nGenEvents );
    setNormIntMatrix( m_accMCVecs.m_pdIntegralMatrix );
  }
  m_emptyNormIntCache = false;
}

#endif

void
NormIntInterface::exportNormIntCache( const string& fileName ) const
{
  ofstream out( fileName.c_str() );
  out.precision( 15 );
  exportNormIntCache( out );
}

void
NormIntInterface::exportNormIntCache( ostream& out ) const
{
  
#ifndef __ACLIC__
  if( m_emptyNormIntCache || m_emptyAmpIntCache ) forceCacheUpdate();
#endif
  
  out << m_nGenEvents << "\t" << m_sumAccWeights << endl;
  
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

      out << value << "\t";
    }

    out << endl;
  }
  
  // then accepted matrix
  for( int i = 0; i < n ; ++i ){
    for( int j = 0; j < n; ++j ) {
      
      complex< double > value( m_normIntCache[2*i*n+2*j],
                               m_normIntCache[2*i*n+2*j+1] );
            
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
  
  m_normIntCache = new double[m_cacheSize];
  m_ampIntCache = new double[m_cacheSize];
  
  memset( m_normIntCache, 0, m_cacheSize * sizeof( double ) );
  memset( m_ampIntCache, 0, m_cacheSize * sizeof( double ) );
}

void
NormIntInterface::setAmpIntMatrix( const double* input ) const {
  
  memcpy( m_ampIntCache, input, m_cacheSize * sizeof( double ) );
  m_emptyAmpIntCache = false;
}

void
NormIntInterface::setNormIntMatrix( const double* input ) const {
  
  memcpy( m_normIntCache, input, m_cacheSize * sizeof( double ) );
  m_emptyNormIntCache = false;
}

#ifndef __ACLIC__
map< DataReader*, AmpVecs* > NormIntInterface::m_uniqueDataSets;
#endif

