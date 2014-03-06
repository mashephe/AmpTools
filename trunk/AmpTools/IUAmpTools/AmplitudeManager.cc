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
#include <string.h>

#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/NormIntInterface.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

AmplitudeManager::AmplitudeManager( const vector< string >& reaction,
                                    const string& reactionName) :
IntensityManager( reaction, reactionName ),
m_optimizeParIteration( false )
{
  cout << endl << "## AMPLITUDE MANAGER INITIALIZATION ##" << endl;
  cout << " Creating amplitude manager for reaction:  " << reactionName << endl;
  
  // a list of index switches needed to generate the symmetrized amplitude
  // group the switches by particle type
  // dump out some information
  map< string, vector< pair< int, int > > > swapsByType;
  for( unsigned int i = 0; i < reaction.size(); ++i ){
    
    cout << "\t particle index assignment:  " << reaction[i] << " -->> " << i << endl;
    
    for( unsigned int j = i + 1; j < reaction.size(); ++j ){
      
      if( reaction[i] == reaction[j] ){
        
        swapsByType[reaction[i]].push_back( pair< int, int >( i, j ) );
      }
    }
  }
  
  // determine how many combinations exist by taking the the product of the
  // number of swaps for each unique final state particle
  int numberOfCombos = 1;
  for( map< string, vector< pair< int, int > > >::iterator partItr = swapsByType.begin();
      partItr != swapsByType.end(); ++partItr ){
    
    // don't forget the option of leaving the particles unchanged
    int partCombos = partItr->second.size() + 1;
    
    cout << "There are " << partCombos
    << " ways of rearranging particles of type: "
    << partItr->first << endl;
    
    numberOfCombos *= partCombos;
  }
  
  // setup the vector of symmetric combinations -- first initialize
  // it with numberOfCombos copies of the default ordering
  // then go in and make the swaps
  vector< int > defaultOrder( reaction.size() );
  for( unsigned int i = 0; i < reaction.size(); ++i ){
    
    defaultOrder[i] = i;
  }
  
  // strip away the particle names to build a vector of vector of swaps
  // analagous to the map from particle -> vector of swaps
  vector< vector< pair< int, int > > > swaps;
  for( map< string, vector< pair< int, int > > >::iterator itr = swapsByType.begin();
      itr != swapsByType.end(); ++itr ){
    
    swaps.push_back( itr->second );
    
    // for each type of particle, add a swap that leaves particles
    // as they are in the default ordering
    // (this switches particle 0 with 0 which leaves order unchanged
    // for any particle in consideration)
    swaps.back().push_back( pair< int, int >( 0, 0 ) );
  }
  
  // now use a recursive algorithm to step through the list of swaps for
  // each unique final state particle and add to the vector of symmetric
  // combinations
  generateSymmetricCombos( vector< pair< int, int > >( 0 ),
                          swaps, defaultOrder );
  
  if( m_symmCombos.size() > 1 ){
    
    cout << "The following " << numberOfCombos << " orderings of the particles are" << endl
    << "indistinguishable and will be permuted when computing amplitudes." << endl;
    
    for( unsigned int i = 0; i < m_symmCombos.size(); ++i ){
      
      for( unsigned int j = 0; j < reaction.size(); ++j ){
        
        cout << "\t" << m_symmCombos[i][j];
      }
      cout << endl;
    }
  }
}

AmplitudeManager::~AmplitudeManager() {
  
  for( map< string, Amplitude* >::iterator mapItr = m_registeredFactors.begin();
      mapItr != m_registeredFactors.end(); ++mapItr ){
    
    delete mapItr->second;
  }
  
  for( map< string, vector< const Amplitude* > >::iterator mapItr = m_mapNameToAmps.begin();
      mapItr != m_mapNameToAmps.end(); ++mapItr ){
    
    for( vector< const Amplitude* >::iterator vecItr = mapItr->second.begin();
        vecItr != mapItr->second.end(); ++vecItr ){
      
      delete *vecItr;
    }
  }
  
}

int
AmplitudeManager::termFactorStoragePerEvent() const {
  
  vector< string > ampNames = getTermNames();
  
  int nAmps = ampNames.size();

  int nAmpFactorsAndPerms = 0;
  
  for( int i = 0; i < nAmps; i++ )
  {
    int iNPermutations = getPermutations( ampNames[i] ).size();
    int iNFactors = getFactors( ampNames[i] ).size();
    
    assert( iNPermutations*iNFactors );
    
    nAmpFactorsAndPerms += iNPermutations*iNFactors;
  }
  
  // for each factor and permutation we need to store
  // a complex number, which is two numbers
  
  return 2 * nAmpFactorsAndPerms;
}

int
AmplitudeManager::termStoragePerEvent() const {
  
  // for each amplitude we need to store a complex
  // number -- that is two numbers
  
  return 2 * getTermNames().size();
}


bool
AmplitudeManager::hasTermWithFreeParam() const {
  
  for( vector< bool >::const_iterator isFixed = m_vbIsAmpFixed.begin();
      isFixed != m_vbIsAmpFixed.end();
      ++isFixed ){
    
    if( !(*isFixed) ) return true;
  }
  
  return false;
}

void
AmplitudeManager::calcTerms( AmpVecs& a, bool bIsFirstPass, bool useMC ) const
{
  
#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcTerms" );
#endif

  timeval tStart, tIStart, tStop,tSpan;
  double dTime;
  
  if( bIsFirstPass ) gettimeofday( &(tIStart), NULL );
  
#ifdef GPU_ACCELERATION
  if( bIsFirstPass ) initGPU( a, useMC );
  GPUManager& gpuMan = ( useMC ? m_mcGPUManGTX : m_dataGPUManGTX );
#endif
  
  vector< string > ampNames = getTermNames();
  
  int iNAmps = ampNames.size();
  
  assert( iNAmps && a.m_iNEvents && a.m_iNTrueEvents );
  assert( a.m_pdAmps && a.m_pdAmpFactors);
  
  int iAmpIndex, iAmpFactOffset=0;
  for( iAmpIndex = 0; iAmpIndex < iNAmps; iAmpIndex++ )
  {
    
    map< string, vector< vector< int > > >::const_iterator permItr =
    m_ampPermutations.find( ampNames[iAmpIndex] );
    assert( permItr != m_ampPermutations.end() );
    const vector< vector< int > >& vvPermuations = permItr->second;
    int iNPermutations = vvPermuations.size();
    
    vector< const Amplitude* > vAmps =
    m_mapNameToAmps.find(ampNames.at(iAmpIndex))->second;
    
    int iLocalOffset = 0;
    int iFactor, iNFactors = vAmps.size();
    
    // if it is not the first pass through and this particular
    // amplitude is fixed, then we can skip to the next amplitude
    // and avoid all of the symmetrization computation as well
    // as checking each individual factor for free parameters.
    // HOWEVER we need to be sure to adjust the iAmpFactOffset correctly
    // so it is ready for the next amplitude.
    if( !bIsFirstPass && m_vbIsAmpFixed[iAmpIndex] ){
      
      iAmpFactOffset += 2 * a.m_iNEvents * iNPermutations * iNFactors;
      continue;
    }
    
    // cout << "Recomputing amplitude: " << m_ampNames[iAmpIndex] << endl;
    
    // calculate all the factors that make up an amplitude for
    // for all events serially on CPU or in parallel on GPU
    
    const Amplitude* pCurrAmp = 0;
    for( iFactor=0; iFactor < iNFactors;
        iFactor++, iLocalOffset += 2 * a.m_iNEvents * iNPermutations ){
      
      pCurrAmp = vAmps.at( iFactor );
      
      if( !bIsFirstPass && !pCurrAmp->containsFreeParameters() )
        continue;
      
      // now check to see if the value of this amplitudes parameters
      // are the same as they were the last time they were evaluated
      // for this particular dataset -- if so don't reevaluate
      
      if( !bIsFirstPass && m_optimizeParIteration &&
          m_dataAmpIteration[&a][pCurrAmp] == m_ampIteration[pCurrAmp] )
        continue;
      
      m_dataAmpIteration[&a][pCurrAmp] = m_ampIteration[pCurrAmp];
      
      // cout << "Recomputing factor: " << pCurrAmp->name() << endl;
      
      if( bIsFirstPass ) gettimeofday( &(tStart), NULL );
      
#ifndef GPU_ACCELERATION
      pCurrAmp->
      calcAmplitudeAll( a.m_pdData,
                       a.m_pdAmpFactors + iAmpFactOffset + iLocalOffset,
                       a.m_iNEvents, &vvPermuations );
#else
      gpuMan.calcAmplitudeAll( pCurrAmp,
                              a.m_pdAmpFactors + iAmpFactOffset + iLocalOffset,
                              &vvPermuations );
#endif//GPU_ACCELERATION
      
      if( bIsFirstPass ){
        gettimeofday( &(tStop), NULL );
        timersub( &(tStop), &(tStart), &tSpan );
        dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
        
        //cout << "-> Seconds spent calculating "
        //<< pCurrAmp->name() << " for " << permItr->first << ":  "
        //<< dTime << endl << flush;
      }
            
    }
    
    // now assemble all the factors in an amplitude into a single
    // symmetrized amplitude -- do this for all events for this
    // amplitude (on the CPU)
    
    GDouble dSymmFactor = 1.0f/sqrt( iNPermutations );
    GDouble dAmpFacRe, dAmpFacIm, dTRe, dTIm;
    int iEvent, iPerm;
    int iOffsetA, iOffsetP, iOffsetF;
    
    // re-ordering of data will be useful to not fall out of (CPU) memory cache!!!
    
    // zeroing out the entire range
    memset( (void*)( a.m_pdAmps + 2 * a.m_iNEvents * iAmpIndex ), 0,
           2 * a.m_iNEvents * sizeof(GDouble) );
    
    // only sum over the true events from data and skip paddings
    for( iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ )
    {
      iOffsetA = 2 * a.m_iNEvents * iAmpIndex + 2 * iEvent;
      
      for( iPerm = 0; iPerm < iNPermutations; iPerm++ )
      {
        iOffsetP = iAmpFactOffset + 2 * a.m_iNEvents * iPerm + 2 * iEvent;
        
        dAmpFacRe = a.m_pdAmpFactors[iOffsetP];
        dAmpFacIm = a.m_pdAmpFactors[iOffsetP+1];
        
        for( iFactor = 1; iFactor < iNFactors; iFactor++ )
        {
          iOffsetF = iOffsetP + 2 * a.m_iNEvents * iNPermutations * iFactor;
          
          dTRe = dAmpFacRe;
          dTIm = dAmpFacIm;
          
          dAmpFacRe = dTRe * a.m_pdAmpFactors[iOffsetF] -
          dTIm * a.m_pdAmpFactors[iOffsetF+1];
          dAmpFacIm = dTRe * a.m_pdAmpFactors[iOffsetF+1] +
          dTIm * a.m_pdAmpFactors[iOffsetF];
        }
        
        a.m_pdAmps[iOffsetA]   += dAmpFacRe;
        a.m_pdAmps[iOffsetA+1] += dAmpFacIm;
      }
      
      a.m_pdAmps[iOffsetA]   *= dSymmFactor;
      a.m_pdAmps[iOffsetA+1] *= dSymmFactor;
    }
    
    //Update the global offset of ampfactor array
    iAmpFactOffset += iLocalOffset;
  }
  
  if( bIsFirstPass ){
    
    gettimeofday( &(tStop), NULL );
    
    timersub( &(tStop), &(tIStart), &tSpan );
    dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
    //cout << "\t-->> Seconds spent calculating all amplitudes:  " << dTime << endl;
  }
  
}

double
AmplitudeManager::calcIntensities( AmpVecs& a, bool bIsFirstPass ) const
{

#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcIntensities" );
#endif

  // check to be sure destination memory has been allocated
  assert( a.m_pdIntensity );
  
  double maxInten = 0;
  
  //First update Amplitudes if needed
  calcTerms( a, bIsFirstPass );
  
  vector< string > ampNames = getTermNames();
  
  int iNAmps = ampNames.size();
  
  //Now pre-calculate ViVj* and include factor of 2 for off-diagonal elements
  double* pdViVjRe = new double[iNAmps*(iNAmps+1)/2];
  double* pdViVjIm = new double[iNAmps*(iNAmps+1)/2];
  
  int i,j;
  complex< double > cTmp;
  for( i = 0; i < iNAmps; i++ ){
    for( j = 0; j <= i; j++ ){
      
      cTmp = productionFactor( i ) * conj( productionFactor( j ) );
            
      if( termsAreRenormalized() ){
        
        cTmp /= sqrt( normInt()->ampInt( ampNames[i], ampNames[i] ) *
                      normInt()->ampInt( ampNames[j], ampNames[j] ) );
      }
      
      pdViVjRe[i*(i+1)/2+j] = cTmp.real();
      pdViVjIm[i*(i+1)/2+j] = cTmp.imag();
      
      if( i != j ) {
        
        pdViVjRe[i*(i+1)/2+j] *= 2;
        pdViVjIm[i*(i+1)/2+j] *= 2;
      }
    }
  }
  
  double dIntensity;
  double cAiAjRe,cAiAjIm;
  
  //Re-ordering of data will be useful to not fall out of (CPU) memory cache!!!
  //Only sum over the true events from data and skip paddings
  int iEvent;
  for( iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ )
  {
    dIntensity = 0;
    for( i = 0; i < iNAmps; i++ ){
      for( j = 0; j <= i; j++ ){
        
        // remove cross terms from incoherent amplitudes
        if( !m_sumCoherently[i][j] ) continue;
        
        //AiAj*
        cAiAjRe = a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent] *
        a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent] +
        a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent+1] *
        a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent+1];
        
        cAiAjIm= -a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent] *
        a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent+1] +
        a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent+1] *
        a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent];
        
        dIntensity += pdViVjRe[i*(i+1)/2+j] * cAiAjRe -
        pdViVjIm[i*(i+1)/2+j] * cAiAjIm;
      }
    }
    
    dIntensity *= a.m_pdWeights[iEvent];
    
    a.m_pdIntensity[iEvent] = dIntensity;
    if( dIntensity > maxInten ) maxInten = dIntensity;
  }
  
  delete[] pdViVjRe;
  delete[] pdViVjIm;
  
  return maxInten;
}


double
AmplitudeManager::calcSumLogIntensity( AmpVecs& a, bool bIsFirstPass ) const
{
  
#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcSumLogIntensity" );
#endif

  // this may be inefficienct since there are two
  // loops over events, one here and one in the
  // calculation of intensities -- however, this
  // streamlines the code a little
  // this may be a place for optimization later
  
  double dSumLogI = 0;
  
#ifndef GPU_ACCELERATION
  
  calcIntensities( a, bIsFirstPass );
  
  for( int iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ ){
    
    // here divide out the weight that was put into the intensity calculation
    // and weight the log -- in practice this just contributes an extra constant
    // term in the likelihood equal to sum -w_i * log( w_i ), but the division
    // helps avoid problems with negative weights, which may be used
    // in background subtraction
    dSumLogI += a.m_pdWeights[iEvent] *
    G_LOG( a.m_pdIntensity[iEvent] / a.m_pdWeights[iEvent] );
  }
  
#else
  
  // need to compute the production coefficients with all scale factors
  // taken into account
  
  vector< string > ampNames = getTermNames();
  
  vector< complex< double > > gpuProdPars( ampNames.size() );
  
  for( int i = 0; i < ampNames.size(); ++i ){
    
    gpuProdPars[i] = productionFactor( ampNames[i] );
    
    if( termsAreRenormalized() ){
      
      gpuProdPars[i] /= sqrt( normInt()->ampInt( ampNames[i], ampNames[i] ) );
    }
  }
  
  // need to explicitly do amplitude calculation
  // since intensity and sum is done directly on GPU

  if( bIsFirstPass || hasTermWithFreeParam() ){

    calcTerms( a, bIsFirstPass );
    
    // this operation is only done on data
    m_dataGPUManGTX.copyAmpsToGPU( a );
  }
  
  dSumLogI = m_dataGPUManGTX.calcSumLogIntensity( gpuProdPars, m_sumCoherently );
  
#endif
  
  return( dSumLogI );
}


void
AmplitudeManager::calcIntegrals( AmpVecs& a, int iNGenEvents, bool bIsFirstPass ) const
{

#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcIntegrals" );
#endif

  GDouble* integralMatrix = a.m_pdIntegralMatrix;
  
  // this method could be made more efficient by caching a table of
  // integrals associated with each AmpVecs object and then, based on the
  // variables bIsFirstPass and m_vbIsAmpFixed data compute only
  // those terms that could have changed
  
  // amp -> amp* -> value
  assert( iNGenEvents );
  calcTerms( a, bIsFirstPass, true );
  
  int iNAmps = a.m_iNTerms;

  // zero out the matrix to start
  memset( integralMatrix, 0, 2*iNAmps*iNAmps*sizeof( GDouble ) );
  
  int i, j, iEvent;
  for( i = 0; i < iNAmps;i++ )
  {
    
    for( j = 0; j <= i; j++ )
    {
     
      // if the amplitude isn't floating and it isn't the first pass
      // through these data, then its integral can't change
      
      if( !bIsFirstPass && m_vbIsAmpFixed[i] && m_vbIsAmpFixed[j] )
        continue;
      
      // if two amps don't interfere the relevant integral is zero
      if( m_sumCoherently[i][j] ){
        
        for( iEvent = 0; iEvent < a.m_iNTrueEvents; iEvent++ )
        {
          //AiAj*
          integralMatrix[2*i*iNAmps+2*j] += a.m_pdWeights[iEvent] *
          ( a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent] *
           a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent] +
           a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent+1] *
           a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent+1] );
          
          integralMatrix[2*i*iNAmps+2*j+1] += a.m_pdWeights[iEvent] *
          ( -a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent] *
           a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent+1] +
           a.m_pdAmps[2*a.m_iNEvents*i+2*iEvent+1] *
           a.m_pdAmps[2*a.m_iNEvents*j+2*iEvent] );
        }
        
        //Normalize
        integralMatrix[2*i*iNAmps+2*j] /= static_cast< double >( iNGenEvents );
        integralMatrix[2*i*iNAmps+2*j+1] /= static_cast< double >( iNGenEvents );
      }
      
      //Complex conjugate
      if( i != j ) {
        
        integralMatrix[2*j*iNAmps+2*i] = integralMatrix[2*i*iNAmps+2*j];
        integralMatrix[2*j*iNAmps+2*i+1] = -integralMatrix[2*i*iNAmps+2*j+1];
      }
    }
  }
}

const vector< vector< int > >&
AmplitudeManager::getPermutations( const string& name ) const {
  
  map< string, vector< vector< int > > >::const_iterator mapItr =
  m_ampPermutations.find( name );
  
  // check to be sure the amplitude is there:
  assert( mapItr != m_ampPermutations.end() );
  
  return( mapItr->second );
}

const vector< const Amplitude* >&
AmplitudeManager::getFactors( const string& name ) const {
  
  map< string, vector< const Amplitude* > >::const_iterator mapItr =
  m_mapNameToAmps.find( name );
  
  // check to be sure the amplitude is there:
  assert( mapItr != m_mapNameToAmps.end() );
  
  return( mapItr->second );
}

void
AmplitudeManager::addAmpFactor( const string& name,
                                const string& factorName,
                                const vector< string >& args,
                                const string& sum,
                                const string& scale ){
  
  map< string, Amplitude* >::iterator defaultAmp =
     m_registeredFactors.find( factorName );
  
  if( defaultAmp == m_registeredFactors.end() ){
    
    cout << "ERROR: amplitude factor with name " << factorName
    << " has not been registered." << endl;
    assert( false );
  }
  
  Amplitude* newAmp = defaultAmp->second->newAmplitude( args );
  
  // check to see if this is a new term and do some setup if it is
  if( !hasTerm( name ) ){
    
    addTerm( name, scale );
    
    m_ampSum.push_back( sum );
    m_vbIsAmpFixed.push_back( true );
    
    m_mapNameToAmps[name] = vector< const Amplitude* >( 0 );
    
    // check to see if permutations have already been added for this
    // amplitude (before the amplitude was added itself) if so, add
    // the set of permutations that comes from permuting identical
    // particles
    if( m_ampPermutations.find( name ) != m_ampPermutations.end() )
    {
      // permutations have already been added
      for( vector< vector< int > >::iterator vecItr = m_symmCombos.begin();
          vecItr != m_symmCombos.end(); ++vecItr )
      {
        m_ampPermutations[name].push_back( *vecItr );
      }
    }
    else
    {
      // start the set of permutations with those that include
      // just identical particles
      m_ampPermutations[name] = m_symmCombos;
    }
    
    // adjust the matrix that determines which amplitudes add coherently
    // by looking at this sum and other sums
    int nAmps = m_ampSum.size();
    vector< bool > lastRow;
    
    // simultaneously build the last column and last row
    // (since matrix is symmetric)
    for( int i = 0; i < nAmps; ++i ){
      
      bool coh = ( m_ampSum[i] == sum );
      
      if( i < nAmps - 1 ){
        
        m_sumCoherently[i].push_back( coh );
        lastRow.push_back( coh );
      }
    }
    // this is the lower right element on the diagonal -- always true
    lastRow.push_back( true );
    m_sumCoherently.push_back( lastRow );
  }

  m_mapNameToAmps[name].push_back( newAmp );
  
  //Enable a short-cut if no factors are variable in the amplitude
  m_vbIsAmpFixed[termIndex(name)] =
  m_vbIsAmpFixed[termIndex(name)] && !newAmp->containsFreeParameters();
}

void
AmplitudeManager::addAmpPermutation( const string& ampName, const vector< int >& permutation ){
  
  map< string, vector< vector< int > > >::iterator mapItr = m_ampPermutations.find( ampName );
  
  if( mapItr == m_ampPermutations.end() ){
    
    cout << "WARNING:  adding permutation for nonexistent amplitude " << ampName
    << endl;
    
    m_ampPermutations[ampName] = vector< vector< int > >( 0 );
    m_ampPermutations[ampName].push_back( permutation );
    
  }
  else{
    
    bool foundPermutation = false;
    
    for( vector< vector< int > >::const_iterator vecItr = mapItr->second.begin();
        vecItr != mapItr->second.end();
        ++vecItr ){
      
      if( permutation == *vecItr ) foundPermutation = true;
    }
    
    if( !foundPermutation ){
      
      cout << "Adding a new permutation for " << ampName << ":  ";
      for( vector< int >::const_iterator itr = permutation.begin();
          itr != permutation.end();
          ++itr ){
        
        cout << *itr << " ";
      }
      cout << endl;
      
      mapItr->second.push_back( permutation );
    }
    else{
      
      cout << "The permutation ";
      for( vector< int >::const_iterator itr = permutation.begin();
          itr != permutation.end();
          ++itr ){
        
        cout << *itr << " ";
      }
      cout << "already exists for " << ampName << endl;
    }
  }
}

void
AmplitudeManager::setupFromConfigurationInfo( const ConfigurationInfo* configInfo ){
  
  vector< string > sumName;
  
  // loop over amplitudes in the ConfigurationInfo
  vector<AmplitudeInfo*> ampInfoVector = configInfo->amplitudeList(reactionName());
  for (unsigned int i = 0; i < ampInfoVector.size(); i++){
    
    string ampName = ampInfoVector[i]->fullName();
    string sumName = ampInfoVector[i]->sumName();
    string scale   = ampInfoVector[i]->scale();
    
    // add amplitudes
    vector< vector<string> > ampFactors = ampInfoVector[i]->factors();
    for (unsigned int j = 0; j < ampFactors.size(); j++){
      string factorName = ampFactors[j][0];
      vector<string> ampParameters = ampFactors[j];
      ampParameters.erase(ampParameters.begin());
      addAmpFactor( ampName, factorName, ampParameters, sumName, scale );
    }
    
    // add permutations
    vector< vector<int> > permutations = ampInfoVector[i]->permutations();
    for (unsigned int j = 0; j < permutations.size(); j++){
      addAmpPermutation( ampName, permutations[j] );
    }
    
    // add production amplitudes
    setDefaultProductionFactor(ampName, ampInfoVector[i]->value());
    
    // if the amplitude has parameters we should go ahead and set their
    // values, otherwise they will be set to the default value as defined
    // by the AmpParameter class -- in a fit, the ParameterManager will
    // later reset these pointers to point to floating parameters in MINUIT
    vector< ParameterInfo* > pars = ampInfoVector[i]->parameters();
    for( vector< ParameterInfo* >::iterator parItr = pars.begin();
        parItr != pars.end();
        ++parItr ){
      
      setParValue( ampName, (**parItr).parName(), (**parItr).value() );
    }
    
    // finally initialize the amplitudes
    vector< const Amplitude* > ampVec = m_mapNameToAmps[ampName];
    for( vector< const Amplitude* >::iterator amp = ampVec.begin();
        amp != ampVec.end();
        ++amp ) {
      
      // init needs to be non-const or else the user has to deal with
      // mutable data -- in reality we really want some aspects of the
      // amplitude like the name, arguments, etc. to never change and
      // other aspects to be mutable, but this seems to put an extra
      // burden on the user
      
      const_cast< Amplitude* >(*amp)->init();
    }
  }
}

void
AmplitudeManager::setParPtr( const string& name, const string& parName,
                            const double* ampParPtr ){
  
  IntensityManager::setParPtr( name, parName, ampParPtr );
  
  // now look for the parameter as part of the amplitude factors
  
  for( vector< const Amplitude* >::iterator factorItr = m_mapNameToAmps[name].begin();
      factorItr != m_mapNameToAmps[name].end();
      ++factorItr ){
    
    if( (**factorItr).setParPtr( parName, ampParPtr ) ){
      
      m_vbIsAmpFixed[termIndex(name)] = false;
    }
  }
}

void
AmplitudeManager::setParValue( const string& name, const string& parName,
                               double val ){
  
  IntensityManager::setParValue( name, parName, val );
  
  // we will redetermine the status of this variable in the loop below
  
  m_vbIsAmpFixed[termIndex(name)] = true;
  
  // now loop through the amplitude factors looking for the parameter
  
  for( vector< const Amplitude* >::iterator factorItr = m_mapNameToAmps[name].begin();
      factorItr != m_mapNameToAmps[name].end();
      ++factorItr ){
    
    (**factorItr).setParValue( parName, val );
    
    m_vbIsAmpFixed[termIndex(name)] = m_vbIsAmpFixed[termIndex(name)] &&
    !(**factorItr).containsFreeParameters();
  }
}

void
AmplitudeManager::updatePar( const string& parName ) const {
  
  for( map< string, vector< const Amplitude* > >::const_iterator mapItr = m_mapNameToAmps.begin();
      mapItr != m_mapNameToAmps.end();
      ++mapItr ){
    
    for( vector< const Amplitude* >::const_iterator ampItr = mapItr->second.begin();
        ampItr != mapItr->second.end();
        ++ampItr ){
      
      // if we find an amplitude with the parameter update the iteration
      // counter; this may result in multiple increments over one fuction
      // call but that is OK -- iteration numbers just need to be unique
      if( (**ampItr).updatePar( parName ) ){
        
        ++m_ampIteration[*ampItr];
      }
    }
  }
}

void
AmplitudeManager::registerAmplitudeFactor( const Amplitude& amplitude ){
  
  m_registeredFactors[amplitude.name()] = amplitude.clone();
}

// private functions

void
AmplitudeManager::generateSymmetricCombos( const vector< pair< int, int > >& prevSwaps,
                                          vector< vector< pair< int, int > > > remainingSwaps,
                                          const vector< int >& defaultOrder ){
  
  if( remainingSwaps.size() == 0 ){
    
    // we've reached the bottom of the list of swaps
    
    // make the swaps to the default ordering
    vector< int > swappedOrder = defaultOrder;
    for( vector< pair< int, int > >::const_iterator swapItr = prevSwaps.begin();
        swapItr != prevSwaps.end(); ++swapItr ){
      
      int temp = swappedOrder[swapItr->first];
      swappedOrder[swapItr->first] = swappedOrder[swapItr->second];
      swappedOrder[swapItr->second] = temp;
    }
    
    // and push the combination on the list
    m_symmCombos.push_back( swappedOrder );
  }
  else{
    
    // pop off a set of remaining swaps
    vector< pair< int, int > > currentSwaps = remainingSwaps.back();
    remainingSwaps.pop_back();
    
    // loop through them pushing each swap onto the list of prevSwaps
    for( vector< pair< int, int > >::iterator itr = currentSwaps.begin();
        itr != currentSwaps.end(); ++itr ){
      
      vector< pair< int, int > > newPrevSwaps = prevSwaps;
      newPrevSwaps.push_back( *itr );
      generateSymmetricCombos( newPrevSwaps, remainingSwaps, defaultOrder );
    }
  }
}


#ifdef GPU_ACCELERATION
void
AmplitudeManager::initGPU( const AmpVecs& a, bool useMC ) const {
  
  // check to be sure both data is loaded
  // and amplitudes have been allocated
  assert( a.m_pdAmps && a.m_iNTrueEvents );
  
  GPUManager& gpuMan = ( useMC ? m_mcGPUManGTX : m_dataGPUManGTX );
  
  gpuMan.init( a, useMC );
  gpuMan.copyDataToGPU( a );
}
#endif

