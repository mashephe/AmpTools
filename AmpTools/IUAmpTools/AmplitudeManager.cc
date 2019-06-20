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
#include <string.h>
#include <set>

using namespace std;

#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/NormIntInterface.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

AmplitudeManager::AmplitudeManager( const vector< string >& reaction,
                                    const string& reactionName) :
IntensityManager( reaction, reactionName ),
m_needsUserVarsOnly( true ),
m_optimizeParIteration( false ),
m_flushFourVecsIfPossible( false ),
m_forceUserVarRecalculation( false )
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
   
    // this will clean up new amplitudes generated from this amplitude
    delete mapItr->second;
  }
}

unsigned int
AmplitudeManager::maxFactorStoragePerEvent() const {
  
  vector< string > ampNames = getTermNames();
  
  unsigned int nAmpFactorsAndPerms = 0;
  
  for( int i = 0; i < getTermNames().size(); i++ ) {
    
    unsigned int iNPermutations = getPermutations( ampNames[i] ).size();
    unsigned int iNFactors = getFactors( ampNames[i] ).size();
    
    assert( iNPermutations*iNFactors );
    
    if( ( iNPermutations * iNFactors ) > nAmpFactorsAndPerms )
      nAmpFactorsAndPerms = iNPermutations * iNFactors;
  }
  
  // for each factor and permutation we need to store
  // a complex number, which is two numbers
  
  return 2 * nAmpFactorsAndPerms;
}

unsigned int
AmplitudeManager::termStoragePerEvent() const {
  
  // for each amplitude we need to store a complex
  // number -- that is two numbers
  
  return 2 * getTermNames().size();
}

unsigned int
AmplitudeManager::userVarsPerEvent() const {
  
  set< string > countedStaticAmps;
  set< string > countedUniqueAmps;
  
  vector< string > ampNames = getTermNames();
  
  unsigned int userStorage = 0;
  
  for( int i = 0; i < getTermNames().size(); i++ ) {

    unsigned int iNPermutations = getPermutations( ampNames[i] ).size();
    vector< const Amplitude* > factorVec =
      m_mapNameToAmps.find( ampNames[i] )->second;

    for( int j = 0; j < factorVec.size(); ++j ){

      if( factorVec[j]->areUserVarsStatic() ){
       
        // for factors that have static data, we only want to
        // count the allocation once -- if we have counted
        // it already then skip to the next factor
        if( countedStaticAmps.find( factorVec[j]->name() ) ==
            countedStaticAmps.end() )
          countedStaticAmps.insert( factorVec[j]->name() );
        else continue;
      }
      else{
        // user data is not static, so
        // check to see if we have seen an instance of this
        // same amplitude that would behave in the same way
        // (i.e., has the same arguments) and if it exists, we
        // will use user data block from it instead
        if( countedUniqueAmps.find( factorVec[j]->identifier() ) ==
           countedUniqueAmps.end() )
          countedUniqueAmps.insert( factorVec[j]->identifier() );
          else continue;
      }
      
      userStorage += iNPermutations * factorVec[j]->numUserVars();
    }
  }
  
  return userStorage;
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
AmplitudeManager::calcUserVars( AmpVecs& a ) const
{

#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcUserVars" );
#endif

  const vector< string >& ampNames = getTermNames();
  
  int iNAmps = ampNames.size();
  
  assert( iNAmps && a.m_iNEvents && a.m_iNTrueEvents && a.m_pdUserVars);

  int iAmpIndex;
  unsigned long long iUserVarsOffset = 0;
  for( iAmpIndex = 0; iAmpIndex < iNAmps; iAmpIndex++ )
  {
    
    map< string, vector< vector< int > > >::const_iterator permItr =
    m_ampPermutations.find( ampNames[iAmpIndex] );
    assert( permItr != m_ampPermutations.end() );
    const vector< vector< int > >& vvPermuations = permItr->second;
    int iNPerms = vvPermuations.size();
    
    vector< const Amplitude* > vAmps =
    m_mapNameToAmps.find(ampNames.at(iAmpIndex))->second;
    
    int iFactor, iNFactors = vAmps.size();
    
    const Amplitude* pCurrAmp = 0;
    for( iFactor=0; iFactor < iNFactors; iFactor++ ){
      
      pCurrAmp = vAmps.at( iFactor );
      
      if( pCurrAmp->numUserVars() == 0 ) continue;

      // this is the number of variables for the data set
      int iNVars = pCurrAmp->numUserVars();
      int iNData = iNVars * a.m_iNEvents * iNPerms;
      
      // we will set this based on the algorithm below
      unsigned long long thisOffset = 0;

      if( pCurrAmp->areUserVarsStatic() ){
        
        // the user variables are static, so let's look at the
        // list of data pointers stored in the relevant AmpVecs
        // object and see if there is one associated with this
        // amplitude name
        
        map< string, unsigned long long >::const_iterator offsetItr =
        a.m_userVarsOffset.find( pCurrAmp->name() );
        
        if( offsetItr == a.m_userVarsOffset.end() ){
          
          // set the offset to where the calculation will end up
          a.m_userVarsOffset[pCurrAmp->name()] = iUserVarsOffset;
          
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
        a.m_userVarsOffset.find( pCurrAmp->identifier() );
        
        if( offsetItr == a.m_userVarsOffset.end() ){
          
          // set the offset to where the calculation will end up
          a.m_userVarsOffset[pCurrAmp->identifier()] = iUserVarsOffset;
          
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
      
      pCurrAmp->
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

  return;
}


bool
AmplitudeManager::calcTerms( AmpVecs& a ) const
{
  
#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcTerms" );
#endif
  
  // on the first pass through this data set be sure to calculate
  // the user data first, if needed, before doing term calculations
  if( !a.m_termsValid && a.m_userVarsPerEvent > 0 ){
    
    calcUserVars( a );
    if( m_needsUserVarsOnly && m_flushFourVecsIfPossible ) a.clearFourVecs();
  }
  
  bool modifiedTerm = false;
  
  const vector< string >& ampNames = getTermNames();
  
  int iNAmps = ampNames.size();
  
  assert( iNAmps && a.m_iNEvents && a.m_iNTrueEvents );
#ifndef GPU_ACCELERATION
  assert( a.m_pdAmps && a.m_pdAmpFactors);
#endif
  
  int iAmpIndex;
  for( iAmpIndex = 0; iAmpIndex < iNAmps; iAmpIndex++ )
  {
    
    map< string, vector< vector< int > > >::const_iterator permItr =
    m_ampPermutations.find( ampNames[iAmpIndex] );
    assert( permItr != m_ampPermutations.end() );
    const vector< vector< int > >& vvPermuations = permItr->second;
    int iNPermutations = vvPermuations.size();
    
    vector< const Amplitude* > vAmps =
    m_mapNameToAmps.find(ampNames.at(iAmpIndex))->second;
    
    int iFactor, iNFactors = vAmps.size();
    
    // if it is not the first pass through and this particular
    // amplitude is fixed, then we can skip to the next amplitude
    // and avoid all of the symmetrization computation as well
    // as checking each individual factor for free parameters.
    if( a.m_termsValid && m_vbIsAmpFixed[iAmpIndex] ) continue;
    
    // now figure out if we need to recalculate a factors for an amplitude
    // we may not have changed a parameter since the last call to
    // calcTerms --
    // a parameter in the last iteration that is related to a factor in
    // this amplitude

    bool recalculateFactors = false;

    const Amplitude* pCurrAmp = 0;
    
    for( iFactor=0; iFactor < iNFactors; iFactor++ ){
      
      pCurrAmp = vAmps.at( iFactor );
      
      // now check to see if the value of this amplitudes parameters
      // are the same as they were the last time they were evaluated
      // for this particular dataset -- if not, recalculate

      if( !( a.m_termsValid && m_optimizeParIteration &&
            m_dataAmpIteration[&a][pCurrAmp] == m_ampIteration[pCurrAmp] ) ){

        recalculateFactors = true;
      }
    }
    
    if( !recalculateFactors ) continue;
    
    // if we get to here, we are changing the stored factors of the
    // amplitude
    
    modifiedTerm = true;
    
    // calculate all the factors that make up an amplitude for
    // for all events serially on CPU or in parallel on GPU
    unsigned long long iLocalOffset = 0;
    for( iFactor=0; iFactor < iNFactors;
         iFactor++, iLocalOffset += 2 * a.m_iNEvents * iNPermutations ){
      
      pCurrAmp = vAmps.at( iFactor );
      
      // if we have static user data, look up the location in the data array
      // if not, then look up by identifier
      unsigned long long uOffset =
      ( pCurrAmp->areUserVarsStatic() ?
        a.m_userVarsOffset[pCurrAmp->name()] :
        a.m_userVarsOffset[pCurrAmp->identifier()] );

#ifndef GPU_ACCELERATION
      pCurrAmp->
      calcAmplitudeAll( a.m_pdData,
                        a.m_pdAmpFactors + iLocalOffset,
                        a.m_iNEvents, &vvPermuations,
                        a.m_pdUserVars + uOffset );
#else
      a.m_gpuMan.calcAmplitudeAll( pCurrAmp, iLocalOffset,
                                   &vvPermuations,
                                   uOffset );
#endif//GPU_ACCELERATION
    }
    
    
    // now assemble all the factors in an amplitude into a single
    // symmetrized amplitude for each event
    
#ifndef GPU_ACCELERATION
    
    GDouble dSymmFactor = 1.0f/sqrt( iNPermutations );
    GDouble dAmpFacRe, dAmpFacIm, dTRe, dTIm;
    int iEvent, iPerm;
    unsigned long long iOffsetA, iOffsetP, iOffsetF;
    
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
        iOffsetP = 2 * a.m_iNEvents * iPerm + 2 * iEvent;
        
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
    
#else
    // on the GPU the terms are assembled and never copied out
    // of GPU memory
    
    a.m_gpuMan.assembleTerms( iAmpIndex, iNFactors, iNPermutations );
#endif
    
  }
  
  a.m_termsValid = true;

  // finally make a loop and record the iteration of parameters
  // that was used to make this calculation -- this cannot
  // be done above because some Amplitude instances may
  // be reused for a number of amplitudes and it is important
  // to recalculate all amplitudes that depend on a changed
  // parameter
  for( iAmpIndex = 0; iAmpIndex < iNAmps; iAmpIndex++ )
  {
    vector< const Amplitude* > vAmps =
    m_mapNameToAmps.find(ampNames.at(iAmpIndex))->second;
    
    int iFactor, iNFactors = vAmps.size();
    
    for( iFactor=0; iFactor < iNFactors; iFactor++ ){

      const Amplitude* pCurrAmp = vAmps.at( iFactor );
      m_dataAmpIteration[&a][pCurrAmp] = m_ampIteration[pCurrAmp];
    }
  }

  return modifiedTerm;
}

double
AmplitudeManager::calcIntensities( AmpVecs& a ) const
{

#ifdef VTRACE
  VT_TRACER( "AmplitudeManager::calcIntensities" );
#endif

  // check to be sure destination memory has been allocated
  assert( a.m_pdIntensity );
  
  double maxInten = 0;
  
  // first update the amplitudes
  calcTerms( a );
  
  // In GPU running mode amplitudes are maintained on the GPU and
  // the sum of the logs of intensities are calculated directly.
  // This is a CPU calculation that was likely called from the
  // parent IntensityManager with a single Kinematics object or
  // perhaps during MC generation. Copy the amplitudes out of the
  // GPU and compute the intensities (on the CPU).  There is no
  // GPU accelerated intensity calculation, just a GPU accelerated
  // log( intensity ) calculation.
  
#ifdef GPU_ACCELERATION
  if( a.m_pdAmps == NULL ) a.allocateCPUAmpStorage( *this );
  a.m_gpuMan.copyAmpsFromGPU( a );
#endif

  const vector< string >& ampNames = getTermNames();
  
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
AmplitudeManager::calcSumLogIntensity( AmpVecs& a ) const
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
  
  calcIntensities( a );
  
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

  if( !a.m_termsValid || hasTermWithFreeParam() ){

    calcTerms( a );
  }
  
  dSumLogI = a.m_gpuMan.calcSumLogIntensity( gpuProdPars, m_sumCoherently );
  
#endif
  
  return( dSumLogI );
}


void
AmplitudeManager::calcIntegrals( AmpVecs& a, int iNGenEvents ) const
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
  bool termChanged = calcTerms( a );
  
  // if nothing changed and it isn't the first pass, return
  if( !termChanged && a.m_integralValid ) return;
    
  int iNAmps = a.m_iNTerms;
  
  int i, j, iEvent;
  for( i = 0; i < iNAmps;i++ ) {
    for( j = 0; j <= i; j++ ) {
     
      // if the amplitude isn't floating and it isn't the first pass
      // through these data, then its integral can't change
      if( a.m_integralValid && m_vbIsAmpFixed[i] && m_vbIsAmpFixed[j] ){
        
        // if the amplitude isn't floating and it isn't the first pass
        // through these data, then its integral can't change
        
        continue;
      }
      else{
        
        // otherwise zero it out and recalculate it
        
        integralMatrix[2*i*iNAmps+2*j] = 0;
        integralMatrix[2*i*iNAmps+2*j+1] = 0;
      }
      
      // if two amps don't interfere the relevant integral is zero
      if( m_sumCoherently[i][j] ){
	        
#ifndef GPU_ACCELERATION
	
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
        // normalize
        integralMatrix[2*i*iNAmps+2*j] /= static_cast< GDouble >( iNGenEvents );
        integralMatrix[2*i*iNAmps+2*j+1] /= static_cast< GDouble >( iNGenEvents );	

#else
        a.m_gpuMan.calcIntegral( &(integralMatrix[2*i*iNAmps+2*j]), i, j, iNGenEvents );
#endif
      }
      
      // complex conjugate
      if( i != j ) {
        
        integralMatrix[2*j*iNAmps+2*i] = integralMatrix[2*i*iNAmps+2*j];
        integralMatrix[2*j*iNAmps+2*i+1] = -integralMatrix[2*i*iNAmps+2*j+1];
      }
    }
  }
  
  a.m_integralValid = true;
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
  
  m_needsUserVarsOnly = m_needsUserVarsOnly && newAmp->needsUserVarsOnly();
  
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

