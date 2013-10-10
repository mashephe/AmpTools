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

#include <map>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/Kinematics.h"

void 
Amplitude::calcAmplitudeAll( GDouble* pdData, GDouble* pdAmps, int iNEvents,
                            const vector< vector< int > >* pvPermutations ) const
{
  complex< GDouble > cRes;
  
  int iPermutation, iNPermutations = pvPermutations->size();
  assert( iNPermutations );
  
  int iNParticles = pvPermutations->at(0).size();
  assert( iNParticles );
  
  GDouble** pKin = new GDouble*[iNParticles];
  
  /*
   if( m_registeredParams.size() != 0 ) {
   
   cout << "Current values of parameters: " << endl;
   for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
   parItr != m_registeredParams.end();
   ++parItr ){
   
   cout << "\t" << (**parItr).name() << ":  " << (**parItr) << endl;
   }
   }
   */
  
  int i, iEvent;
  for( iEvent=0; iEvent<iNEvents; iEvent++ ){        
    
    for( iPermutation = 0; iPermutation < iNPermutations; iPermutation++ ){
      
      m_currentPermutation = (*pvPermutations)[iPermutation];
      
      for( i = 0; i < iNParticles; i++ ){
        
        int j = (*pvPermutations)[iPermutation][i];
        pKin[i] = &(pdData[4*iNParticles*iEvent+4*j]);
        
      }
      
      // pKin is an array of pointers to the particle four-momentum
      // that gets reordered for each permutation so the user
      // doesn't need to deal with permutations in their calcAmplitude
      // routine
      
      cRes = calcAmplitude( pKin );
      
      pdAmps[2*iNEvents*iPermutation+2*iEvent] = cRes.real();
      pdAmps[2*iNEvents*iPermutation+2*iEvent+1] = cRes.imag();
    }
  }
  
  delete[] pKin;
}


complex< GDouble >
Amplitude::calcAmplitude( const Kinematics* pKin ) const {
  
  vector<int> permutation;
  
  vector<HepLorentzVector> particleList = pKin->particleList();
  
  for (int i = 0; i < particleList.size(); i++){
    permutation.push_back(i);
  }
  
  return calcAmplitude( pKin, permutation );
  
}


complex< GDouble >
Amplitude::calcAmplitude( const Kinematics* pKin, const vector< int >& permutation) const {
  
  vector<HepLorentzVector> particleList = pKin->particleList();
  
  GDouble** pData = new GDouble*[particleList.size()];
  
  if (particleList.size() != permutation.size()) assert(false);
  
  for (int i = 0; i < particleList.size(); i++){
    pData[i] = new GDouble[4];
    pData[i][0] = particleList[permutation[i]].e();
    pData[i][1] = particleList[permutation[i]].px();
    pData[i][2] = particleList[permutation[i]].py();
    pData[i][3] = particleList[permutation[i]].pz();
  }
  
  complex< GDouble > value = calcAmplitude(pData);
  
  for (int i = 0; i < particleList.size(); i++){
    delete[] pData[i];
  }
  delete[] pData;
  
  return value;
  
}


#ifdef GPU_ACCELERATION 
void
Amplitude::calcAmplitudeGPU( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                            const vector< int >& perm ) const {
  
  m_currentPermutation = perm;
  launchGPUKernel( dimGrid, dimBlock, GPU_AMP_ARGS );
}
#endif


bool
Amplitude::containsFreeParameters() const {
  
  bool hasFreeParam = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).hasExternalPtr() ) hasFreeParam = true;
  }
  
  return hasFreeParam;
}

bool
Amplitude::setParPtr( const string& name, const double* ptr ) const {
  
  bool foundPar = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).name().compare( name ) == 0 ){
      
      foundPar = true;
      (**parItr).setExternalValue( ptr );
      
      // pass in the name here to
      // use the const member function here so we only have one const-cast
      // that calls the non-const user function
      updatePar( (**parItr).name() );
    }
  }
  
  return foundPar;
}

bool
Amplitude::setParValue( const string& name, double val ) const {
  
  bool foundPar = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).name().compare( name ) == 0 ){
      
      foundPar = true;
      (**parItr).setValue( val );
      
      // pass in the name here to
      // use the const member function here so we only have one const-cast
      // that calls the non-const user function
      updatePar( (**parItr).name() );
    }
  }
  
  return foundPar;
}

bool
Amplitude::updatePar( const string& name ) const {
  
  bool foundPar = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).name().compare( name ) == 0 ){
      
      // The const_cast is a little bit undesirable here.  It can be removed
      // at the expensive of requiring the user to declare all member data in
      // the Amplitude class that is updated on a parameter update "mutable."
      // Since we are trying to maximize user-friendliness, for now we will
      // remove this potential annoyance.
      
      const_cast< Amplitude* >(this)->updatePar( **parItr );
      foundPar = true;
    }
  }
  
  return foundPar;
}

void
Amplitude::registerParameter( AmpParameter& par ){
  
  m_registeredParams.push_back( &par );
}

