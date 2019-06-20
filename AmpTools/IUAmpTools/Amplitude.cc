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

#ifdef VTRACE
#include "vt_user.h"
#endif

string
Amplitude::identifier() const {
  
  string id = name();
  
  // make the identifer unique from a likely name in the case of no args
  id += "%%";
  
  for( vector< string >::const_iterator arg = m_args.begin();
      arg != m_args.end(); ++arg ){
    
    id += *arg;
    id += " ";
  }
  
  return id;
}

void
Amplitude::calcUserVarsAll( GDouble* pdData, GDouble* pdUserVars, int iNEvents,
                            const vector< vector< int > >* pvPermutations ) const
{
  
#ifdef VTRACE
  string info = name();
  info += "::calcUserVarsAll";
  VT_TRACER( info.c_str() );
#endif
  
//  cout << "Caculating user data for " << name() << endl;
  
  unsigned int numVars = numUserVars();
  
  // exit immediately if there is nothing to compute
  if( !numVars ) return;
  
  int iPermutation, iNPermutations = pvPermutations->size();
  assert( iNPermutations );
  
  int iNParticles = pvPermutations->at(0).size();
  assert( iNParticles );
  
  GDouble** pKin = new GDouble*[iNParticles];
  
  int i, iEvent;
  for( iEvent=0; iEvent<iNEvents; iEvent++ ){
    
    for( iPermutation = 0; iPermutation < iNPermutations; iPermutation++ ){
      
      m_currentPermutation = (*pvPermutations)[iPermutation];

      // pKin is an array of pointers to the particle four-momentum
      // that gets reordered for each permutation so the user
      // doesn't need to deal with permutations in their calcAmplitude
      // routine
      
      for( i = 0; i < iNParticles; i++ ){
        
        int j = (*pvPermutations)[iPermutation][i];
        pKin[i] = &(pdData[4*iNParticles*iEvent+4*j]);
        
      }
      
      unsigned int userIndex = iNEvents*iPermutation*numVars + iEvent*numVars;
      calcUserVars( pKin, &(pdUserVars[userIndex]) );
    }
  }
  
  delete[] pKin;
}

void
Amplitude::calcAmplitudeAll( GDouble* pdData, GDouble* pdAmps, int iNEvents,
                            const vector< vector< int > >* pvPermutations,
                             GDouble* pdUserVars ) const
{
  
#ifdef VTRACE
  string info = name();
  info += "::calcAmplitudeAll";
  VT_TRACER( info.c_str() );
#endif

  complex< GDouble > cRes;
  
  unsigned int numVars = numUserVars();

  int iPermutation, iNPermutations = pvPermutations->size();
  assert( iNPermutations );
  
  int iNParticles = pvPermutations->at(0).size();
  assert( iNParticles );
  
  GDouble** pKin = new GDouble*[iNParticles];
  
  int i, iEvent;
  for( iEvent = 0; iEvent < iNEvents; iEvent++ ){
    
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
      
      unsigned int userIndex = iNEvents*iPermutation*numVars + iEvent*numVars;

      if( numVars != 0 ){
      
        cRes = calcAmplitude( pKin, &(pdUserVars[userIndex]) );
      }
      else{
        cRes = calcAmplitude( pKin );
      }
      
      pdAmps[2*iNEvents*iPermutation+2*iEvent] = cRes.real();
      pdAmps[2*iNEvents*iPermutation+2*iEvent+1] = cRes.imag();
    }
  }
  
  delete[] pKin;
}


complex< GDouble >
Amplitude::calcAmplitude( const Kinematics* pKin, GDouble* userVars ) const {
  
  vector<int> permutation;
  
  vector<TLorentzVector> particleList = pKin->particleList();
  
  for (int i = 0; i < particleList.size(); i++){
    permutation.push_back(i);
  }
  
  return calcAmplitude( pKin, permutation, userVars );
  
}

complex< GDouble >
Amplitude::calcAmplitude( GDouble** pKin, GDouble* userVars ) const {
  
  cout
  << "***********************************************************\n"
  << "ERROR in the construction of the class that defines\n"
  << "the Amplitude named " << name() << ".\n"
  << "One of the following two cases result in this error.\n\n"
  << "(1) The numUserVars() method of the class indicates that\n"
  << "    at least one user-defined variable will be calculated,\n"
  << "    but the calcAmplitude method hasn't been defined such\n"
  << "    that it can accept a pointer to the user-defined data\n"
  << "    block.  Please define the function:\n"
  << "      " << name() << "::\n"
  << "         calcAmplitude( GDouble** pKin, GDouble* userVars )\n\n"
  << "(2) No calcAmplitude function (with or without user data\n"
  << "    data pointer is defined in the class.\n"
  << "***********************************************************\n"
  << endl;
  
  assert( false );
  
}

complex< GDouble >
Amplitude::calcAmplitude( GDouble** pKin ) const {

  // It is possible to end up here if the user has
  // defined calcAmplitude such that it takes two
  // arguments and the number of user variables
  // to calculate is zero. (This is the else clause
  // in the next method.)  In this case try to call
  // the user's calcAmplitude function by passing
  // in a NULL pointer to the user data block.
  // If that isn't defined either then the error
  // above will print and the program will exit.
  
  return calcAmplitude( pKin, NULL );
}



complex< GDouble >
Amplitude::calcAmplitude( const Kinematics* pKin,
                          const vector< int >& permutation,
                          GDouble* userVars ) const {

#ifdef VTRACE
  string info = name();
  info += "::calcAmplitude";
  VT_TRACER( info.c_str() );
#endif

  vector<TLorentzVector> particleList = pKin->particleList();
  
  GDouble** pData = new GDouble*[particleList.size()];
  
  if (particleList.size() != permutation.size()) assert(false);
  
  for (int i = 0; i < particleList.size(); i++){
    pData[i] = new GDouble[4];
    pData[i][0] = particleList[permutation[i]].E();
    pData[i][1] = particleList[permutation[i]].Px();
    pData[i][2] = particleList[permutation[i]].Py();
    pData[i][3] = particleList[permutation[i]].Pz();
  }
  
  complex< GDouble > value;
  
  if( userVars != NULL ){
  
    value = calcAmplitude( pData, userVars );
  }
  else{

    // this call should ensure backwards compatibility with
    // older Amplitude classes that just took one argument
    
    value = calcAmplitude( pData );
  }
  
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
#ifdef VTRACE
  string info = name();
  info += "::calcAmplitudeGPU";
  VT_TRACER( info.c_str() );
#endif

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
  
#ifdef VTRACE
  string info = (*this).name();
  info += "::updatePar [";
  info += name.c_str();
  info += "]";
  VT_TRACER( info.c_str() );
#endif

  
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
