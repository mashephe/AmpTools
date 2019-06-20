//******************************************************************************
// This file is part of AmpTools, a package for performing Term Analysis
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

#include <vector>
#include <cassert>

#include "IUAmpTools/Term.h"
#include "IUAmpTools/AmpParameter.h"

using namespace std;

void
Term::registerParameter( AmpParameter& par ){
  
  m_registeredParams.push_back( &par );
}


bool
Term::containsFreeParameters() const {
  
  bool hasFreeParam = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).hasExternalPtr() ) hasFreeParam = true;
  }
  
  return hasFreeParam;
}

bool
Term::setParPtr( const string& name, const double* ptr ) const {
  
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
Term::setParValue( const string& name, double val ) const {
  
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
Term::updatePar( const string& name ) const {
  
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
      // the Term class that is updated on a parameter update "mutable."
      // Since we are trying to maximize user-friendliness, for now we will
      // remove this potential annoyance.
      
      const_cast< Term* >(this)->updatePar( **parItr );
      foundPar = true;
    }
  }
  
  return foundPar;
}

string
Term::identifier() const {
  
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
Term::calcUserVarsAll( GDouble* pdData, GDouble* pdUserVars, int iNEvents,
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
  
  // if the user variables are static, add a pointer to the data set for
  // which we just did the calculation;  that we know later
  // whether it has been calculated
  if( areUserVarsStatic() ) m_staticUserVarsCalculated.insert( pdData );
  m_userVarsCalculated.insert( pdData );
  
  delete[] pKin;
}

set< GDouble* > Term::m_staticUserVarsCalculated = set< GDouble* >();

