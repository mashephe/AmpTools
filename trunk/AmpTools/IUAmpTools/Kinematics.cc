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

#include <vector>
#include <cassert>
#include <iostream>

#include "IUAmpTools/Kinematics.h"
#include "TLorentzVector.h"

using namespace std;

const vector< TLorentzVector >&
Kinematics::particleList() const {
  
  return m_particleList;
}

const TLorentzVector&
Kinematics::particle( unsigned int index ) const {
  
  if( index >= m_particleList.size() ){
    
    cout << "ERROR:  Request for particle index " << index << 
    ", but only " << m_particleList.size() << " particles exist." << endl;
    
    assert( false );
  }
  
  return m_particleList[index];
}

void
Kinematics::setParticleList( const vector< TLorentzVector >& particleList ){
  
  assert( particleList.size() <= kMaxParticles );
  
  m_particleList = particleList;
}

int Kinematics::m_globalEventID = 0;
