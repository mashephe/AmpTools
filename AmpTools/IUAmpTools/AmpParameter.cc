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

#include <cassert>
#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "IUAmpTools/AmpParameter.h"

AmpParameter::AmpParameter( const string& arg ) :
m_valPtr( &m_defaultValue ),
m_defaultValue( 1E9 ),
m_name("-"),
m_hasExternalPtr( false )
{
    int length = arg.size();
    assert( length >= 1 );
    
    if( strncmp( &(arg[0]), "[", 1 ) == 0 ){
        
        if( strncmp( &(arg[length-1]), "]", 1 ) != 0 ){
            
            cerr << "ERROR in AmpParameter argument: " << arg << endl;
            assert( false );
        }
        
        int nameLength = length - 2;
        m_name = arg.substr( 1, nameLength );
    }
    else{
        
        m_defaultValue = atof( arg.c_str() );
    }
}

AmpParameter::AmpParameter( const AmpParameter& ampPar ) : 
m_defaultValue( ampPar ),
m_name( ampPar.name() ) {

  if( ampPar.hasExternalPtr() ){
    
    m_valPtr = ampPar.valPtr();
    m_hasExternalPtr = true;
  }
  else{
    
    m_valPtr = &m_defaultValue;
    m_hasExternalPtr = false;
  }
}


AmpParameter&
AmpParameter::operator=( const AmpParameter& ampPar ){
  m_defaultValue = ampPar;
  m_name = ampPar.name();

  if( ampPar.hasExternalPtr() ){

    m_valPtr = ampPar.valPtr();
    m_hasExternalPtr = true;
  }
  else{
    
    m_valPtr = &m_defaultValue;
    m_hasExternalPtr = false;
  }

  return *this;
}


bool
AmpParameter::operator==( const AmpParameter& otherPar ) const {
  
  if( otherPar.name() != m_name ) return false;
  
  if( hasExternalPtr() && otherPar.hasExternalPtr() ){
    
    return( m_valPtr == otherPar.valPtr() );
  }
  
  if( !hasExternalPtr() && !otherPar.hasExternalPtr() ){
    
    return( (*m_valPtr) == otherPar );
  }
  
  // if we get here one parameter has an external pointer and the
  // other one doesn't -- return false
  
  return false;
}

void 
AmpParameter::setExternalValue( const double* ptr ){
  
  m_valPtr = ptr;
  m_hasExternalPtr = true;
}

void
AmpParameter::setValue( double val ){
  
  m_defaultValue = val;
  m_valPtr = &m_defaultValue;
  m_hasExternalPtr = false;
}
