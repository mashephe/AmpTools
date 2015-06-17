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
#include <string>
#include <complex>
#include <cassert>

#include "IUAmpTools/ComplexParameter.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MinuitParameterManager.h"

ComplexParameter::ComplexParameter( const string& name,
                                   MinuitMinimizationManager& fitManager,
                                   complex< double > initialValue,
                                   bool purelyReal ) :
MIObserver(),
m_name( name ),
m_value( initialValue ),
m_purelyReal( purelyReal ),
m_realPar( NULL ),
m_imPar( NULL )
{ 
    
  // be sure we have a valid reference to MinuitMinimizationManager before
  // registering the new parameters -- otherwise this object looks
  // just like complex number
  
  m_realPar = new MinuitParameter( name + "_re", 
                                  fitManager.parameterManager(),
                                  real( m_value ) );
  
  m_realPar->attach( this );
  
  if( !purelyReal ){
    
    m_imPar = new MinuitParameter( name + "_im", 
                                  fitManager.parameterManager(),
                                  imag( m_value ) );
    m_imPar->attach( this );
  }
  
}

ComplexParameter::~ComplexParameter()
{
  //  not yet supported un UpRootMinuit
  // m_realMagPar->unregister();
  // m_imPhasePar->unregister();
  
  // problems with double delete?
  // delete m_realPar;
  // delete m_imPar;
}

void 
ComplexParameter::update( const MISubject* callingSubject ){
  
  if( callingSubject == m_realPar ){
    
    m_value = complex< double >( m_realPar->value(), imag( m_value ) );
  }
  else if( callingSubject == m_imPar ){
    
    m_value = complex< double >( real( m_value ), m_imPar->value() );
  }
}

void
ComplexParameter::setValue( complex< double > value ){
  
  // the set value calls below will trigger the notify() method
  // of MISubject, which will call the update method above
  // to set the member data m_value to the new value
  
  m_realPar->setValue( real( value ) );
  if( !m_purelyReal ) m_imPar->setValue( imag( value ) );
}

void
ComplexParameter::fix(){
  
  m_realPar->fix();
  if( !m_purelyReal ) m_imPar->fix();
}

void
ComplexParameter::free(){
  
  m_realPar->free();
  if( !m_purelyReal ) m_imPar->free();
}

bool
ComplexParameter::isFixed() const {

  return !m_realPar->floating();
}

