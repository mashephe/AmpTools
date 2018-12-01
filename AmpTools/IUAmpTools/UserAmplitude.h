#if !(defined USERAMPLITUDE)
#define USERAMPLITUDE

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

#include <string>
#include <vector>
#include "IUAmpTools/Amplitude.h"

using namespace std;

/**
 * This class handles the creation and cloning of the user's amplitudes.
 *
 * It is intended that the user writes a functional amplitude that
 * inherits from this class (and thus the Amplitude base class through this
 * class) and also defines the necessary virtual member functions.
 *
 * \ingroup IUAmpTools
 */

template< class T >
class UserAmplitude : public Amplitude
{
  
public:
  
  /**
   * This is the default constructor.  It should be called in the default
   * constructor of the user's derived class.
   */
  UserAmplitude< T >() : Amplitude(), m_dataCache( NULL ) { }
  
  
  /**
   * This constructor takes a list of arguments and stores them.  There should
   * be a corresponding constructor in the user's derived class that calls
   * this constructor.
   */
  UserAmplitude< T >( const vector< string >& args ) :
  Amplitude( args ), m_dataCache( NULL ) { }
  
  
  /**
   * This is the destructor.
   */
  virtual ~UserAmplitude< T >() { }
  
  
  /**
   * This method can create a new amplitude (of the derived type).
   */
  Amplitude* newAmplitude( const vector< string >& args ) const{
    return new T( args );
  }
  
  
  /**
   * This method can create a clone of an amplitude (of the derived type).
   */
  Amplitude* clone() const{
    return ( isDefault() ? new T() : new T( arguments() ) );
  }
  
#ifdef GPU_ACCELERATION
  
  // for amplitudes with no free parameters, we can use the user
  // data block to facilitate a generic GPU-based repeated
  // calculation of the amplitude -- in that case we should
  // intercept the following function calls
  
  void calcUserData( GDouble** pKin, GDouble* userData ) const;
  unsigned int numUserVars() const;
  bool needsUserDataOnly() const;
  
  // if the user defines a custom GPU execution kernel
  // for the amplitude then this function should be
  // overriden to return true
  
  virtual bool isGPUEnabled() const { return false; }
  
#endif // GPU_ACCELERATION
  
private:
  
  mutable GDouble* m_dataCache;
  
};


#ifdef GPU_ACCELERATION

template< class T >
void UserAmplitude<T>::calcUserData( GDouble** pKin, GDouble* userData ) const {
  
  if( containsFreeParameters() || isGPUEnabled() ){
    
    // just pass the call through and return
    dynamic_cast< const T* >(this)->calcUserData( pKin, userData );
    return;
  }
  
  // if the user amplitude does not contain free parameters, we
  // can use the user data block to define a default GPU accelerated
  // version of the amplitude by calculating the amplitude and
  // and storing the result as user data and reading it back
  // with a standardized GPU amplitude kernel
  
  unsigned int numVars = dynamic_cast< const T* >(this)->numUserVars();
  
  if( m_dataCache == NULL ){
    
    m_dataCache = new GDouble[numVars];
  }
  
  complex< GDouble > result;
  
  if( numVars > 0 ){
   
    dynamic_cast< const T* >(this)->calcUserData( pKin, m_dataCache );
    result = calcAmplitude( pKin, m_dataCache );
  }
  else{
    
    result = calcAmplitude( pKin );
  }

  cout << result << endl;
  
  // and store the result as user data in the framework
  
  userData[0] = result.real();
  userData[1] = result.imag();
}

template< class T >
unsigned int UserAmplitude<T>::numUserVars() const {
  
  return ( containsFreeParameters() || isGPUEnabled() ?
	   dynamic_cast< const T* >(this)->numUserVars() : 2 );
}

template< class T >
bool UserAmplitude<T>::needsUserDataOnly() const {
  
  return (containsFreeParameters() || isGPUEnabled() ?
          dynamic_cast< const T* >(this)->needsUserDataOnly() : true );
}

#endif // GPU_ACCELERATION


#endif
