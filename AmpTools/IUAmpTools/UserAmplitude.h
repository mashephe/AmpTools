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
  UserAmplitude< T >() : Amplitude() { }
  
  
  /**
   * This constructor takes a list of arguments and stores them.  There should
   * be a corresponding constructor in the user's derived class that calls
   * this constructor.
   */
  UserAmplitude< T >( const vector< string >& args ) : Amplitude( args ) { }
  
  
  /**
   * This is the destructor.
   */
   ~UserAmplitude< T >() {
    
    for( typename map< string, T* >::iterator amp = m_ampInstances.begin();
         amp != m_ampInstances.end(); ++amp ){
      
       delete amp->second;
    }
  }
  
  
  /**
   * This method can create a new amplitude (of the derived type).
   * First check to see if an existing instance will function in exactly
   * the same way (has the same arguments).  If it will, use that,
   * otherwise make a new one and save it.  (Reusing amplitude pointers
   * allows for memory optimization in the IntensityManager, which
   * is why this funny-looking roundabout method is here.)
   */
  Amplitude* newAmplitude( const vector< string >& args ) const{
    
    // we need the identifier of the new amplitude first:
    T* newAmp = new T( args );
    
    string ident = newAmp->identifier();
    
    if( m_ampInstances.find( ident ) == m_ampInstances.end() ){
      
      m_ampInstances[ident] = newAmp;
    }
    else{
      
      // already have a functional instance, so delete this one
      delete newAmp;
    }
    
    return m_ampInstances[ident];
  }
  
  
  /**
   * This method can create a clone of an amplitude (of the derived type).
   */
  Amplitude* clone() const {
    
    return ( isDefault() ? new T() : new T( arguments() ) );
   }
  
private:
  
  mutable map< string, T* > m_ampInstances;
  
};

#endif
