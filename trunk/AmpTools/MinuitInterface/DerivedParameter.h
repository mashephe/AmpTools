#if !defined(MINUITINTERFACE_DERIVEDPARAMETER_H)
#define MINUITINTERFACE_DERIVEDPARAMETER_H

// This file is a part of MinuitInterface - a front end for the Minuit minimization
//       package (Minuit itself was authored by Fred James, of CERN)
// 
// 
// Copyright Cornell University 1993, 1996, All Rights Reserved.
// 
// This software written by Lawrence Gibbons, Cornell University.
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
// obtained from Cornell University.
// 
// CORNELL MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By way
// of example, but not limitation, CORNELL MAKES NO REPRESENTATIONS OR
// WARRANTIES OF MERCANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
// THE USE OF THIS SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS,
// COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS.  Cornell University shall not be
// held liable for any liability with respect to any claim by the user or any
// other party arising from use of the program.
//

// The DerivedParameter is a parameter that depends on one or
// more other parameters via a known functional form.  The
// template class T communicates the set of "parent"
// parameters upon which this class depends.  The only
// constraint on T is 
//   1) an "attach()" member function, which causes this 
//      parameter to be attached to the base MISubject
//      of all of the parent Parameters.  That will cause this
//      DerivedParameter to be updated whenever any of the 
//      parent parameters is updated.
//   2) similarly, a "detach()" member function
//
// The update is handled by a ConversionFunction, which implements
// an operator()( T& ) member function, with a pure virual definition
// in the ConversionFunction base class.

#include <string>

#include "MinuitInterface/MIObserver.h"
#include "MinuitInterface/Parameter.h"
#include "MinuitInterface/ConversionFunction.h"

template< class T>
class DerivedParameter : public Parameter, public MIObserver
{
public:
   
   DerivedParameter( const std::string& name, 
                     T& parentParameters, 
                     ConversionFunction<T>& converter ); 
   
   virtual ~DerivedParameter();
   
   void update( const MISubject* updatedParameter );
   
protected:
   // DerivedParameter only updates via changes in any of the base parameters
   // upon which it depends
      void setValue( double newValue ) { Parameter::setValue( newValue ); }
   void setError( double newError, bool notify = false ) { Parameter::setError( newError, notify ); }
   void setValueError( double newValue, double newError ) { Parameter::setValueError( newValue, newError ); }
   
private:
   // ------------ unused copy/assignment constructors
   DerivedParameter(); // stop default
   DerivedParameter( const DerivedParameter& ); // stop default

   // ------------ member data ----------
   ConversionFunction<T>& m_converter;
   T& m_baseParameters;
};

template<class T>
DerivedParameter<T>::DerivedParameter(  const std::string& name, 
                                        T& parentParameters,
                                        ConversionFunction<T>& converter ) :
Parameter(name),
MIObserver(),
m_converter( converter ),
m_baseParameters( parentParameters ) 
{
   m_baseParameters.attach( this );   
   setValue( m_converter( m_baseParameters ) );
}

template<class T>
DerivedParameter<T>::~DerivedParameter() {
   m_baseParameters.detach( this );
}

template<class T>
void 
DerivedParameter<T>::update( const MISubject* updatedParameter ) {
   setValueError( m_converter(m_baseParameters), m_converter.error(m_baseParameters) );
}
#endif
