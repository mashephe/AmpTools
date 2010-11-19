
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

#include <ios>
#include <iomanip>

#include "MinuitInterface/MinuitParameter.h"
#include "MinuitInterface/MinuitParameterManager.h"

using namespace std;

MinuitParameter::MinuitParameter(  const string& name,
                                   MinuitParameterManager& aManager,
                                   double initialValue,
				   bool bounded,
				   double lowerBound,
				   double upperBound ) :
   Parameter( name, initialValue, 0 ),
   m_minuitId( kUnknownMinuitId ),
   m_floating( true ),
   m_lowerBound( lowerBound ),
   m_upperBound( upperBound ),
   m_bounded( bounded ),
   m_validErrors( false ),
   m_asymmetricErrors(),
   m_parameterManager( aManager )
{
      registerWithManager();
}

MinuitParameter::~MinuitParameter()
{}

void 
MinuitParameter::registerWithManager() {
   m_parameterManager.registerParameter( this );
}

void MinuitParameter::unregister() {
   m_parameterManager.removeParameter( this );
}

void
MinuitParameter::fix() {
   m_floating = false;
}

void
MinuitParameter::free() {
   m_floating = true;
}

void
MinuitParameter::bound( double lowerBound, double upperBound ) {
   m_lowerBound = lowerBound;
   m_upperBound = upperBound;
   m_bounded = true;
}

void
MinuitParameter::unbound() {
   m_bounded = false;
   m_lowerBound = 0;
   m_upperBound = 0;
}

unsigned int 
MinuitParameter::minuitId() const {
   return m_minuitId;
}

bool
MinuitParameter::floating() const {
   return m_floating;
}

double
MinuitParameter::error() const {
   if ( !m_validErrors ) { m_parameterManager.updateErrors(); }
   return Parameter::error();
}

const AsymmetricError& 
MinuitParameter::asymmetricErrors() const {
   if ( !m_validErrors ) { m_parameterManager.updateErrors(); }
   return m_asymmetricErrors;
}

double 
MinuitParameter::globalCorrelationCoefficient() const {
   return m_globalCorrelationCoefficient;
}

bool
MinuitParameter::bounded() const {
   return m_bounded;
}

double 
MinuitParameter::lowerBound() const {
   return m_lowerBound;
}

double 
MinuitParameter::upperBound() const {
   return m_upperBound;
}

void
MinuitParameter::invalidateErrors() {
   m_validErrors = false;
}

void
MinuitParameter::validateErrors() {
   m_validErrors = true;
}

void
MinuitParameter::setAsymmetricErrors( const pair<double,double>& newAsymmetricErrors ) {
   m_asymmetricErrors = newAsymmetricErrors;
}

void
MinuitParameter::setGlobalCorrelation( double globalCorrelationCoefficient ) {
   m_globalCorrelationCoefficient = globalCorrelationCoefficient;
}

void
MinuitParameter::setMinuitId( unsigned int anId ) {
   m_minuitId = anId;
}

ostream&
MinuitParameter::dump( ostream& aStream ) const {
   const AsymmetricError asymmErr = asymmetricErrors();
   
   aStream << name();
   aStream.setf(ios::floatfield);
   aStream.precision(5);
   aStream << " " << setw(9) << value();
   aStream << " " << setw(9) << error();
   aStream << "   " << setw(9) << asymmErr.lower() << " " << setw(9) << asymmErr.upper();
   if ( bounded() ) {
      aStream << "   bounds: " << setw(9) << lowerBound() << " " << setw(9) << upperBound();
   }
   
   return aStream;
}
