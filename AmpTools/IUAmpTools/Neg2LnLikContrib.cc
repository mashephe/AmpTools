#include <map>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

#include "IUAmpTools/Neg2LnLikContrib.h"

MinuitMinimizationManager* Neg2LnLikContrib::m_minManager;

double Neg2LnLikContrib::operator()(){

  return neg2LnLikelihood();
}

double Neg2LnLikContrib::neg2LnLikelihood(){

  return 0.;
}

string
Neg2LnLikContrib::identifier() const {
  
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

bool
Neg2LnLikContrib::setParPtr( const string& name, const double* ptr ) const {
  
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
Neg2LnLikContrib::setParValue( const string& name, double val ) const {
  
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
Neg2LnLikContrib::updatePar( const string& name ) const {
  
  bool foundPar = false;
  
  for( vector< AmpParameter* >::const_iterator parItr = m_registeredParams.begin();
      parItr != m_registeredParams.end();
      ++parItr ){
    
    if( (**parItr).name().compare( name ) == 0 ){
      
      // The const_cast is a little bit undesirable here.  It can be removed
      // at the expensive of requiring the user to declare all member data in
      // the class that is updated on a parameter update "mutable."
      // Since we are trying to maximize user-friendliness, for now we will
      // remove this potential annoyance.
      
      const_cast< Neg2LnLikContrib* >(this)->updatePar( **parItr );
      foundPar = true;
    }
  }
  
  return foundPar;
}

void
Neg2LnLikContrib::registerParameter( AmpParameter& par ){
  
  m_registeredParams.push_back( &par );
}
