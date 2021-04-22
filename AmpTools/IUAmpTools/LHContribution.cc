#include <map>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

#include "IUAmpTools/LHContribution.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

MinuitMinimizationManager* LHContribution::m_minManager;
BinnedData LHContribution::m_data; 

double LHContribution::operator()(){
  cout << "CALLING THE OPERATOR() in " << __FILE__ << " " << __LINE__ << endl;
  return neg2LnLikelihood();
  cout << "AFTER CALL TO OPERATOR() in " << __FILE__ << " " << __LINE__ << endl;
}

double LHContribution::neg2LnLikelihood(){
  return 0.;
}

string
LHContribution::identifier() const {
  
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

double LHContribution::calcLHContribution(double x) const {
  return 0.;
}

#ifdef GPU_ACCELERATION 
void
LHContribution::calcLHContributionGPU( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
#ifdef VTRACE
  string info = name();
  info += "::calcLHContributionGPU";
  VT_TRACER( info.c_str() );
#endif
  launchGPUKernel( dimGrid, dimBlock, GPU_AMP_ARGS );
}
#endif

bool
LHContribution::setParPtr( const string& name, const double* ptr ) const {
  
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
LHContribution::setParValue( const string& name, double val ) const {
  
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
LHContribution::updatePar( const string& name ) const {
  
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
      
      const_cast< LHContribution* >(this)->updatePar( **parItr );
      foundPar = true;
    }
  }
  
  return foundPar;
}

void
LHContribution::registerParameter( AmpParameter& par ){
  m_registeredParams.push_back( &par );
}