#if !defined(MINUITINTERFACE_MINUITPARAMETER_H)
#define MINUITINTERFACE_MINUITPARAMETER_H

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

#include "MinuitInterface/Parameter.h"
#include "MinuitInterface/AsymmetricError.h"

class MinuitParameterManager;

class MinuitParameter : public Parameter {
  
  friend class MinuitParameterManager;
  enum { kUnknownMinuitId = 0 };
  
public:
  MinuitParameter( const std::string& name,
                  MinuitParameterManager& aManager,
                  double initialValue = 0,
                  bool bounded = false,
                  double lowerBound = 0,
                  double upperBound = 0 ); 
  
  virtual ~MinuitParameter();
  
  // user allowed changes in the status of the parameter:
  void fix();   // keep parameter constant during minimization
  void free(); // allow parameter to float during minimization
  void bound( double lowerBound, double upperBound ); // limit parameter range in minimization
  void unbound();
  
  unsigned int minuitId() const;
  bool floating() const;
  
  double error() const; // symmetric (parabolic) error
  const AsymmetricError& asymmetricErrors() const;
  double globalCorrelationCoefficient() const;
  
  bool bounded() const;
  double lowerBound() const;
  double upperBound() const;
  
  void registerWithManager();
  void unregister();
  
  std::ostream& dump( std::ostream& aStream ) const;
  
protected:
  // status changes managed by MinuitParameterManager
  void invalidateErrors();
  void validateErrors();
  void setParabolicError( double newParabolicError );
  void setAsymmetricErrors( const std::pair<double,double>& newAsymmetricErrors );
  void setGlobalCorrelation( double globalCorrelationCoefficient );
  void setMinuitId( unsigned int minuitId );
  
private:
  MinuitParameter(); // stop default
  MinuitParameter( const MinuitParameter& ); // stop default
  const MinuitParameter& operator=( const MinuitParameter& rhs ); // stop default
  
  // --------- member data ----------
  
  unsigned int m_minuitId;
  bool m_floating;
  
  double m_lowerBound;
  double m_upperBound;
  bool m_bounded;
  
  bool m_validErrors;
  double m_parabolicError;
  AsymmetricError m_asymmetricErrors;
  
  double m_globalCorrelationCoefficient;
  
  MinuitParameterManager& m_parameterManager;
};
inline  std::ostream& operator<<( std::ostream& aStream, 
                                 const MinuitParameter& aParameter ) {
  return aParameter.dump( aStream );
}
#endif
