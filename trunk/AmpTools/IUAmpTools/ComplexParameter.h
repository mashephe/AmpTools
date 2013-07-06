#if !defined(COMPLEXPARAMETER)
#define COMPLEXPARAMETER

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
#include <complex>

#include "MinuitInterface/MIObserver.h"
#include "MinuitInterface/MinuitParameter.h"

class MinuitMinimizationManager;
class MISubject;

using namespace std;

/**
 * This class manages the conversion of complex parameters
 * to real parameters that minuit minimizes.  
 * 
 * If a valid MinuitMinimizationManager reference is passed in it
 * it registers two parameters with MinuitParameterManager when created.
 * It is an also MIObserver -- on update calls it maintains a
 * complex< double > that represents the value of the parameter.
 *
 * It is also useful to be able to create ComplexParameter objects with
 * no connection to MinuitMinimizationManager in MPI implementations of
 * the fitter.  These objects are valid, but are not registered to or updated 
 * by the minimization manager.
 *
 * \ingroup IUAmpTools
 */

class ComplexParameter : public MIObserver
{
  
public:
  
  /**
   * This constructs the ComplexParameter.
   * 
   * \param[in] name the string giving the name of the parameter, real and
   * imaginary parts of the parameter will be identified in MINUIT with this
   * string plus "_re" or "_im" appened.
   *
   * \param[in] aManager an instance of the MinuitMinimizationManager with
   * which to create and register the parameters
   *
   * \param[in] initialValue an optional initial value
   *
   * \param[in] purelyReal an optional argument that when set true forces
   * the imaginary part of the complex number to zero
   */
  
  ComplexParameter( const string& name,
                   MinuitMinimizationManager& aManager,
                   complex< double > initialValue = complex< double >( 1, 0 ),
                   bool purelyReal = false );
  
  /**
   * This is the ComplexParameter destructor.
   */
  ~ComplexParameter();
  
  /**
   * This is called by the MinuitInterface and triggers the updating of the
   * values of the real and imaginary parts of the complex parameter.
   *
   * \param[in] callingSubject pointer to the MinuitParameter representing
   * either the real or imaginary part of the value
   */
  void update( const MISubject* callingSubject );
  
  /**
   * This sets the value of the parameter.
   *
   * \param[in] value the desired value of the complex parameter
   */
  void setValue( complex< double > value );
  
  /**
   * This fixes the parameter at its current value.
   *
   */
  void fix();
  
  /**
   * This allows the parameter to float.
   *
   */
  void free();
  
  /**
   * A function that returns the name of the ComplexParameter
   */
  const string& name() const { return m_name; }
  
  /**
   * A function that returns the name of the parameter representing the real
   * part of the complex number.
   */
  string realName() const { return m_name + "_re"; }
  
  /**
   * A function that returns the name of the parameter representing the imaginary
   * part of the complex number.
   */
  string imagName() const { return m_name + "_im"; }
  
  /**
   * A function that returns a boolean indicating whether the parameter is
   * purely real or not.
   */
  bool isPurelyReal() const { return m_purelyReal == true; }
  
  /**
   * A function to see if the parameter is floating.  
   */
  bool isFixed() const;
   
  /**
   * A function that returns the current value of the complex parameter.
   */
  complex< double > value() const { return m_value; }
  
  /**
   * A function that returns a const pointer to the value of the parameter.
   * This is useful, for example, when telling the AmplitudeManager to 
   * use an external parameter for a production amplitude.
   *
   * \see AmplitudeManager::setExternalProductionAmplitude
   */
  const complex< double >* constValuePtr() const { return &m_value; }
  
  /**
   * A function that returns a pointer (not const) to the value of the
   * parameter.
   */
  complex< double >* valuePtr() { return &m_value; }
  
private:
  
  // prevent default construction or copying
  ComplexParameter();
  ComplexParameter( const ComplexParameter& );
  
  string m_name;
  complex< double > m_value;
  bool m_purelyReal;
  
  // parameters to hold the two numbers that represent this complex parameter
  MinuitParameter* m_realPar;
  MinuitParameter* m_imPar;
};

#endif
