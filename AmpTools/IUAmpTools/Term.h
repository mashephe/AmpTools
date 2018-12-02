#if !defined( TERM )
#define TERM

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

#include <vector>

#include "IUAmpTools/Term.h"

using namespace std;

class AmpParameter;

class Term {
  
public:

  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  Term( ) : m_isDefault(true) { }

  /**
   * This constructor takes a list of arguments to inititialize an
   * amplitude and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  Term( const vector< string >& args ) : m_isDefault(false),
  m_args(args){ }

  /**
   * The default destructor.
   */
  
  virtual ~Term() {}
  
  /**
   * Must be overridden by the user to provide the name of the amplitude.
   * This is necessary to connect ConfigurationInfo to specific class
   * instances when the AmplitudeManager is being setup.
   */
  virtual string name() const = 0;
  
  /**
   * Returns a boolean to indicate if this amplitude contains at least one
   * floating parameter.
   */
  bool containsFreeParameters() const;
  
  /**
   * A function that indicates if this amplitude was created with the
   * default constructor.
   */
  bool isDefault() const { return ( m_isDefault == true ); }
  
  /**
   * Returns the list of arguments that was passed to the constructor.
   */
  vector<string> arguments() const { return m_args; }
  
  /**
   * The user can override this to do specific one-time tasks that
   * before the fit begins, but after the parameters have been initialized.
   * In general these tasks should go into the this init routine
   * rather than the consructor.
   */
  virtual void init(){}
  
  /**
   * This tells the AmpParameter with the indicated name to obtain its
   * value from some external location in memory.  This functionality is
   * is useful when fitting as the ParameterManager can automatically
   * updated these external values.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this amplitude, false otherwise.  The function calls
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] ptr the location in memory to find the value of the AmpParameter
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParPtr
   */
  bool setParPtr( const string& name, const double* ptr ) const;
  
  /**
   * This tells the AmpParameter with the indicated name to set its
   * value to the specified value.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this amplitude, false otherwise. Function calls
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] val the value to set AmpParameter to
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParValue
   */
  bool setParValue( const string& name, double val ) const;
  
  /**
   * The user may override this function to recalculate member data in the
   * amplitude class whenever a parameter changes.  For example, the user's
   * amplitude calculation may have expensive integration or other functions
   * that only need to be computed whenever a parameter changes, rather than
   * being computed on the fly for every event.
   *
   * \param[in] par a const reference to the parameter that has been updated
   */
  virtual void updatePar( const AmpParameter& par ) {}
  
  /**
   * \overload
   *
   * A function used by the framework to signal that a parameter with a
   * particular name has changed.  This function searches the registered
   * parameters to see if this amplitude contains that parameter.  If so
   * it calls the virtual function updatePar with the AmpParameter.  The
   * function returns true of a parameter with the specified name is found.
   *
   * \param[in] name the name of the updated parameter
   */
  
  bool updatePar( const string& name ) const;
  
protected:
  
  /**
   * Any user-defined class derived from Amplitude that has a parameter
   * in the amplitude should register the parameter using this routine.  This
   * shoudl be done for each parameter in the constructor of the user's
   * Amplitude class.
   *
   * \param[in] par a reference to the AmpParmeter object
   */
  void registerParameter( AmpParameter& par );
  
private:
  
  bool m_isDefault;
  vector<string> m_args;
  
  vector< AmpParameter* > m_registeredParams;
};


#endif // defined TERM
