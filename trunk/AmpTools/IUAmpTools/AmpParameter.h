#if !defined(AMPPARMETER)
#define AMPPARMETER

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

#include<string>
#include<iostream>

using namespace std;

/**
 * This class is used to track parameters that may be floating in Amplitudes.
 * Each parameter must be a floating point number and have a name associated
 * with it.  The class defines a standard conversion to type double and
 * allow the value for the paramter be defined by some external memory.
 *
 * \ingroup IUAmpTools
 */

class AmpParameter {
  
public:
  
  /**
   * The constructor for the AmpParameter.  This takes an argument that is
   * constructed by the ConfigFileParser and it parses that string argument
   * according to the following rules.  If the first character of arg is
   * "[" and the last character is "]" the substring inside of the square
   * braces is interpreted as the name of the parameter and sets the default
   * value to an arbitrary number (1E9 currently).  If the string
   * is not started and ended with braces it converts the string to a double
   * and sets the converted value to the default value of the parameter with
   * the name being the empty string.
   *
   * \param[in] arg string argument that is interpreted as described above
   *
   * \see setValue
   */
  AmpParameter( const string& arg );
  
  /**
   * The default constructor.
   */
  AmpParameter(){}
  
  /**
   * A copy constructor.  If the AmpParameter being copied points to external
   * memory, the copy will also point to the same external memory.
   */
  AmpParameter( const AmpParameter& ampPar );

  /**
   * The assignment operator (similar to the copy constructor).
   */
  AmpParameter& operator=( const AmpParameter& ampPar );

  /**
   * A function to return the string name of this parameter.
   */
  string name() const { return m_name; }
  
  /**
   * An equivalence operator.  For parameters with external values, this
   * compares the values of the pointers.
   */
  bool operator==( const AmpParameter& otherPar ) const;

  /**
   * An operator for interpreting the class as a double.
   */
  inline operator double() const { return *m_valPtr; }
  
  /**
   * A function that returns a true if the AmpParameter is currently obtaining
   * its value from some external source.
   */
  bool hasExternalPtr() const { return m_hasExternalPtr; }

  /**
   * A function that returns a const pointer to the value of the parameter.
   */
  const double* valPtr() const { return m_valPtr; }
  
  /**
   * A function to tell the AmpParmaeter to find its value from some external
   * location in memory.  This is useful for fitting where the ParameterManager
   * maintains this external memory.
   *
   * \param[in] ptr the memory location of the intended value of the parameter
   *
   * \see Amplitude::setParPtr
   * \see AmplitudeManager::setAmpParPtr
   */
  void setExternalValue( const double* ptr );

  /**
   * A function to set the value of the parameter.
   *
   * \param[in] val the desired value of the parameter
   *
   * \see Amplitude::setParValue
   * \see AmplitudeManager::setAmpParValue
   */
  void setValue( double val );
  
  /**
   * A function to change the name of the parmeter.
   *
   * \param[in] name the new name of the parameter
   */
  void setName( const string& name ) { m_name = name; }
  
private:
  
  const double* m_valPtr;
  double m_defaultValue;
  string m_name;
  
  bool m_hasExternalPtr;
};

#endif

