#if !(defined DATAREADER)
#define DATAREADER

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

using namespace std;

class Kinematics;

/**
 * This base class provides an interface to data.
 *
 * It is intended that the user writes a functional data reader that
 * inherits indirectly from this class (through the UserDataReader class)
 * and defines the necessary virtual member functions for accessing data.
 *
 * \ingroup IUAmpTools
 */

class DataReader
{
  
public:
  
  /**
   * This is the default constructor.
   */
  DataReader( ) : 
  m_isDefault(true) { }
  
  /**
   * This constructor takes a list of arguments and stores them.
   */
  DataReader( const vector< string >& args ) :
  m_isDefault(false),
  m_args(args) { }
  
  /**
   * This is the destructor.
   */
  virtual ~DataReader() {}
  
  /**
   * The user should override this function with one that provides
   * a pointer to a Kinematics object.  When the end of a source is
   * reached, a null pointer should be returned.
   *
   * Return of a pointer is done to save memory copy and allocation time.  
   * The class the receives the pointer is responsible for deleting memory
   * when it is finished with it.
   *
   * This method should be virtual in the user's class if the DataReaderMPI 
   * template is to be used.
   *
   * \see Kinematics
   * \see DataReaderMPI
   */
  virtual Kinematics* getEvent() = 0;
  
  /**
   * The user should override this function with one that resets the source
   * so that the next call to getEvent() begins reading from the first 
   * event in the source.
   *
   * This method should be virtual in the user's class if the DataReaderMPI 
   * template is to be used.
   *
   * \see DataReaderMPI
   */
  virtual void resetSource() = 0; 
  
  /**
   * The user should override this function with one that returns the number
   * of events in the source.  This is used for allocating memory to store
   * the data and amplitudes.
   *
   * This method should be virtual in the user's class if the DataReaderMPI 
   * template is to be used.
   *
   * \see DataReaderMPI
   */
  virtual unsigned int numEvents() const = 0;
  
  /**
   * The user should override this function with one that returns the 
   * class name of the derived data reader.
   */
  virtual string name() const = 0;
  
  /**
   * This method is overridden by the UserDataReader class and does not
   * need to be defined (or used) by the user.
   *
   * \see UserDataReader
   */
  virtual DataReader* newDataReader( const vector< string >& args ) const = 0;
  
  /**
   * This method is overridden by the UserDataReader class and does not
   * need to be defined (or used) by the user.
   *
   * \see UserDataReader
   */
  virtual DataReader* clone() const = 0;
  
  /**
   * Returns the list of arguments that was passed to the constructor.
   */
  virtual vector<string> arguments() const { return m_args; }
  
  /**
   * Returns true if this instance was created using the default constructor
   * and returns false otherwise.
   */
  virtual bool isDefault() const { return ( m_isDefault == true ); }
  
  
private:
  
  bool m_isDefault;
  
  vector<string> m_args;
  
  
};

#endif
