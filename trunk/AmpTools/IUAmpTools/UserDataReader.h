#if !(defined USERDATAREADER)
#define USERDATAREADER

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
#include "IUAmpTools/DataReader.h"

using namespace std;

/**
 * This class handles the creation and cloning of the user's data reader.
 *
 * It is intended that the user writes a functional data reader that
 * inherits from this class (and thus the DataReader base class through this 
 * class) and also defines the necessary virtual member functions for 
 * accessing data.
 *
 * \ingroup IUAmpTools
 */

template< class T >
class UserDataReader : public DataReader
{

public:

  /**
   * This is the default constructor.  It should be called in the default
   * constructor of the user's derived class.
   */
        UserDataReader< T >() : DataReader() { }


  /**
   * This constructor takes a list of arguments and stores them.  There should
   * be a corresponding constructor in the user's derived class that calls
   * this constructor.
   */
        UserDataReader< T >( const vector< string >& args ) : DataReader( args ) { }


  /**
   * This is the destructor.
   */
        virtual ~UserDataReader< T >() { }


  /**
   * This method can create a new data reader (of the derived type).
   */
        virtual DataReader* newDataReader( const vector< string >& args ) const{
          return new T( args );
        }


  /**
   * This method can create a clone of a data reader (of the derived type).
   */
        virtual DataReader* clone() const{
          return ( isDefault() ? new T() : new T( arguments() ) );
        }

};

#endif
