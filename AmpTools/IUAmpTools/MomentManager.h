#if !defined( MOMENTMANAGER )
#define MOMENTMANAGER

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


#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <complex>

#include "IUAmpTools/IntensityCalculator.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/AmpVecs.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "GPUManager/GPUCustomTypes.h"

#ifdef GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif //GPU_ACCELERATION

class Kinematics;
class NormIntInterface;

using namespace std;

/**
 * A class for managing the construction of a sum of moments.
 *
 * This is currently a placeholder for future development that will
 * allow moment-based fitting.
 *
 * There should be one instance of the MomentManager for each reaction or
 * set of unique final state particles in the fit.
 *
 * \ingroup IUAmpTools
 */

class MomentManager : public IntensityCalculator
{
  
public:
  
  /** Constructor.
   * Constructs an MomentManager
   *
   * \param[in] finalState a vector of strings, one to identify each final state
   * argument.  Particles with identical string identifiers will be automatically
   * permuted when the intesnity is calculated.
   *
   * \param[in] reactionName an optional name for the reaction
   *
   * \see getPermutations
   * \see addAmpPermutation
   */
  MomentManager( const vector< string >& finalState,
                 const string& reactionName = "");
  
  /** Destructor.
   */
  ~MomentManager();
  
  
private:
 
  
  
};

#endif
