#if !defined(LIKELIHOODCALCULATOR)
#define LIKELIHOODCALCULATOR

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


#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/DataReader.h"

#include "MinuitInterface/MIFunctionContribution.h"

class MinuitMinimizationManager;
class MinuitParameter;

/**
 * This class calculates -2 ln( likelihood ) for the fit.
 *
 * This class utilizes the data, amplitude manager, and the interface
 * to the normalization integrals to compute -2 ln likelihood for the fit.
 * The result of this computation is added to any other other contributions
 * (e.g., from other LikelihoodCalculators) and the sum is minimized by 
 * the MinuitInterface by varying the parameters.
 *
 * This class reads and caches the data internally.  
 *
 * \ingroup IUAmpTools
 */

class LikelihoodCalculator : public MIFunctionContribution
{

 public:
	
  LikelihoodCalculator( const AmplitudeManager& ampManager,
                        const NormIntInterface& normInt,
                        DataReader& dataReader,
                        const ParameterManager& parManager );
  
  ~LikelihoodCalculator(){}
	
  // where the likelihood calculation gets "assembled"
  virtual double operator()();

  // this function will be called on a parameter change -- do expensive
  // computations of the likelihood here that go into assembly of the
  // final likelihood above.  Dividing task this way allows for better
  // implementation of parallel computations
  virtual void update( const MISubject* );
  
  // need to override contribution method in MIFunctionContribution to
  // allow LikelihoodCalculator class to handle update calls correctly
  virtual double contribution() { return operator()(); }
   
protected:
	
  // helper functions -- also useful for pulling parts of the
  // likelihood calculation
  double dataTerm();
  double normIntTerm();
    
private:
			
  const AmplitudeManager& m_ampManager;
  const NormIntInterface& m_normInt;
  DataReader& m_dataReader;

  bool m_functionEvaluated;

  // calculate this "expensive" number once when parameters are updated
  double m_sumLnI;

  // The flat array of kinematics and amplitudes 
  AmpVecs m_ampVecs;
		
};

#endif
