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

#include <string>

#include "IUAmpTools/AmpVecs.h"
#include "IUAmpTools/IntensityManager.h"

#include "MinuitInterface/MIFunctionContribution.h"

class NormIntInterface;
class ParameterManager;
class DataReader;
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

using namespace std;

class LikelihoodCalculator : public MIFunctionContribution
{
  
public:
  
  LikelihoodCalculator( const IntensityManager& intenManager,
                        const NormIntInterface& normInt,
                        DataReader* dataReaderSignal,
                        DataReader* dataReaderBkgnd,
                        const ParameterManager& parManager );
  
  virtual ~LikelihoodCalculator();
  
  string reactionName() const { return m_intenManager.reactionName(); }
  
  // this method delivers the likelihood
  double operator()();
  
protected:
  
  // helper functions -- also useful for pulling parts of the
  // likelihood calculation in parallel implementations
  double dataTerm();
  double normIntTerm();
  
  // these are useful for MPI implementations since there are a
  // few sums that must be maintained across all processes to properly
  // compute the normalization integral terms of the likelihood
  
  double sumBkgWeights() const { return m_sumBkgWeights; }
  double numBkgEvents()  const { return m_numBkgEvents;  }
  double numDataEvents() const { return m_numDataEvents; }
  
  void setSumBkgWeights( double sum ) { m_sumBkgWeights = sum; }
  void setNumBkgEvents ( double num ) { m_numBkgEvents  = num; }
  void setNumDataEvents( double num ) { m_numDataEvents = num; }
  
private:
  
  bool m_hasBackground;
  
  const IntensityManager& m_intenManager;
  const NormIntInterface& m_normInt;

  DataReader* m_dataReaderSignal;
  DataReader* m_dataReaderBkgnd;
  
  bool m_firstNormIntCalc;
  bool m_firstDataCalc;
  
  double* m_prodFactorArray;
  const double* m_normIntArray;
  const double* m_ampIntArray;
  
  // The flat array of kinematics and amplitudes 
  AmpVecs m_ampVecsSignal;
  AmpVecs m_ampVecsBkgnd;
  
  double m_sumBkgWeights;
  double m_numBkgEvents;
  double m_numDataEvents;
};

#endif
