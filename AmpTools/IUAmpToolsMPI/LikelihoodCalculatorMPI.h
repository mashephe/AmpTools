#if !defined(LIKELIHOODCALCULATORMPI)
#define LIKELIHOODCALCULATORMPI

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

#include "IUAmpToolsMPI/MPITag.h"
#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpTools/LikelihoodCalculator.h"

class LikelihoodManagerMPI;
class MISubject;

/**
 * This class is implements the parallel likelihood calculation in the MPI
 * enviornment.  The general idea is that MPI communications are confined
 * to this class, which utilizes the base class for doing the actual
 * computation.  This avoids replication of likelihood computation in this
 * class and confines MPI impliementation to just the methods in this class.
 *
 * \ingroup IUAmpToolsMPI
 */

class LikelihoodCalculatorMPI : public LikelihoodCalculator
{
  
  friend class LikelihoodManagerMPI;

 public:

  /**
   * This enum allows multiple instances of the likelihoodCalculatorMPI to
   * synchronize themselves across many processes.  One instance of the
   * likelihoodCalculatorMPI is created for every reaction in the fit.  The
   * code uses the MPI tags to send commands (via LikelihoodManagerMPI) to
   * each instance.  Therefore the instance numbering starts with the first
   * tag available after the commands listed in MPITags.
   */
  
  enum { kFirstId = MPITag::kMaxTags };

  /**
   * This is the constructor.  On all nodes setupMPI() is called first to
   * define the rank and number of processes.  On the woker nodes, each
   * instance of LikelihoodCalulatorMPI registers itself with the 
   * LikelihoodManagerMPI.  Each instance on the worker nodes has a unique
   * ID which is derived from the cuerrent value of the static member data
   * m_idCounter.
   *
   * Arguments to this constructor are similar to those for LikelihoodCalculator.
   * This helps to maintain a familar interface for the user.
   * 
   * \param[in] ampManager A reference to the amplitude manager that this
   * likelihood calculator should use for amplitude calculations
   * \param[in] normInt A reference to the appropriate normalization integral
   * interface
   * \param[in] dataReader A reference to the DataReader that will supply
   * the data
   * \param[in] parManager A reference to the ParameterManagerMPI instance
   * that has been created to manage parameters on the node.
   *
   * \see LikelihoodCalculator
   * \see LikelihoodManagerMPI
   */
  
  LikelihoodCalculatorMPI( const IntensityManager& intenManager,
                           const NormIntInterface& normInt,
                           DataReader* dataReaderSignal,
                           DataReader* dataReaderBkgnd,
                           ParameterManagerMPI& parManager );
  
  /**
   * This is the destructor.  When the instance of LikelihoodCalculatorMPI
   * is destroyed on the master node, it sends the exit flag (via the
   * LikelihoodManagerMPI) to the corresponding LikelihoodCalculatorMPI that
   * exists on the worker nodes.
   */
  
  ~LikelihoodCalculatorMPI();

  /**
   * The following operator should only ever be called on the master.  It 
   * overrides the operator()() that is called from MIFunctionContribution class.
   * Using the LikelihoodManagerMPI, it directs all of the workers to compute
   * and then send partial contributions of the sum of log intensities.  The
   * routine on the master collects and sums the contributions from the workers.
   * It then, if necessary, triggers an update of the normalization integral 
   * calculculation on the workers.  Finally it provides the new 
   * -2 ln( likelihood ) for the fit.
   */
  double operator()();
  
private:

  // the following functions are used by the LikelihoodManager to trigger
  // portions of the likeihood calculation -- they should only be called
  // on the worker nodes
  void updateParameters();
  void updateAmpParameter();
  void computeLikelihood();

  static int m_idCounter;

  void setupMPI();
  
  const IntensityManager& m_intenManager;

  ParameterManagerMPI& m_parManager;
  int m_thisId;

  int m_rank;
  int m_numProc;
  bool m_isMaster;
  bool m_firstPass;
};

#endif
