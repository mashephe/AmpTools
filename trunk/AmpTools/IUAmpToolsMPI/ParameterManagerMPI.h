#if !defined(PARAMETERMANAGERMPI)
#define PARAMETERMANAGERMPI

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
#include <map>
#include <complex>

#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpToolsMPI/LikelihoodManagerMPI.h"
#include "MinuitInterface/MIObserver.h"

using namespace std;

class ParameterManagerMPI : public ParameterManager
{
  
public:
  
  enum { kMaxNameLength = 200 };
  
  // since there is only one MinuitMinimizationManager we must add
  // an additional constructor for the worker nodes
  
  // this constructor is called on the head node
  ParameterManagerMPI( MinuitMinimizationManager* minuitManager,
                      IntensityManager* intenManager );
  ParameterManagerMPI( MinuitMinimizationManager* minuitManager,
                      const vector< IntensityManager* >& intenManagers );
  
  // this constructor should be called on the worker nodes
  ParameterManagerMPI( IntensityManager* intenManager );
  ParameterManagerMPI( const vector< IntensityManager* >& intenManagers );
  
  ~ParameterManagerMPI();
  
  // override these functions to do the appropriate thing depending
  // on whether this instance is on the head or worker node
  void addProductionParameter( const string& termName, bool real = false, bool fixed = false );
  void addAmplitudeParameter( const string& termName, const ParameterInfo* parInfo );
  
  // the likelihood calculator will need to call this routine in order
  // to update the parameters in advance of the likelihood calculation
  void updateParameters();
  
  // the likelihood calculator on the worker nodes will call this routine
  // in order to update a changed ampitude parameter -- the update routine
  // on the master (function below) will call this directly with the
  // parameter name
  void updateAmpParameter( const string& parName = "" );
  
protected:
  
  // this overrides the base class function
  void update( const string& parName );
  
private:
  
  void setupMPI();
  
  vector<IntensityManager*> m_intenManagers;
  
  int m_rank;
  int m_numProc;
  bool m_isMaster;
  
  map< string, complex< double >* > m_prodParMap;
  map< string, double* > m_ampParMap;
  
};

#endif
