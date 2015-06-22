#if !defined( PARAMETERMANAGER )
#define PARAMETERMANAGER

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
#include <vector>
#include <string>

#include "IUAmpTools/ComplexParameter.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/GaussianBound.h"
#include "MinuitInterface/MIObserver.h"
#include "MinuitInterface/MISubject.h"
#include "IUAmpTools/IntensityManager.h"

class ConfigurationInfo;
class ParameterInfo;

class ParameterManager : MIObserver
{
  
public:
  
  // any negative integer -- used note fixed parameters
  // in an array of indices
  enum { kFixedIndex = -1 };
  
  ParameterManager( MinuitMinimizationManager* minuitManager,
                    IntensityManager* intenManager );
  
  ParameterManager( MinuitMinimizationManager* minuitManager,
                   const vector<IntensityManager*>& intenManagers );
  
  ~ParameterManager();
  
  MinuitMinimizationManager* fitManager() const { return m_minuitManager; }
  
  void setupFromConfigurationInfo( ConfigurationInfo* cfgInfo );
    
  // these functions provide a list of all known parameters, including those that are
  // constrained to other parameters in addition to a covariance matrix that
  // incorporates those constraints
  
  vector< double > parameterValues() const { return m_parValues; }
  vector< string > parameterList() const { return m_parList; }
  map< string, int > parameterIndex() const { return m_parIndex; }
  vector< vector< double > > covarianceMatrix() const { return m_covMatrix; }
  
  bool hasConstraints(const string& ampName) const;
  bool hasParameter(const string& ampName) const;
  
  // this gets called whenever an amplitude parameter changes
  void update( const MISubject* parPtr );
  
protected:
  
  // MPI implementations on the worker nodes need to be able
  // to create a ParameterManager without attaching it to
  // a MinuitMinimizationManager - no one else should be using
  // these constructors.  Use of these constructors requires
  // overriding the other functions below to avoid dereferencing
  // a null pointer to the MinuitMinimizationManager.

  ParameterManager( IntensityManager* intenManager );
  ParameterManager( const vector<IntensityManager*>& intenManager );
  
  // these functions need to be virtual so that parallel implementations
  // can override their functionality correctly since they are called
  // from within setupFromConfigurationInfo
  
  virtual void addProductionParameter( const string& ampName, bool real = false, bool fixed = false );
  virtual void addAmplitudeParameter( const string& ampName, const ParameterInfo* parInfo );
  
  // useful for MPI implementations of ParameterManager
  complex< double >* getProdParPtr( const string& ampName );
  double* getAmpParPtr( const string& parName );
  
  virtual void update( const string& parName );

 private:

  // stop default and copy
  ParameterManager();
  ParameterManager( const ParameterManager& );
  
  ComplexParameter* findParameter(const string& ampName) const;
  void updateParCovariance();
  
  MinuitMinimizationManager* m_minuitManager;
  
  vector< IntensityManager* > m_intenManagers;
  
  vector< double > m_parValues;
  vector< string > m_parList;
  map< string, int > m_parIndex;
  vector< vector< double > > m_covMatrix;
  
  map< string, ComplexParameter* > m_prodParams;
  vector< ComplexParameter* > m_prodPtrCache;
  
  map< string, MinuitParameter* > m_ampParams;
  vector< MinuitParameter* > m_ampPtrCache;
  
  vector< GaussianBound* > m_boundPtrCache;
  
  map <string, vector<string> > m_constraintMap;
};

#endif
