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

#include "IUAmpTools/ComplexParameter.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/GaussianBound.h"
#include "MinuitInterface/MIObserver.h"
#include "MinuitInterface/MISubject.h"
#include "IUAmpTools/AmplitudeManager.h"

class ConfigurationInfo;

class ParameterManager : MIObserver
{
	
 public:

  ParameterManager( MinuitMinimizationManager& minuitManager,
                   AmplitudeManager* ampManager );

  ParameterManager( MinuitMinimizationManager& minuitManager,
		    const vector<AmplitudeManager*>& ampManager );
	
  ~ParameterManager();

  MinuitMinimizationManager& fitManager() const { return m_minuitManager; }
    
  void setupFromConfigurationInfo( ConfigurationInfo* cfgInfo );
  
  // these functions need to be virtual so that parallel implementations
  // can override their functionality correctly since they are called
  // from within setupFromConfigurationInfo
  virtual void addProductionParameter( const string& ampName, bool real = false );  
  virtual void addAmplitudeParameter( const string& ampName, const ParameterInfo* parInfo );
  
  void writeParameters( ofstream& file ) const;

  void addConstraintMap(const map<string, vector<string> >& constraintMap)
  {m_constraintMap = constraintMap;}

  bool hasConstraints(const string& ampName) const;

  bool hasParameter(const string& ampName) const;

  ComplexParameter* findParameter(const string& ampName) const;
  
  void update( const MISubject* parPtr );

 protected:

  // useful for MPI implementations of ParameterManager
  complex< double >* getProdParPtr( const string& ampName );
  double* getAmpParPtr( const string& parName );
	
 private:
	
  // stop default
  ParameterManager();
  ParameterManager( const ParameterManager& );
		
  MinuitMinimizationManager& m_minuitManager;

  vector< AmplitudeManager* > m_ampManagers;

  map< string, ComplexParameter* > m_prodParams;
  vector< ComplexParameter* > m_prodPtrCache;
  
  map< string, MinuitParameter* > m_ampParams;
  vector< MinuitParameter* > m_ampPtrCache;
  
  vector< GaussianBound* > m_boundPtrCache;
  
  map <string, vector<string> > m_constraintMap;
  
  double m_scale;
};

#endif
