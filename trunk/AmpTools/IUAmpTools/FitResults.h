#ifndef FITRESULTS
#define FITRESULTS

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
#include <string>
#include <map>

class LikelihoodCalculator;
class ParameterManager;
class MinuitMinimizationManager;
class AmplitudeManager;
class NormIntInterface;
class ConfigurationInfo;

using namespace std;

/**
 * A class for storing the information from a fit so that subsequent
 * analysis can be performed.
 *
 * This class captures a snapshot of the relevant information from the
 * AmplitudeManager, LikelihoodCalculator, MinuitMinimizationManager
 * and ParameterManager after a fit.  In addition it has the capability
 * to store and redeliver NormIntInterface objects to provide the
 * normalization integrals for the amplitude.
 */

class FitResults
{
  
public:
  
  FitResults( ConfigurationInfo* cfgInfo,
              vector< AmplitudeManager* > ampManVec,
              map< string, LikelihoodCalculator* > likCalcMap,
              map< string, NormIntInterface* > normIntMap,
              MinuitMinimizationManager* minManager,
              ParameterManager* parManager );

  FitResults( const string& inFile );
  
  ~FitResults();

  void saveResults();
  void writeResults( const string& fileName );
  void loadResults( const string& fileName );
  
  const NormIntInterface* normInt( const string& reactionName ) const;
  const ConfigurationInfo* configInfo() const { return m_cfgInfo; }
    
  double likelihood() const { return m_likelihoodTotal; }
//  double likelihood( const string& reactionName ) const;
  
  int eMatrixStatus() const { return m_eMatrixStatus; }
  int lastMinuitCommandStatus() const { return m_lastCommandStatus; }
  int lastMinuitCommand() const { return m_lastCommand; }
  double minuitPrecision() const { return m_precision; }
  int minuitStrategy() const { return m_strategy; }
  double estDistToMinimum() const { return m_estDistToMin; }
  double bestMinimum() const { return m_bestMin; }
  
  const vector< string >& parList() const { return m_parNames; }
//  double parValue( const string& parName ) const; 
//  double parError( const string& parName ) const;

  const vector< vector< double > >& errorMatrix() const { return m_covMatrix; }
//  double covariance( const string& par1, const string& par2 );
  
  const vector< string >& reactionList() const { return m_reactionNames; }
//  const vector< string >& ampList( const string& reaction ) const;
//  double ampScale( const string& reaction, const string& amplitude );
//  string ampScaleName( const string& reaction, const string& amplitude );
  
protected:
  
  // disable the copy constuctor and assignment operator until memory
  // handling issues are fixed
  
  FitResults( const FitResults& fitResults );
  FitResults();
  FitResults& operator=( const FitResults& results );
  
private:
  
  void recordAmpSetup();
  void recordLikelihood();
  void recordFitStats();
  void recordParameters();

  // amplitude manager info
  
  int m_numReactions;
  vector< string > m_reactionNames;
  vector< int > m_numAmps;
  vector< vector< string > > m_ampNames;
  vector< vector< string > > m_ampScaleNames;
  vector< vector< double > > m_ampScaleValues;
  
  map< string, int > m_reacIndex;
  vector< map< string, int > > m_ampIndex;
  
  // likelihood info
  
  map< string, double > m_likelihoodMap;
  double m_likelihoodTotal;
  
  // info from the fitter
  
  int m_eMatrixStatus;
  int m_lastCommandStatus;
  int m_lastCommand;
  double m_precision;
  int m_strategy;
  double m_estDistToMin;
  double m_bestMin;
  
  // info about the parameters
  
  vector< string > m_parNames;
  vector< double > m_parValues;
  vector< vector< double > > m_covMatrix;
  
  map< string, int > m_parIndex;
  
  // hold pointers to the classes that are needed to fetch the results
  
  ConfigurationInfo* m_cfgInfo;
  vector< AmplitudeManager* > m_ampManVec;
  map< string, LikelihoodCalculator* > m_likCalcMap;
  map< string, NormIntInterface* > m_normIntMap;
  MinuitMinimizationManager* m_minManager;
  ParameterManager* m_parManager;
  
  bool m_createdFromFile;
  
};


#endif
