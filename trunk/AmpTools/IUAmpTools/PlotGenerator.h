#if !(defined PLOTGENERATOR)
#define PLOTGENERATOR

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
#include <complex>
#include <utility>

#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigurationInfo.h"

using namespace std;

class PlotGenerator
{
	
public:
  
  enum PlotType { kData = 0, kAccMC, kGenMC, kNumPlotTypes };
  
  PlotGenerator( const ConfigurationInfo* cfgInfo, 
                 const string& parFile );
	
	virtual ~PlotGenerator();
  
  // get intensity for all amplitudes that the plot generator
  // has turned on
  pair< double, double > intensity( bool accCorrected = true ) const;
  // get intensity for a subset of amps **within the same mode**
	pair< double, double > intensity( const vector< string >& amps,
                                    bool accCorrected = true ) const;
  
    // this computes the phase differences between sums of production
    // coefficients -- it is not clear if this has any meaningful signficance
    
	double phaseDiff( const vector< string >& amps1, 
                   const vector< string >& amps2 );
  
  bool haveAmp( const string& amp ) const;
  
  const ConfigurationInfo* cfgInfo() const { return m_cfgInfo; }
  
  // this method overridden by dervied class
  virtual const vector< string >& availablePlots() const = 0;
  
  // the list with reactions and sums
  const vector< string >& fullAmplitudes()  const { return m_fullAmplitudes;  }
  // the list of unique amplitudes
  const vector< string >& uniqueAmplitudes()  const { return m_uniqueAmplitudes;  }
  // the list of unique sums
  const vector< string >& uniqueSums() const { return m_uniqueSums; }
  
  // the modes themselves
  vector< string > reactions() const;
  
  const Histogram& projection( unsigned int projectionIndex, string fsName,
                              PlotType type );
  
  void disableReaction( const string& fsName );
  void enableReaction( const string& fsName );
  
  // these args are the index in the list of unique amplitudes and sums
  void disableAmp( unsigned int uniqueAmpIndex );
  void enableAmp( unsigned int uniqueAmpIndex );
  void disableSum( unsigned int uniqueSumIndex );
  void enableSum( unsigned int uniqueSumIndex );
  
  bool isReactionEnabled( const string& reactName ) const;
  bool isAmpEnabled( unsigned int uniqueAmpIndex ) const;
  bool isSumEnabled( unsigned int uniqueSumIndex ) const;
  
protected:
  
	unsigned int getAmpIndex( const string& ampName ) const;
  unsigned int getErrorMatrixIndex( const string& ampName ) const;

  const AmplitudeManager& ampManager( const string& fsName );
  const NormIntInterface& normIntInterface( const string& fsName );
  
  // this needs to be called by the derived class so that the 
  // base class can setup amp managers and norm ints
  void initialize();
  
private:
  
  // this should be overridden by derived class
  virtual vector< Histogram > fillProjections( const string& fsName, PlotType type ) = 0;
  virtual void registerPhysics( AmplitudeManager* ampManager ) = 0;
  
  void recordConfiguration();
  void buildUniqueAmplitudes();
  
  double phase( const complex< double >& num ) const;
  
  vector< string >
  stringSplit(const string& str, const string& delimiters = " ") const;
  
	const ConfigurationInfo* m_cfgInfo;
  
	map< string, NormIntInterface* > m_normIntMap;
  map< string, AmplitudeManager* > m_ampManagerMap;
  
  // dimension of these vectors is the number of amplitudes (constrained + free)
	vector< complex< double > > m_fitProdAmps;
  vector< complex< double > > m_zeroProdAmps;	
  vector< complex< double > > m_prodAmps;	
  
  // the fit error matrix, dimension:  2 * nFreeProdAmps + nAmpPars
  // imaginary parts of fixed phase production amps are padded with zeros
  vector< vector< double > > m_errorMatrix;
  // a map from amplitude / parameter name to the error matrix row/column
  map< string, unsigned int > m_eMatIndex;
  
  // a vector of amplitude parameters
  map< string, double > m_ampParameters;
  
  // these amplitudes include reaction name separated
  // by a '::'
  vector< string > m_fullAmplitudes;
  
  // these lists are independent of reaction
  vector< string > m_uniqueAmplitudes;
  vector< string > m_uniqueSums;
  
  // index of production parameter for a particular (full) amplitude
	map< string, unsigned int > m_ampIndex;
  
  // keep track of which amplitudes and final states are enabled
  map< string, bool > m_ampEnabled;
  map< string, bool > m_sumEnabled;
  map< string, bool > m_reactEnabled;
  
  // a map from ampConfiguration -> final state -> histogram cache
	mutable map< string, map< string, vector< Histogram > > > m_accMCHistCache;
	mutable map< string, map< string, vector< Histogram > > > m_genMCHistCache;
  // note that for data the ampConfiguration is useless, but keep the
  // same structure for ease of coding -- bypass first string
  // in the code that fetches projections
  mutable map< string, map< string, vector< Histogram > > > m_dataHistCache;
  
  string m_currentConfiguration;
  Histogram m_emptyHist;
};

#endif
