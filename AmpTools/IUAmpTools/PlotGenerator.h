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
#include <map>
#include <string>
#include <complex>
#include <utility>

#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Histogram2D.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigurationInfo.h"

class AmpToolsInterface;
class FitResults;

using namespace std;

class PlotGenerator
{
  
public:
  
  enum { kData = 0, kBkgnd, kGenMC, kAccMC, kNumTypes };
  
  PlotGenerator( const FitResults& fitResults );
  PlotGenerator( ); // empty constructor

  virtual ~PlotGenerator();
  
  // get intensity for all amplitudes that are turned on
  pair< double, double > intensity( bool accCorrected = true ) const;
  
  bool haveAmp( const string& amp ) const;
  
  const ConfigurationInfo* cfgInfo() const { return m_cfgInfo; }
  
  const vector< string >& availablePlots() const { return m_histTitles; }
  
  // the list with reactions and sums
  const vector< string >& fullAmplitudes()  const { return m_fullAmplitudes;  }
  // the list of unique amplitudes
  const vector< string >& uniqueAmplitudes()  const { return m_uniqueAmplitudes;  }
  // the list of unique sums
  const vector< string >& uniqueSums() const { return m_uniqueSums; }
  
  // the reactions themselves
  vector< string > reactions() const;
  
  Histogram* projection( unsigned int projectionIndex, string fsName,
                        unsigned int type );
  
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
   
  Histogram* getHistogram( int index );

protected:
  
  // this function maintained to support older code 
  void bookHistogram( int index, const string& title, Histogram* hist );
  void bookHistogram( int index, Histogram* hist );

  void fillHistogram( int index, double value );
  void fillHistogram( int index, double valueX,double valueY );
  void fillHistogram( int index, vector <double> &data, double weight = 1);
  
  unsigned int getAmpIndex( const string& ampName ) const;
 
private:
 
  // this function should be overridden by the derived class
  // that is written by the user
  virtual void projectEvent( Kinematics* kin ) = 0;
  
  void clearHistograms();
  void fillProjections( const string& reactName, unsigned int type );
  
  void recordConfiguration();
  void buildUniqueAmplitudes();
    
  vector< string >
     stringSplit(const string& str, const string& delimiters = " ") const;
  
  const FitResults& m_fitResults;
  AmpToolsInterface m_ati;
  const ConfigurationInfo* m_cfgInfo;
  
  map< string, NormIntInterface* > m_normIntMap;
  map< string, IntensityManager* > m_intenManagerMap;
  
  // dimension of these vectors is the number of amplitudes (constrained + free)
  vector< complex< double > > m_fitProdAmps;
  vector< complex< double > > m_zeroProdAmps; 
  vector< complex< double > > m_prodAmps; 
  
  // a vector of amplitude parameters
  map< string, double > m_ampParameters;
  
  // these amplitudes include reaction name separated by a '::'
  vector< string > m_fullAmplitudes;
  
  // these lists are independent of reaction
  vector< string > m_uniqueAmplitudes;
  vector< string > m_uniqueSums;
  
  // index of production parameter for a particular (full) amplitude
  map< string, unsigned int > m_ampIndex;
  
  // index of a particular reaction
  map< string, unsigned int > m_reactIndex;
  
  // keep track of which amplitudes and final states are enabled
  map< string, bool > m_ampEnabled;
  map< string, bool > m_sumEnabled;
  map< string, bool > m_reactEnabled;
    
  // a map from ampConfiguration -> final state -> histogram cache
  mutable map< string, map< string, vector< Histogram* > > > m_accMCHistCache;
  mutable map< string, map< string, vector< Histogram* > > > m_genMCHistCache;
  // note that for data and bkgnd the ampConfiguration is useless, but keep the
  // same structure for ease of coding -- bypass first string
  // in the code that fetches projections
  mutable map< string, map< string, vector< Histogram* > > > m_dataHistCache;
  mutable map< string, map< string, vector< Histogram* > > > m_bkgndHistCache;
  
  vector< string > m_histTitles;
  mutable vector< Histogram*> m_histVect,m_histVect_clone;
  mutable double m_currentEventWeight;
  
  string m_currentConfiguration;
};

#endif
