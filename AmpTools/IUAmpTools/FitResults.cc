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

#include <iostream>

#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"

FitResults::FitResults( ConfigurationInfo* cfgInfo,
                        vector< AmplitudeManager* > ampManVec,
                        map< string, LikelihoodCalculator* > likCalcMap,
                        map< string, NormIntInterface* > normIntMap,
                        MinuitMinimizationManager* minManager,
                        ParameterManager* parManager) :
m_likelihoodTotal( 0 ),
m_cfgInfo( cfgInfo ),
m_ampManVec( ampManVec ),
m_likCalcMap( likCalcMap ),
m_normIntMap( normIntMap ),
m_minManager( minManager ),
m_parManager( parManager ),
m_createdFromFile( false ),
m_warnedAboutFreeParams( false ),
m_isValid( true ){
    
}

FitResults::FitResults( const string& inFile ) :
m_createdFromFile( true ),
m_warnedAboutFreeParams( false ),
m_isValid( false ){
  
  loadResults( inFile );
}

FitResults::~FitResults() {
  
  if( m_createdFromFile ) {
   
    for( map< string, NormIntInterface* >::iterator mapItr = m_normIntMap.begin();
         mapItr != m_normIntMap.end();
         ++mapItr ){
    
      delete mapItr->second;
    }
  }
}

double
FitResults::likelihood( const string& reaction ) const {

  map< string, double >::const_iterator reac = m_likelihoodMap.find( reaction );
  
  if( reac == m_likelihoodMap.end () ){
    
    cout << "FitResults ERROR:  request for likelihood of unknown reaction: " << reaction << endl;
    return sqrt( -1 );
  }
  
  return reac->second;
}

pair< double, double >
FitResults::intensity( bool accCorrected ) const {

  // return the intensity for all amplitudes
  return intensity( ampList(), accCorrected );
}

pair< double, double >
FitResults::intensity( const vector< string >& amplitudes, bool accCorrected ) const {
  
  // intensity = sum_a sum_a' s_a s_a' V_a V*_a' NI( a, a' )
    
  // a subset of the larger error matrix
  vector< vector< double > > errorMatrix;

  // a vector for the derivatives of the intensity with respect to the
  // real and imaginary parts of the production amplitudes and the
  // scale parameter for the amplitude
  vector< double > deriv( 3 * amplitudes.size() );
  vector< double > parName( 3 * amplitudes.size() );

  double intensity = 0;

  // check for free parameters and print an appropriate warning message
  // about the accuracy of the intensity errors since we don't yet have
  // the ability to numerically compute the derivatives of the
  // normalization integrals with respect to the free parameters
 
  map< string, double > ampParameters = ampParMap();
  for( map< string, double >::iterator par = ampParameters.begin();
       par != ampParameters.end();
       ++par ){
    
    // there are parameters -- check errors to see if they
    // were floating in the fit
    int parIndex = m_parIndex.find( par->first )->second;
    if( fabs( m_covMatrix[parIndex][parIndex] ) > 1E-20 ){
    
      if( !m_warnedAboutFreeParams )
      
      cout
      << "***************************************************************" << endl
      << "* WARNING!:  You are calculating an intensity that depends on *" << endl
      << "*   parameters that were floating in the fit.  Unless the par-*" << endl
      << "*   ameters are amplitude scale factors only, THE ERROR ON    *" << endl
      << "*   THE INTENSITY SHOULD BE CONSIDERED UNRELIABLE because the *" << endl
      << "*   software is assuming that the derivatives of the normali- *" << endl
      << "*   zation integrals with respect to the parameters are equal *" << endl
      << "*   to zero.  (Fixing this involves numercially computing the *" << endl
      << "*   derivatives, which will be implemented in future versions *" << endl
      << "*   of AmpTools.)  For now it is recommended that you repeat  *" << endl
      << "*   the fit fixing the free parameter(s) at +- 1 sigma of the *" << endl
      << "*   parameter uncertainty to systematically probe how the     *" << endl
      << "*   uncertainty on the intensity depends on the uncertainty   *" << endl
      << "*   of the parameter.                                         *" << endl
      << "***************************************************************" << endl
      << endl;
      
      m_warnedAboutFreeParams = true;
    }
  }
  
  // now calculate the vector of partial derivatives of the intensity
  // with respect to the real, imaginary, and scale factors for each
  // production amplitude
  
  for( vector< string >::const_iterator amp = amplitudes.begin();
      amp != amplitudes.end(); ++amp ){

    // the vector of amplitudes should be "full" amplitude names
    // that include the reaction, sum, and amplitude name
    vector<string> ampNameParts = stringSplit( *amp, "::" );
    assert( ampNameParts.size() == 3 );
    string reaction = ampNameParts[0];
    
    // we need to build a small covriance matrix for the calculation
    // use lower case to denote the index in the small matrix while
    // upper case denotes the index in the full matrix
    
    int iRe = 3 * ( amp - amplitudes.begin() );
    int iIm = iRe + 1;
    int iScale = iRe + 2;
    
    int IRe    = m_parIndex.find( realProdParName( *amp ) )->second;
    int IIm    = m_parIndex.find( imagProdParName( *amp ) )->second;
    
    // if the name of the scale parameter doesn't appear in the list of
    // free parameters set the index to -1 and watch for this below
    map< string, int >::const_iterator idx1 = m_parIndex.find( ampScaleName( *amp ) );
    int IScale = ( idx1 == m_parIndex.end() ? -1 : idx1->second );
    
    errorMatrix.push_back( vector< double >( 3 * amplitudes.size() ) );
    errorMatrix.push_back( vector< double >( 3 * amplitudes.size() ) );
    errorMatrix.push_back( vector< double >( 3 * amplitudes.size() ) );
    
    deriv[iRe]    = 0;
    deriv[iIm]    = 0;
    deriv[iScale] = 0;
    
    for( vector< string >::const_iterator conjAmp = amplitudes.begin();
        conjAmp != amplitudes.end(); ++conjAmp ) {

      vector<string> conjAmpNameParts = stringSplit( *conjAmp, "::" );
      assert( conjAmpNameParts.size() == 3 );
      string conjReaction = ampNameParts[0];

      int jRe = 3 * ( conjAmp - amplitudes.begin() );
      int jIm = jRe + 1;
      int jScale = jRe + 2;

      int JRe    = m_parIndex.find( realProdParName( *conjAmp ) )->second;
      int JIm    = m_parIndex.find( imagProdParName( *conjAmp ) )->second;
      
      // if the name of the scale parameter doesn't appear in the list of
      // free parameters set the index to -1 and watch for this below
      map< string, int >::const_iterator idx2 = m_parIndex.find( ampScaleName( *conjAmp ) );
      int JScale = ( idx2 == m_parIndex.end() ? -1 : idx2->second );
    
      // if the two amplitudes are from different reactions then they do not add coherently
      // and we should set the integral of A A* to zero
      // for amplitudes from different sums within the same reaction, this happens
      // automatically in the normalization integral interface
      complex< double > ampInt;
      if( strcmp( reaction.c_str(), conjReaction.c_str() ) == 0 ){
        
        if (accCorrected)  ampInt = m_normIntMap.find(reaction)->second->ampInt( *amp, *conjAmp );
        else               ampInt = m_normIntMap.find(reaction)->second->normInt( *amp, *conjAmp );
      }
      else{
        
        ampInt = complex< double >( 0 , 0 );
      }
      
      // copy a 3 x 3 block into the small error matrix
      
      errorMatrix[iRe][jRe] = m_covMatrix[IRe][JRe];
      errorMatrix[iRe][jIm] = m_covMatrix[IRe][JIm];
      errorMatrix[iRe][jScale] =
         ( JScale == -1 ? 0 : m_covMatrix[IRe][JScale] );
      errorMatrix[iIm][jRe] = m_covMatrix[IIm][JRe];
      errorMatrix[iIm][jIm] = m_covMatrix[IIm][JIm];
      errorMatrix[iIm][jScale] =
         ( JScale == -1 ? 0 : m_covMatrix[IIm][JScale] );
      errorMatrix[iScale][jRe] =
         ( IScale == -1 ? 0 : m_covMatrix[IScale][JRe] );
      errorMatrix[iScale][jIm] =
         ( IScale == -1 ? 0 : m_covMatrix[IScale][JIm] );
      errorMatrix[iScale][jScale] =
         ( IScale == -1 || JScale == -1 ? 0 : m_covMatrix[IScale][JScale] );
      
      deriv[iRe]    += 2 * ampScale( *amp ) * ampScale( *conjAmp ) *
        ( m_parValues[JRe] * real( ampInt ) + m_parValues[JIm] * imag( ampInt ) );
      deriv[iIm]    += 2 * ampScale( *amp ) * ampScale( *conjAmp ) *
        ( m_parValues[JIm] * real( ampInt ) - m_parValues[JRe] * imag( ampInt ) );
  
      double intensityContrib = ampScale( *amp ) * ampScale( *conjAmp ) *
        ( ( m_parValues[IRe] * m_parValues[JRe] + m_parValues[IIm] * m_parValues[JIm] ) * real( ampInt ) -
          ( m_parValues[IIm] * m_parValues[JRe] - m_parValues[IRe] * m_parValues[JIm] ) * imag( ampInt ) );
      
      deriv[iScale] += ( intensityContrib / ampScale( *amp ) );
      
      intensity += intensityContrib;
    }
  }
  
  // now compute the error
  double variance = 0;
  for( unsigned int i = 0; i < deriv.size(); ++i ){
    for( unsigned int j = 0; j < deriv.size(); ++j ){
      
      variance += deriv[i] * deriv[j] * errorMatrix[i][j];
    }
  }
  
  return pair< double, double >( intensity, sqrt( variance ) );
}

pair< double, double >
FitResults::phaseDiff( const string& amp1, const string& amp2 ) const {
  
  vector<string> ampNameParts = stringSplit( amp1, "::" );
  assert( ampNameParts.size() == 3 );
  string reaction1 = ampNameParts[0];
  string sum1 = ampNameParts[1];
  
  ampNameParts = stringSplit( amp2, "::" );
  assert( ampNameParts.size() == 3 );
  string reaction2 = ampNameParts[0];
  string sum2 = ampNameParts[1];

  if( ( strcmp( reaction1.c_str(), reaction2.c_str() ) != 0 ) ||
      ( strcmp( sum1.c_str(), sum2.c_str() ) != 0 ) ){
    
    cout << "FitResults WARNING:: request to compute phase difference of " << endl
         << "                     amplitudes from different sums or different" << endl
         << "                     reactions which is not meaningful, returning 0. " << endl
         << "    amp1: " << amp1 << endl
         << "    amp2: " << amp2 << endl;
    
    return pair< double, double >( 0, 0 );
  }

  
  // The phase difference depends only on two parameters, the real and imaginary
  // components of the two amplitudes.  It is independent of the scale of the
  // amplitudes.
  vector< int > idx( 4 );
  idx[0] = m_parIndex.find( amp1 + "_re" )->second;
  idx[1] = m_parIndex.find( amp1 + "_im" )->second;
  idx[2] = m_parIndex.find( amp2 + "_re" )->second;
  idx[3] = m_parIndex.find( amp2 + "_im" )->second;

  // this makes the code a little easier to read
  double a1Re = m_parValues[idx[0]];
  double a1Im = m_parValues[idx[1]];
  double a2Re = m_parValues[idx[2]];
  double a2Im = m_parValues[idx[3]];

  vector< double > pDeriv( 4 );
  pDeriv[0] = ( -a1Im / ( a1Re * a1Re + a1Im * a1Im ) );
  pDeriv[1] = (  a1Re / ( a1Re * a1Re + a1Im * a1Im ) );
  pDeriv[2] = (  a2Im / ( a2Re * a2Re + a2Im * a2Im ) );
  pDeriv[3] = ( -a2Re / ( a2Re * a2Re + a2Im * a2Im ) );
  
  double pVar = 0;
  for( unsigned int i = 0; i < pDeriv.size(); ++i ){
    for( unsigned int j = 0; j < pDeriv.size(); ++j ){
      
      pVar += pDeriv[i] * pDeriv[j] * m_covMatrix[idx[i]][idx[j]];
    }
  }
  
  return pair< double, double >( arg( complex< double > ( a1Re, a1Im ) ) -
                                 arg( complex< double > ( a2Re, a2Im ) ),
                                 sqrt(pVar) );
}

string
FitResults::realProdParName( const string& amplitude ) const {
  
  string parName = amplitude + "_re";
  
  // be sure the parameter actually exists before returning its name
  // if this fails, then a bogus amplitude name was passed in
  assert( m_parIndex.find( parName ) != m_parIndex.end() );
  
  return parName;
}

string
FitResults::imagProdParName( const string& amplitude ) const {
  
  string parName = amplitude + "_im";
  
  // be sure the parameter actually exists before returning its name
  // if this fails, then a bogus amplitude name was passed in
  assert( m_parIndex.find( parName ) != m_parIndex.end() );
  
  return parName;
}

string
FitResults::ampScaleName( const string& amplitude ) const {
  
  // the vector of amplitudes should be "full" amplitude names
  // that include the reaction, sum, and amplitude name
  
  vector<string> ampNameParts = stringSplit( amplitude, "::" );
  assert( ampNameParts.size() == 3 );
  string reaction = ampNameParts[0];
  
  map< string, int >::const_iterator reactIndexPair = m_reacIndex.find( reaction );
  
  if( reactIndexPair == m_reacIndex.end() ){
    
    cout << "FitResults::ampScaleName ERROR:: no such reaction: " << reaction << endl;
    assert( false );
  }
  
  map< string, int >::const_iterator ampIndexPair =
     m_ampIndex.at(reactIndexPair->second).find( amplitude );
  
  if( ampIndexPair == m_ampIndex.at(reactIndexPair->second).end() ){
    
    cout << "FitResults::ampScaleName ERROR:: no such amplitude: " << amplitude << endl;
    assert( false );
  }
  
  return m_ampScaleNames.at( reactIndexPair->second ).at( ampIndexPair->second );
}

map< string, complex< double > >
FitResults::ampProdParMap() const {
  
  map< string, complex< double > > ampMap;
  
  vector< string > amps = ampList();
  for( vector< string >::iterator amp = amps.begin();
       amp != amps.end();
       ++amp ){

    ampMap[*amp] = productionParameter( *amp );
  }

  return ampMap;
}

map< string, double >
FitResults::ampScaleParMap() const {
  
  map< string, double > ampMap;
  
  vector< string > amps = ampList();
  for( vector< string >::iterator amp = amps.begin();
      amp != amps.end();
      ++amp ){
    
    ampMap[*amp] = ampScale( *amp );
  }
  
  return ampMap;
}

map< string, double >
FitResults::ampParMap() const {
  
  map< string, double > parMap;
  
  vector< ParameterInfo* > parList = m_cfgInfo->parameterList();
  
  for( vector< ParameterInfo* >::iterator par = parList.begin();
       par != parList.end(); ++par ){
    
    if( (**par).fixed() ) continue;
    
    parMap[(**par).parName()] = parValue( (**par).parName() );
  }
  
  return parMap;
}

complex< double >
FitResults::productionParameter( const string& ampName ) const {
    
  int iRe = m_parIndex.find( realProdParName( ampName ) )->second;
  int iIm = m_parIndex.find( imagProdParName( ampName ) )->second;

  return complex< double >( m_parValues[iRe], m_parValues[iIm] );
}


complex< double >
FitResults::scaledProductionParameter( const string& ampName ) const {
  
  return productionParameter( ampName ) * ampScale( ampName );
}

double
FitResults::parValue( const string& parName ) const {
  
  map< string, int >::const_iterator parIndexPair = m_parIndex.find( parName );
  
  if( parIndexPair == m_parIndex.end() ){
    
    cout << "FitResults:: ERROR:  request for unknown parameter " << parName << endl
         << "                     returning nan." << endl;
    
    return sqrt( -1 );
  }
  
  return m_parValues[parIndexPair->second];
}

double
FitResults::parError( const string& parName ) const {

  return sqrt( covariance( parName, parName ) );
}

double
FitResults::covariance( const string& par1, const string& par2 ) const {
  
  map< string, int >::const_iterator par1Pair = m_parIndex.find( par1 );
  map< string, int >::const_iterator par2Pair = m_parIndex.find( par2 );
  
  if( par1Pair == m_parIndex.end() || par2Pair == m_parIndex.end() ){
    
    cout << "FitResults:: ERROR:  request for covaraince of unkown parameters "
         << par1 << ", " << par2 << endl
         << "                     returning nan" << endl;
    return sqrt( -1 );
  }
  
  return m_covMatrix[par1Pair->second][par2Pair->second];
}

double
FitResults::ampScale( const string& amplitude ) const {
  
  // this needs a full amplitude name:  react::sum::amp
  
  vector<string> ampNameParts = stringSplit( amplitude, "::" );
  assert( ampNameParts.size() == 3 );
  string reaction = ampNameParts[0];
  
  map< string, int >::const_iterator reactIndexPair = m_reacIndex.find( reaction );
  
  if( reactIndexPair == m_reacIndex.end() ){
    
    cout << "FitResults::ampScaleName ERROR:: no such reaction: " << reaction << endl;
    assert( false );
  }
  
  map< string, int >::const_iterator ampIndexPair =
  m_ampIndex.at(reactIndexPair->second).find( amplitude );
  
  if( ampIndexPair == m_ampIndex.at(reactIndexPair->second).end() ){
    
    cout << "FitResults::ampScaleName ERROR:: no such amplitude: " << amplitude << endl;
    assert( false );
  }
  
  return m_ampScaleValues.at( reactIndexPair->second ).at( ampIndexPair->second );
}
  
vector< string >
FitResults::ampList() const {
  
  vector< string > list;
  
  for( vector< string >::const_iterator reacItr = m_reactionNames.begin();
      reacItr != m_reactionNames.end();
      ++reacItr ){
    
    vector< string > thisReacList = ampList( *reacItr );
    
    for( vector< string >::iterator ampItr = thisReacList.begin();
         ampItr != thisReacList.end();
         ++ampItr ){
      
      list.push_back( *ampItr );
    }
  }
  
  return list;
}

vector< string >
FitResults::ampList( const string& reaction ) const {
  
  map< string, int >::const_iterator reacIndex = m_reacIndex.find( reaction );
  
  if( reacIndex == m_reacIndex.end() ){
    
    cerr << "FitResults ERROR:: unkown reaction: " << reaction << endl;
    assert( false );
  }
  
  return m_ampNames[reacIndex->second];
}

void
FitResults::saveResults() {
  
  if( !m_createdFromFile ){
    
    recordAmpSetup();
    recordLikelihood();
    recordParameters();
    recordFitStats();
  }
}

void
FitResults::writeResults( const string& outFile ) const {
  
  ofstream output( outFile.c_str() );
  
  output.precision( 15 );
  
  output << "*** DO NOT EDIT THIS FILE - IT IS FORMATED FOR INPUT ***" << endl;
  output << "+++ Reactions, Amplitudes, and Scale Parameters +++" << endl;
  output << "  " << m_numReactions << endl;
  for( int i = 0; i < m_numReactions; ++i ){

    output << "  " <<m_reactionNames[i] << "\t" << m_numAmps[i] << endl;
    
    for( int j = 0; j < m_ampNames[i].size(); ++j ) {
      
      output << "  " <<m_ampNames[i][j] << "\t" 
             << m_ampScaleNames[i][j] << "\t"
             << m_ampScaleValues[i][j] << endl;
    }
  }
  
  output << "+++ Likelihood Total and Partial Sums +++" << endl;
  output << "  " <<m_likelihoodTotal << endl;
  for( int i = 0; i < m_numReactions; ++i ){
    
    output << "  " <<m_reactionNames[i] << "\t" 
           << m_likelihoodMap.find(m_reactionNames[i])->second << endl;
  }
  
  output << "+++ Fitter Information +++" << endl;
  output << "  " << "lastMinuitCommand\t" << m_lastCommand << endl;
  output << "  " << "lastMinuitCommandStatus\t" << m_lastCommandStatus << endl;
  output << "  " << "eMatrixStatus\t" << m_eMatrixStatus << endl;
  output << "  " << "minuitPrecision\t" << m_precision << endl;
  output << "  " << "minuitStrategy\t" << m_strategy << endl;
  output << "  " << "estDistToMinimum\t" << m_estDistToMin << endl;
  output << "  " << "bestMinimum\t" << m_bestMin << endl;
  
  output << "+++ Parameter Values and Errors +++" << endl;
  output << "  " << m_parNames.size() << endl;
  for( int i = 0; i < m_parNames.size(); ++i ){
    
    output << "  " <<m_parNames[i] << "\t" << m_parValues[i] << endl;
  }
  
  for( int i = 0; i < m_parNames.size(); ++i ){
    for( int j = 0; j < m_parNames.size(); ++j ){
      
      output << "  " << m_covMatrix[i][j] << "\t";
    }
    
    output << endl;
  }
  
  // here we will use the NormIntInterface rather than replicating the
  // functionality that is already found there

  output << "+++ Normalization Integrals +++" << endl;

  for( int i = 0; i < m_reactionNames.size(); ++i ){
   
    string reac = m_reactionNames[i];
    output << "  " << reac << endl;
    
    const NormIntInterface* ni = m_normIntMap.find(reac)->second;
    
    if( m_createdFromFile || !ni->hasAccessToMC() ){
      
      ni->exportNormIntCache( output );
    }
    else{
      ni->forceCacheUpdate();
      ni->exportNormIntCache( output,
                              m_ampManVec[i]->ampsAreRenormalized() );
    }
  }
  
  output << "+++ Below these two lines is a config file that is   +++" << endl;
  output << "+++ functionally equivalent to that used in the fit. +++" << endl;
  output << *m_cfgInfo;
  
  output.close();
}

void
FitResults::loadResults( const string& inFile ){
  
  enum { kMaxLine = 256 };
  char line[kMaxLine];
  string tmp;
  
  ifstream input( inFile.c_str() );
  
  if( input.fail() ){
    
    cerr << "ERROR::  FitResults file does not exist: " << inFile << endl;
    return;
  }
  
  input.getline( line, kMaxLine ); // top message
  input.getline( line, kMaxLine ); // amp manager heading
  
  input >> m_numReactions;
  
  m_reactionNames.resize( m_numReactions );
  m_numAmps.resize( m_numReactions );
  m_ampNames.resize( m_numReactions );
  m_ampScaleNames.resize( m_numReactions );
  m_ampScaleValues.resize( m_numReactions );
  m_ampIndex.resize( m_numReactions );
  
  for( int i = 0; i < m_numReactions; ++i ){
    
    input >> m_reactionNames[i] >> m_numAmps[i];
        
    m_reacIndex[m_reactionNames[i]] = i;
    
    m_ampNames[i].resize( m_numAmps[i] );
    m_ampScaleNames[i].resize( m_numAmps[i] );
    m_ampScaleValues[i].resize( m_numAmps[i] );
    
    for( int j = 0; j < m_ampNames[i].size(); ++j ) {
      
      input >> m_ampNames[i][j] 
            >> m_ampScaleNames[i][j]
            >> m_ampScaleValues[i][j];
      
      m_ampIndex[i][m_ampNames[i][j]] = j;
    }
  }
  
  // one getline clears the newline waiting in the buffer
  input.getline( line, kMaxLine ); // likelihood heading
  input.getline( line, kMaxLine ); // likelihood heading
  
  input >> m_likelihoodTotal;
  
  double val;
  for( int i = 0; i < m_numReactions; ++i ){
    
    input >> tmp >> val;
    m_likelihoodMap[tmp] = val;
  }
  
  input.getline( line, kMaxLine ); // fit info heading
  input.getline( line, kMaxLine ); // fit info heading

  input >> tmp >> m_lastCommand;
  input >> tmp >> m_lastCommandStatus;
  input >> tmp >> m_eMatrixStatus;
  input >> tmp >> m_precision;
  input >> tmp >> m_strategy;
  input >> tmp >> m_estDistToMin;
  input >> tmp >> m_bestMin;
  
  input.getline( line, kMaxLine ); // parameters heading
  input.getline( line, kMaxLine ); // parameters heading

  int nPar;
  input >> nPar;
  m_parNames.resize( nPar );
  m_parValues.resize( nPar );
  m_covMatrix.resize( nPar );
  
  for( int i = 0; i < m_parNames.size(); ++i ){
    
    input >> m_parNames[i] >> m_parValues[i];
    m_parIndex[m_parNames[i]] = i;
  }
  
  for( int i = 0; i < m_parNames.size(); ++i ){
    
    m_covMatrix[i].resize( nPar );
    for( int j = 0; j < m_parNames.size(); ++j ){
      
      input >> m_covMatrix[i][j];
    }    
  }
    
  // here we will use the NormIntInterface rather than replicating the
  // functionality that is already found there
  
  input.getline( line, kMaxLine ); // norm int heading
  input.getline( line, kMaxLine ); // norm int heading
  
  for( int i = 0; i < m_numReactions; ++i ){
    
    string reac;
    input >> reac;
    m_normIntMap[reac] = new NormIntInterface();
    input >> (*m_normIntMap[reac]);
  }    

  // now read back in the ConfigurationInfo using the ConfigFileParser
  input.getline( line, kMaxLine ); // cfg info heading
  input.getline( line, kMaxLine ); // cfg info heading
  input.getline( line, kMaxLine ); // cfg info heading
  
  ConfigFileParser cfgParser( input );
  m_cfgInfo = cfgParser.getConfigurationInfo();
  
  input.close();
  
  // set this to true after the load completes - this allows testing
  // of validity when reading results from a file
  
  m_isValid = true;
}

void
FitResults::recordAmpSetup(){
  
  m_numReactions = m_ampManVec.size();
  
  m_numAmps.clear();
  m_reactionNames.clear();
  m_ampNames.clear();
  m_ampScaleNames.clear();
  m_ampScaleValues.clear();
  
  for( vector< AmplitudeManager* >::iterator ampMan = m_ampManVec.begin();
       ampMan != m_ampManVec.end();
       ++ampMan ){
  
    m_reactionNames.push_back( (**ampMan).reactionName() );
    
    vector< string > ampNames = (**ampMan).getAmpNames();
    
    m_ampNames.push_back( ampNames );
    m_numAmps.push_back( ampNames.size() );
    
    vector< string > ampScaleNames( 0 );
    vector< double > ampScaleValues( 0 );
    
    for( int i = 0; i < ampNames.size(); ++i ){
      
      ampScaleNames.push_back( (**ampMan).getScale( ampNames[i] ).name() );
      ampScaleValues.push_back( (**ampMan).getScale( ampNames[i] ) );
    }
    
    m_ampScaleNames.push_back( ampScaleNames );
    m_ampScaleValues.push_back( ampScaleValues );
  }
}

void
FitResults::recordLikelihood(){

  m_likelihoodTotal = 0;
  
  for( map< string, LikelihoodCalculator* >::iterator 
           mapItr = m_likCalcMap.begin();
       mapItr != m_likCalcMap.end();
       ++mapItr ){
    
    double thisLikVal = (*(mapItr->second))();
    
    m_likelihoodTotal += thisLikVal;
    m_likelihoodMap[mapItr->first] = thisLikVal;
  }
}

void
FitResults::recordParameters(){
  
  m_parNames = m_parManager->parameterList();
  m_parValues = m_parManager->parameterValues();
  m_covMatrix = m_parManager->covarianceMatrix();

  for( int i = 0; i < m_parNames.size(); ++i ){
    
    m_parIndex[m_parNames[i]] = i;
  }

}

void
FitResults::recordFitStats(){
  
  m_eMatrixStatus = m_minManager->eMatrixStatus();
  m_lastCommandStatus = m_minManager->status();
  m_lastCommand = m_minManager->lastCommand();
  m_precision = m_minManager->precision();
  m_strategy = m_minManager->strategy();
  m_estDistToMin = m_minManager->estDistToMinimum();
  m_bestMin = m_minManager->bestMinimum();
}

vector< string >
FitResults::stringSplit(const string& str, const string& delimiters ) const
{
  
  vector< string > tokens;
  
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
  {
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
  
  return tokens;
}


