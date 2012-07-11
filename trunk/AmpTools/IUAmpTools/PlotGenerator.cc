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
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>
#include <cassert>

#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/NormIntInterface.h"

PlotGenerator::PlotGenerator( const ConfigurationInfo* cfgInfo,
                              const string& parFile ) :
m_cfgInfo( cfgInfo ),
m_fullAmplitudes( 0 ),
m_uniqueAmplitudes( 0 ),
m_emptyHist()
{
  
	ifstream inPar( parFile.c_str() );
	
	unsigned int nProdPars, nAmpPars;
	inPar >> nProdPars >> nAmpPars;
	
	vector< bool > fixedPhase( nProdPars );
	
  // first fill in the vector of production amps
  // this comes mainly from fit parameters, but also constrained
  // amplitudes are given an index to the correct production amplitude

  int ampIndex = 0;
	for( unsigned int i = 0; i < nProdPars; ++i ){
    
		string ampName;
		complex< double > value;
		inPar >> ampName >> value;
    
		const char* lastChar = ampName.substr( ampName.size() - 1, 1 ).c_str();
    
		if( strcmp( lastChar, "+" ) == 0 ){
			
			fixedPhase[i] = true;
			ampName.erase( ampName.size() - 1, 1 );
		}
		else{
			
			fixedPhase[i] = false;
		}
		
		m_ampIndex[ampName] = ampIndex++;
 
    // need to know where to find the error for this amplitude
    m_eMatIndex[ampName] = i;
    
    // keep vector with original fit values
		m_fitProdAmps.push_back( value );
    
    // a parallel vector of zeros
    m_zeroProdAmps.push_back( complex< double >( 0, 0 ) );        
    
    // and a vector from which to distribute values to the
    // amplitude managers (initially this vector is set to 
    // the original fit values)
    m_prodAmps.push_back( value );

    // take apart the amplitude name
    vector<string> ampNameParts = stringSplit(ampName,"::");
    if (ampNameParts.size() != 3){
      cout << "PlotGenerator ERROR: Amplitude name (" << ampName << ")" << endl;
      cout <<  " from the parameter file (" << parFile << ")" << endl;
      cout << "  has the wrong format (must be rct::sum::amp)." << endl;
      assert(false);
    }
    
    // get the AmplitudeInfo
    AmplitudeInfo* ampInfo = m_cfgInfo->amplitude(ampNameParts[0],ampNameParts[1],ampNameParts[2]);

    // for any constrained amplitudes create additional entries in the vector
    // of production amplitudes -- this allows the plotter/user to be able
    // to independently turn off and on single amplitudes (individually), even if these
    // amplitudes were constrained  to other amplitudes in the fit
    if (!ampInfo){
      
      cout << "PlotGenerator ERROR:  requesting non-existent amplitude " << ampName << endl;
      cout << "  from the config file." << endl;
      assert(false);
    }
    if (ampInfo->constraints().size() > 0){
      
      vector<AmplitudeInfo*> constraints = ampInfo->constraints();
      for (unsigned int icst = 0; icst < constraints.size(); icst++){
        
        m_fitProdAmps.push_back( value );
        m_zeroProdAmps.push_back( complex< double >( 0, 0 ) );
        m_prodAmps.push_back( value );
        
        m_ampIndex[constraints[icst]->fullName()] = ampIndex++;

        // this is the location of the constrained amplitude in the error
        // matrix
        m_eMatIndex[constraints[icst]->fullName()] = i;
      }
    }
	}
    
  // read in any floating parameters in the amplitudes
  for( unsigned int i = 0; i < nAmpPars; ++i ){
   
    string parName;
    double value;

    inPar >> parName >> value;
    m_ampParameters[parName] = value;
  }
	
  // now load up the error matrix 
	// the error matrix should have dimension 2 * nProdPars + nAmpPars
	// fill in zeroes for production parameters with fixed phases
	for( unsigned int i = 0; i < ( 2 * nProdPars ) + nAmpPars; ++i ){
		
		m_errorMatrix.push_back( vector< double >( 0 ) );
		
		// read the real row
		for( unsigned int j = 0; j < ( 2 * nProdPars ) + nAmpPars; ++j ){
			      
			if( ( ( j < ( 2 * nProdPars ) ) && fixedPhase[j/2] && ( j % 2 == 1 ) ) ||
          ( ( i < ( 2 * nProdPars ) ) && fixedPhase[i/2] && ( i % 2 == 1 ) ) ){
				
				m_errorMatrix[i].push_back( 0 );
			}
			else{
				
				double value;
				inPar >> value;
				m_errorMatrix[i].push_back( value );
			}
		}		
	}
  
}

PlotGenerator::~PlotGenerator(){
  
  for (map<string, NormIntInterface*>::iterator 
       mapItr = m_normIntMap.begin();
       mapItr != m_normIntMap.end(); ++mapItr){ 
    
    if (mapItr->second) delete mapItr->second;
  }
  
  for( map< string, AmplitudeManager* >::iterator
      mapItr = m_ampManagerMap.begin();
      mapItr != m_ampManagerMap.end();
      ++mapItr ){
    
    if (mapItr->second) delete mapItr->second;
  }        
}

void
PlotGenerator::initialize() {
  
  // build the normalization integrals and amplitude managers
  // for each of the final states
  // this will be useful for generating fit projections

  vector<ReactionInfo*> rctInfoVector = m_cfgInfo->reactionList();
  
  for( unsigned int i = 0; i < rctInfoVector.size(); i++ ){
    
    const ReactionInfo* rctInfo = rctInfoVector[i];
    
    m_normIntMap[rctInfo->reactionName()] = new NormIntInterface( rctInfo->normIntFile() );
    
    vector< string > fsParticles = rctInfo->particleList();
    vector< string > ampNames;
    vector<AmplitudeInfo*> ampInfoVector = m_cfgInfo->amplitudeList(rctInfo->reactionName());
    for (unsigned int j = 0; j < ampInfoVector.size(); j++){
      ampNames.push_back(ampInfoVector[j]->fullName());
    }
    
    // create a amplitude manager for this final state
    AmplitudeManager* ampManager = new AmplitudeManager( fsParticles, rctInfo->reactionName() );
    
    // and hang onto it
    m_ampManagerMap[rctInfo->reactionName()] = ampManager;
    
    // enable the reaction by default
    m_reactEnabled[rctInfo->reactionName()] = true;
    
    // ask derived class plotters to register their 
    // amplitudes and contractions
    registerPhysics( ampManager );
    
    // now put in specific amplitudes
    ampManager->setupFromConfigurationInfo(m_cfgInfo);
    
    // tell the amplitude manager to look at member data of this class for
    // the production amplitudes
    for( vector< string >::iterator ampName = ampNames.begin();
        ampName != ampNames.end();
        ++ampName ){
      
      m_fullAmplitudes.push_back( *ampName );
      
      if( m_ampIndex.find( *ampName ) == m_ampIndex.end() ){
        
        cout << "ERROR:  cannot find production parameter for: " 
        << *ampName << "\n\tAre fit results and config file consistent?"
        << endl;
        
        assert( false );
      }
      
      complex< double >* prodPtr = &(m_prodAmps[m_ampIndex[*ampName]]);
      ampManager->setExternalProductionAmplitude( *ampName, prodPtr );
      
      for( map< string, double >::const_iterator mapItr = m_ampParameters.begin();
           mapItr != m_ampParameters.end();
          ++mapItr ){
       
        ampManager->setAmpParValue( *ampName, mapItr->first, mapItr->second );
      }
    }
  }

  // buildUniqueAmplitudes will also create an initalize the maps that store
  // the enable/disable status of the amplitudes and sums
  buildUniqueAmplitudes();

  recordConfiguration();
}

pair< double, double >
PlotGenerator::intensity( bool accCorrected ) const {
  
  map< string, vector< string > > reactionAmp;
  
  // first loop over amplitudes and separate into reactions
  for( vector< string >::const_iterator amp = m_fullAmplitudes.begin();
      amp != m_fullAmplitudes.end();
      ++amp ){
    
    vector< string > parts = stringSplit( *amp, "::" );
    reactionAmp[parts[0]].push_back( *amp );
  }
  
  pair< double, double > inten( 0, 0 );
  
  // now iterate over reactions and calcualte intensity
  // *** ASSUMING ERRORS ARE 100% CORRELATED *** 
  // *** THIS IS PROBABLY NOT QUITE CORRECT ***
  for( map< string, vector< string > >::iterator mapItr = reactionAmp.begin();
      mapItr != reactionAmp.end();
      ++mapItr ){
    
    // don't add anything if the reaction is disabled
    if( !m_reactEnabled.find( mapItr->first )->second ) continue;
    
    pair< double, double > thisInten = intensity( mapItr->second, accCorrected );
    inten.first += thisInten.first;
    inten.second += thisInten.second;
  }
  
  return inten;
}

pair< double, double >
PlotGenerator::intensity( const vector< string >& amplitudes, bool accCorrected ) const {
  
	// intensity = sum_a sum_a' V_a V*_a' NI( a, a' )
	
	// these have dimension twice that of ampNames since they hold
	// real and imaginary parts independently
	
	// a subset of the larger error matrix
	vector< vector< double > > errorMatrix;
  
	// a vector for the derivatives of the intensity with respect to the
	// real and imaginary parts of the production amplitudes
	vector< double > deriv( 2 * amplitudes.size() );
  
	double intensity = 0;
  
  // quick check to make sure amplitudes all are from the same final state
  // FIX ME -- I CRASH FOR DIFFERENT FINAL STATES
  string fs("");
  for (int i = 0; i < amplitudes.size(); i++){
    
    // take apart the amplitude name
    vector<string> ampNameParts = stringSplit(amplitudes[i],"::");
    if (ampNameParts.size() != 3){
      cout << "PlotGenerator ERROR: Amplitude name (" << amplitudes[i] << ")" << endl;
      cout << "  has the wrong format (must be rct::sum::amp)." << endl;
      assert(false);
    }

    string fscheck = ampNameParts[0];
    if ((i > 0) && (fs != fscheck) && (fscheck != "")){
      
      cout << "PlotGenerator ERROR: mixing final states." << endl;
      assert(false);  
    }
    
    fs = fscheck;
  }

  
	for( vector< string >::const_iterator amp = amplitudes.begin(); 
      amp != amplitudes.end(); ++amp ){
		
    // here and below we need the index of the production parameter
    // in the expanded array (free + constrained)
    // but we want the error matrix index for the original error matrix
    // which is just the free production amplitudes
    
		int ampIndex =  getAmpIndex(*amp);
		int ampEmatIndex = getErrorMatrixIndex(*amp); 
    
		// index errorMatrix and deriv indicies with lower case
		// and their corresponding full indicies in upper case
		int iRe = 2 * ( amp - amplitudes.begin() );
		int iIm = iRe + 1;
		int IRe = 2 * ampEmatIndex;
		int IIm = IRe + 1;		
		
		errorMatrix.push_back( vector< double >( 2 * amplitudes.size() ) );
		errorMatrix.push_back( vector< double >( 2 * amplitudes.size() ) );
    
		deriv[iRe] = 0;
		deriv[iIm] = 0;
    
		for( vector< string >::const_iterator conjAmp = amplitudes.begin(); 
        conjAmp != amplitudes.end(); ++conjAmp ) {
      
      int conjIndex = getAmpIndex(*conjAmp);
			int conjEmatIndex = getErrorMatrixIndex(*conjAmp);
      
			int jRe = 2 * ( conjAmp - amplitudes.begin() );
			int jIm = jRe + 1;
			int JRe = 2 * conjEmatIndex;
			int JIm = JRe + 1;		
            
			complex< double > ampInt;
      if (accCorrected)  ampInt = m_normIntMap.find(fs)->second->ampInt( *amp, *conjAmp );
      else               ampInt = m_normIntMap.find(fs)->second->normInt( *amp, *conjAmp );
      
			errorMatrix[iRe][jRe] = m_errorMatrix[IRe][JRe];
			errorMatrix[iIm][jIm] = m_errorMatrix[IIm][JIm];
      
			deriv[iRe] += 2 * ( real( m_prodAmps[conjIndex] ) * real( ampInt ) +
                     imag( m_prodAmps[conjIndex] ) * imag( ampInt ) );
			deriv[iIm] += 2 * ( imag( m_prodAmps[conjIndex] ) * real( ampInt ) -
                     real( m_prodAmps[conjIndex] ) * imag( ampInt ) );
      
 			intensity += real( m_prodAmps[ampIndex] * conj( m_prodAmps[conjIndex] ) * ampInt );
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

// this computes the phase differences between sums of production
// coefficients -- it is not clear if this has any meaningful signficance

double
PlotGenerator::phaseDiff( const vector< string >& amps1, 
                         const vector< string >& amps2 ){
	
	complex< double > amp1( 0, 0 );	
	for( vector< string >::const_iterator amp = amps1.begin();
      amp != amps1.end(); ++amp ){
		
		amp1 += m_prodAmps[getAmpIndex(*amp)];
	}
	
	complex< double > amp2( 0, 0 );	
	for( vector< string >::const_iterator amp = amps2.begin();
      amp != amps2.end(); ++amp ){
		
		amp2 += m_prodAmps[getAmpIndex(*amp)];
	}
	
	return( phase( amp1 ) - phase( amp2 ) );
}

double
PlotGenerator::phase( const complex< double >& num ) const {
	
	return arg( num );
}

const Histogram& 
PlotGenerator::projection( unsigned int projectionIndex, string reactName,
                          PlotType type ) {
  
  // return an empty histogram if final state is not enabled
  if( !m_reactEnabled[reactName] ) return m_emptyHist;
  
  // use same ampConfig for data since data histograms are
  // independent of which projections are turned on
  string config = ( type == kData ? "" : m_currentConfiguration );
  
  map< string, map< string, vector< Histogram > > >* cachePtr;
  
  switch( type ){
      
    case kData:
      
      cachePtr = &m_dataHistCache;
      break;
      
    case kAccMC:
      
      cachePtr = &m_accMCHistCache;
      break;
      
    case kGenMC:
      
      cachePtr = &m_genMCHistCache;
      break;
      
    default:
      
      assert( false );
  }
  
  map< string, map< string, vector< Histogram > > >::iterator ampCfg =
  cachePtr->find( config );
  
  // short circuit here so second condition doesn't get a null ptr
  if( ( ampCfg == cachePtr->end() ) || 
     ( ampCfg->second.find( reactName ) == ampCfg->second.end() ) ) {
    
    // histograms don't exist -- fill them
    
    (*cachePtr)[config][reactName] =
    fillProjections( reactName, type );
    
    // renormalize MC
    if( type != kData ){
      
      vector< Histogram >* histVect = &((*cachePtr)[config][reactName]);
      for( vector< Histogram >::iterator hist = histVect->begin();
          hist != histVect->end();
          ++hist ){
        
        switch( type ){
            
          case kAccMC: hist->normalize( intensity( false ).first );
            break;
            
          case kGenMC: hist->normalize( intensity( true ).first );
            break;
            
          default:
            break;
        }
      }
    }
  }
  
  // return the correct histogram:
  return (*cachePtr)[config][reactName][projectionIndex];
}

unsigned int
PlotGenerator::getAmpIndex( const string& ampName) const {
  
  map<string, unsigned int>::const_iterator mapItr = m_ampIndex.find(ampName);
  if (mapItr == m_ampIndex.end()){
    cout << "PlotGenerator ERROR:  Could not find amplitude " << ampName << endl;
    assert (false);
  }
  
  return mapItr->second;
}

unsigned int
PlotGenerator::getErrorMatrixIndex( const string& ampName) const {
  
  map<string, unsigned int>::const_iterator mapItr = m_eMatIndex.find(ampName);
  if (mapItr == m_eMatIndex.end()){
    cout << "PlotGenerator ERROR:  Could not find amplitude " << ampName << endl;
    assert (false);
  }
  
  return mapItr->second;
}


bool
PlotGenerator::haveAmp( const string& amp ) const {
  
  return ( m_ampIndex.find( amp ) != m_ampIndex.end() );
}

void 
PlotGenerator::disableAmp( unsigned int uniqueAmpIndex ){
    
  string amp = m_uniqueAmplitudes[uniqueAmpIndex];
    
  for( map< string, unsigned int >::iterator mapItr = m_ampIndex.begin();
      mapItr != m_ampIndex.end();
      ++mapItr ){
    
    vector< string > ampParts = stringSplit( mapItr->first, "::" );
    if( ampParts[2] != amp ) continue;
    
    unsigned int i = mapItr->second;
    
//     cout << "Setting production amp for " << mapItr->first << " from "
//          << m_prodAmps[i] << " to " << m_zeroProdAmps[i] << endl;
    
    m_prodAmps[i] = m_zeroProdAmps[i];
    m_ampEnabled[amp] = false;
  }
  
  recordConfiguration();
}

void
PlotGenerator::enableAmp( unsigned int uniqueAmpIndex ){
      
  string amp = m_uniqueAmplitudes[uniqueAmpIndex];
    
  for( map< string, unsigned int >::iterator mapItr = m_ampIndex.begin();
      mapItr != m_ampIndex.end();
      ++mapItr ){
          
    vector< string > ampParts = stringSplit( mapItr->first, "::" );
    if( ampParts[2] != amp ) continue;
    
    unsigned int i = mapItr->second;

    m_ampEnabled[amp] = true;

    // on turn back on the amplitude in the amplitude manager if the
    // sum to which it belogins is also enabled
      
    if( m_sumEnabled[ampParts[1]] ){
    
     //   cout << "Setting production amp for " << mapItr->first << " from "
     //         << m_prodAmps[i] << " to " << m_fitProdAmps[i] << endl;
    
      m_prodAmps[i] = m_fitProdAmps[i];
    }
  }
  
  recordConfiguration();
}

void 
PlotGenerator::disableSum( unsigned int uniqueSumIndex ){
  
  string sum = m_uniqueSums[uniqueSumIndex];
  
  for( map< string, unsigned int >::iterator mapItr = m_ampIndex.begin();
      mapItr != m_ampIndex.end();
      ++mapItr ){
    
    vector< string > ampParts = stringSplit( mapItr->first, "::" );
    if( ampParts[1] != sum ) continue;
    
    unsigned int i = mapItr->second;

    m_sumEnabled[sum] = false;

    // cout << "Setting production amp for " << mapItr->first << " from "
    //      << m_prodAmps[i] << " to " << m_zeroProdAmps[i] << endl;
    
    m_prodAmps[i] = m_zeroProdAmps[i];
  }
  
  recordConfiguration();
}

void
PlotGenerator::enableSum( unsigned int uniqueSumIndex ){
  
  string sum = m_uniqueSums[uniqueSumIndex];
    
  for( map< string, unsigned int >::iterator mapItr = m_ampIndex.begin();
      mapItr != m_ampIndex.end();
      ++mapItr ){
    
    vector< string > ampParts = stringSplit( mapItr->first, "::" );
    if( ampParts[1] != sum ) continue;
    
    unsigned int i = mapItr->second;
    
    m_sumEnabled[sum] = true;
    
    // only turn back on this part of the sum if the amplitude
    // is also enabled
    
    if( m_ampEnabled[ampParts[2]] ){
    
      // cout << "Setting production amp for " << mapItr->first << " from "
      //      << m_prodAmps[i] << " to " << m_fitProdAmps[i] << endl;
    
      m_prodAmps[i] = m_fitProdAmps[i];
    }
  }
  
  recordConfiguration();
}


bool
PlotGenerator::isAmpEnabled( unsigned int uniqueIndex ) const {
  
  map< string, bool >::const_iterator itr = 
    m_ampEnabled.find( m_uniqueAmplitudes[uniqueIndex] );

  assert( itr != m_ampEnabled.end() );
  
  return itr->second;
}

bool
PlotGenerator::isSumEnabled( unsigned int uniqueIndex ) const {
  
  map< string, bool >::const_iterator itr = 
    m_sumEnabled.find( m_uniqueSums[uniqueIndex] );
  
  assert( itr != m_sumEnabled.end() );
  
  return itr->second;
}

void 
PlotGenerator::disableReaction( const string& reactName ){
  
  m_reactEnabled[reactName] = false;
}

void 
PlotGenerator::enableReaction( const string& reactName ){
  
  m_reactEnabled[reactName] = true;
}

bool
PlotGenerator::isReactionEnabled( const string& reactName ) const {
  
  map< string, bool >::const_iterator mapItr = m_reactEnabled.find( reactName );
  assert( mapItr != m_reactEnabled.end() );
  return mapItr->second;
}

void
PlotGenerator::recordConfiguration() {
  
  // generate a string to identify the current amplitude configuation
  // by concatenating all amplitudes together -- this string will
  // be used as a key in a cache lookup
  
  m_currentConfiguration = "";
  
  for( vector< string >:: iterator ampItr = m_fullAmplitudes.begin();
      ampItr != m_fullAmplitudes.end();
      ++ampItr ){
    
    vector< string > ampParts = stringSplit( *ampItr, "::" );
    
    if( m_ampEnabled[ampParts[2]] && m_sumEnabled[ampParts[1]] ){
      
      m_currentConfiguration += *ampItr;
    }
  }    
}

void
PlotGenerator::buildUniqueAmplitudes(){
  
  m_uniqueAmplitudes.clear();
  m_uniqueSums.clear();
  
  for( vector< string >::iterator amp = m_fullAmplitudes.begin();
      amp != m_fullAmplitudes.end();
      ++amp ){
    
    vector< string > tok = stringSplit( *amp, "::" );
    
    if( find( m_uniqueAmplitudes.begin(), m_uniqueAmplitudes.end(), tok[2] ) ==
       m_uniqueAmplitudes.end() ){
      
      m_uniqueAmplitudes.push_back( tok[2] );
      
      // if this amp hasn't been added to the ampEnabled list, 
      // then add it and set it to true
      if( m_ampEnabled.find( tok[2] ) == m_ampEnabled.end() ) {
        
        m_ampEnabled[tok[2]] = true;
      }
    }
    
    if( find( m_uniqueSums.begin(), m_uniqueSums.end(), tok[1] ) == 
       m_uniqueSums.end() ){
      
      m_uniqueSums.push_back( tok[1] );

      // if this sum hasn't been added to the sumEnabled list, 
      // then add it and set it to true
      if( m_sumEnabled.find( tok[1] ) == m_sumEnabled.end() ){

        m_sumEnabled[tok[1]] = true;
      }
    }
  }
}

const AmplitudeManager& 
PlotGenerator::ampManager( const string& reactName ){
  
  return *(m_ampManagerMap[reactName]);
}


vector< string >
PlotGenerator::reactions() const{

  vector<string> rcts;
  vector<ReactionInfo*> rctInfoVector = m_cfgInfo->reactionList();
  for (unsigned int i = 0; i < rctInfoVector.size(); i++){
    rcts.push_back(rctInfoVector[i]->reactionName());
  }
  return rcts;
}



vector< string >
PlotGenerator::stringSplit(const string& str, const string& delimiters ) const
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

