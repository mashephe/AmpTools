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
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

PlotGenerator::PlotGenerator( const FitResults& results ) :
m_fitResults( results ),
// in the context of this class we are not going to modify the
// configuration info, but the AmpToolsInterface requires a
// non-const pointer for flexiblity so we will cast away
// the constness although the ConfigurationInfo will really
// remain constant
m_ati( const_cast< ConfigurationInfo* >( results.configInfo() ),
       AmpToolsInterface::kPlotGeneration ),
m_cfgInfo( results.configInfo() ),
m_fullAmplitudes( 0 ),
m_uniqueAmplitudes( 0 ),
m_histVect( 0 ),
m_histTitles( 0 )
{

  vector< string > amps = m_fitResults.ampList();
  
  for( unsigned int i = 0; i < amps.size(); ++i ){
    
    m_ampIndex[amps[i]] = i;
    
    // keep vector with orginal fit values
    m_fitProdAmps.push_back( m_fitResults.productionParameter( amps[i] ) );

    // a parallel vector of zeros
    m_zeroProdAmps.push_back( complex< double >( 0, 0 ) );
    
    // and a vector from which to distribute values to the
    // amplitude managers (initially this vector is set to
    // the original fit values)
    m_prodAmps.push_back( m_fitProdAmps[i] );
  }
  
  // fetch the amplitude parameters
  m_ampParameters = m_fitResults.ampParMap();

  vector<ReactionInfo*> rctInfoVector = m_cfgInfo->reactionList();
  
  for( unsigned int i = 0; i < rctInfoVector.size(); i++ ){
    
    const ReactionInfo* rctInfo = rctInfoVector[i];
    string reactName = rctInfo->reactionName();
    
    m_reactIndex[reactName] = i;
    
    // enable the reaction by default
    m_reactEnabled[rctInfo->reactionName()] = true;
    
    // keep a pointer to the NormalizationIntegralInterface
    m_normIntMap[reactName] = m_ati.normIntInterface( reactName );
    
    // keep a pointer to the AmplitudeManager
    m_intenManagerMap[reactName] = m_ati.intensityManager( reactName );
    
    vector< string > ampNames;
    vector<AmplitudeInfo*> ampInfoVector = m_cfgInfo->amplitudeList( reactName );
    for (unsigned int j = 0; j < ampInfoVector.size(); j++){
      ampNames.push_back(ampInfoVector[j]->fullName());
    }
    
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
      m_intenManagerMap[reactName]->setExternalProductionFactor( *ampName, prodPtr );
      
      for( map< string, double >::const_iterator mapItr = m_ampParameters.begin();
          mapItr != m_ampParameters.end();
          ++mapItr ){
        
        cout << "setting parameter " << mapItr->first << " to " << mapItr->second << endl;
        
        m_intenManagerMap[reactName]->setParValue( *ampName, mapItr->first, mapItr->second );
      }
    }
    
    // now load up the data and MC for that reaction
    m_ati.loadEvents( m_ati.dataReader( reactName ), i * kNumTypes + kData );
    m_ati.loadEvents( m_ati.bkgndReader( reactName ), i * kNumTypes + kBkgnd );
    m_ati.loadEvents( m_ati.accMCReader( reactName ), i * kNumTypes + kAccMC );
    m_ati.loadEvents( m_ati.genMCReader( reactName ), i * kNumTypes + kGenMC );
    
    // calculate the amplitudes and intensities for the accepted and generated MC
    m_ati.processEvents( reactName, i * kNumTypes + kAccMC );
    m_ati.processEvents( reactName, i * kNumTypes + kGenMC );
  }
  
  // buildUniqueAmplitudes will also create an initalize the maps that store
  // the enable/disable status of the amplitudes and sums
  buildUniqueAmplitudes();
  recordConfiguration();
}

/*Empty constructor to contain histograms during event generation*/
PlotGenerator::PlotGenerator( ) :
m_fitResults( *( new FitResults() ) ),
m_histVect( 0 ),
m_histTitles( 0 ),
m_currentEventWeight( 1 )
{ }

/*Delete the histograms in the cache*/
PlotGenerator::~PlotGenerator(){
	map < string, map< string, vector< Histogram* > > >::iterator hit1;
	map< string, vector< Histogram* > >::iterator hit2;
	vector< Histogram* >::iterator hit3;
	
	for( map < string, map< string, vector< Histogram* > > >::iterator hit1 = m_accMCHistCache.begin();hit1 != m_accMCHistCache.end();++hit1 ){
		for(  map< string, vector< Histogram* > >::iterator hit2 = (hit1->second).begin();hit2 != (hit1->second).end();++hit2 ){
			for (vector< Histogram* >::iterator hit3 = (hit2->second).begin();hit3 != (hit2->second).end();++hit3){
				if (*hit3) delete (*hit3);
			}
		}
	}
	for( map < string, map< string, vector< Histogram* > > >::iterator hit1 = m_genMCHistCache.begin();hit1 != m_genMCHistCache.end();++hit1 ){
		for(  map< string, vector< Histogram* > >::iterator hit2 = (hit1->second).begin();hit2 != (hit1->second).end();++hit2 ){
			for (vector< Histogram* >::iterator hit3 = (hit2->second).begin();hit3 != (hit2->second).end();++hit3){
				if (*hit3) delete (*hit3);
			}
		}
	}
	for( map < string, map< string, vector< Histogram* > > >::iterator hit1 = m_dataHistCache.begin();hit1 != m_dataHistCache.end();++hit1 ){
		for(  map< string, vector< Histogram* > >::iterator hit2 = (hit1->second).begin();hit2 != (hit1->second).end();++hit2 ){
			for (vector< Histogram* >::iterator hit3 = (hit2->second).begin();hit3 != (hit2->second).end();++hit3){
				if (*hit3) delete (*hit3);
			}
		}
	}
  for( map < string, map< string, vector< Histogram* > > >::iterator hit1 = m_bkgndHistCache.begin();hit1 != m_bkgndHistCache.end();++hit1 ){
    for(  map< string, vector< Histogram* > >::iterator hit2 = (hit1->second).begin();hit2 != (hit1->second).end();++hit2 ){
      for (vector< Histogram* >::iterator hit3 = (hit2->second).begin();hit3 != (hit2->second).end();++hit3){
        if (*hit3) delete (*hit3);
      }
    }
  }
  
  for( int i=0; i < m_histVect.size(); i++) {

    if( m_histVect[i] ) delete m_histVect[i];
  }
}

pair< double, double >
PlotGenerator::intensity( bool accCorrected ) const {
  
  vector< string > enabledAmps;
  
  // loop over amplitudes and find out what is turned on
  for( vector< string >::const_iterator amp = m_fullAmplitudes.begin();
      amp != m_fullAmplitudes.end();
      ++amp ){
    
    vector< string > parts = stringSplit( *amp, "::" );
    
    // be sure the reaction, sum, and amplitude are enabled
    if( !m_reactEnabled.find( parts[0] )->second ) continue;
    if( !m_sumEnabled.find( parts[1] )->second ) continue;
    if( !m_ampEnabled.find( parts[2] )->second ) continue;
    
    enabledAmps.push_back( *amp );
  }
 
  return m_fitResults.intensity( enabledAmps, accCorrected );
}


Histogram* PlotGenerator::projection( unsigned int projectionIndex, string reactName,
                                     unsigned int type ) {
  
  // return an empty histogram if final state is not enabled
  if( !m_reactEnabled[reactName] ) return NULL;
  
  // use same ampConfig for data since data histograms are
  // independent of which projections are turned on
  string config = ( ( type == kData ) || ( type == kBkgnd ) ?
                   "" : m_currentConfiguration );
  
  map< string, map< string, vector< Histogram* > > > *cachePtr;
  
  switch( type ){
      
    case kData:
      
      cachePtr = &m_dataHistCache;
      break;

    case kBkgnd:
      
      cachePtr = &m_bkgndHistCache;
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
  map< string, map< string, vector< Histogram* > > >::iterator ampCfg = cachePtr->find( config );
  
  // short circuit here so second condition doesn't get a null ptr
  if( ( ampCfg == cachePtr->end() ) || 
      ( ampCfg->second.find( reactName ) == ampCfg->second.end() ) ) {
    
    // histograms don't exist
    // first clear the histogram vector: m_histVect
    clearHistograms();
    
    // this triggers user routines that fill m_histVect
    fillProjections( reactName, type );
    
    //now cache the vector of histograms
    /*Clone m_histVect content and point the cache to the clone address*/
    /*m_histVect is a vector<Histogram*>*/
    m_histVect_clone.clear();
    for( vector< Histogram* >::iterator hist = m_histVect.begin();hist != m_histVect.end();++hist ){
	    m_histVect_clone.push_back((*hist)->Clone());       //The address of the i-th cloned histogram
    }
    (*cachePtr)[config][reactName] = m_histVect_clone;
    
    // renormalize MC
    if( type != kData ){
      
      vector< Histogram* >* histVect = &((*cachePtr)[config][reactName]);
      for( vector< Histogram* >::iterator hist = histVect->begin();
          hist != histVect->end();
          ++hist ){
        
        switch( type ){
            
          case kAccMC: (*hist)->normalize( intensity( false ).first );
            break;
            
          case kGenMC: (*hist)->normalize( intensity( true ).first );
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

Histogram*
PlotGenerator::getHistogram( int index ){
  return (Histogram*)m_histVect[index];
}

void
PlotGenerator::bookHistogram( int index, const string& title, Histogram* hist){
  
  if( index >= m_histVect.size() ){
    
    m_histVect.resize( index + 1 );
    m_histTitles.resize( index + 1 );
  }
  
  m_histVect[index] = hist;
  m_histTitles[index] = title;
}

void
PlotGenerator::bookHistogram( int index, Histogram* hist){

  bookHistogram( index, hist->title(), hist );
}

void
PlotGenerator::clearHistograms(){
  
  for( vector< Histogram*>::iterator hist = m_histVect.begin();
      hist != m_histVect.end();
      ++hist ){
    
    (*hist)->clear();
  }
}

void
PlotGenerator::fillProjections( const string& reactName, unsigned int type ){
  
  bool isDataOrBkgnd = ( ( type == kData ) || ( type == kBkgnd ) ? true : false );
  int dataIndex = m_reactIndex[reactName] * kNumTypes + type;
      
  // calculate intensities for MC:
  if( !isDataOrBkgnd ) m_ati.processEvents( reactName, dataIndex );
        
  // loop over ampVecs and fill histograms
  for( unsigned int i = 0; i < m_ati.numEvents( dataIndex ); ++i ){
    
    // the subsequent calls here allocate new memory for
    // the Kinematics object
    Kinematics* kin = m_ati.kinematics(i, dataIndex);
    m_currentEventWeight = kin->weight();
    
    if( !isDataOrBkgnd ){

      // m_ati.intensity already contains a possible MC-event weight
      m_currentEventWeight = m_ati.intensity( i, dataIndex );
    }
       
    // the user defines this function in the derived class and it
    // calls the fillHistogram method immediately below
    projectEvent( kin );
    
    // cleanup allocated memory
    delete kin;
  }
}

void
PlotGenerator::fillHistogram( int histIndex, double valueX ){
	vector < double > tmp;
	tmp.push_back(valueX);
	m_histVect[histIndex]->fill( tmp, m_currentEventWeight );
}

void
PlotGenerator::fillHistogram( int histIndex, double valueX,double valueY ){
	vector < double > tmp;
	tmp.push_back(valueX);
	tmp.push_back(valueY);
	m_histVect[histIndex]->fill(tmp, m_currentEventWeight );
}

void
PlotGenerator::fillHistogram( int histIndex, vector <double> &data, double weight){
  m_histVect[histIndex]->fill(data, m_currentEventWeight*weight );
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

