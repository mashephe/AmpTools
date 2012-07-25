
#include "DalitzPlot/DalitzPlotGenerator.h"
#include "DalitzAmp/BreitWigner.h"

DalitzPlotGenerator::DalitzPlotGenerator( const ConfigurationInfo* cfgInfo,
                                          const string& parFile ) :
PlotGenerator( cfgInfo, parFile ),
m_histTitles( kNumHist ),
m_histVect( kNumHist )
{
    
  initialize();

  vector< ReactionInfo* > reactions = cfgInfo->reactionList();
  
  for( vector< ReactionInfo* >::iterator react = reactions.begin();
      react != reactions.end();
      ++react ){
    
    string name = (**react).reactionName();
    
    cout << "Caching data for reaction:  " << name << endl;
    DalitzDataReader dataReader( (**react).data().second );
    m_dataMap[name].loadData( &dataReader );
    
    cout << "Caching accepted MC for reaction:  " << name << endl;
    DalitzDataReader accMCReader( (**react).accMC().second );
    m_accMCMap[name].loadData( &accMCReader );
    m_accMCMap[name].allocateAmps( ampManager( name ), true );
    ampManager( name ).calcAmplitudes( m_accMCMap[name] );
    
    cout << "Caching generated MC for reaction:  " << name << endl;
    DalitzDataReader genMCReader( (**react).genMC().second );
    m_genMCMap[name].loadData( &genMCReader );
    m_genMCMap[name].allocateAmps( ampManager( name ), true );
    ampManager( name ).calcAmplitudes( m_genMCMap[name] );

  }
  
  bookHistograms();
  
}


DalitzPlotGenerator::~DalitzPlotGenerator(){
  
  // flush the data cache
  
  for( map< string, AmpVecs >::iterator mapItr = 
      m_dataMap.begin();
      mapItr != m_dataMap.end();
      ++mapItr ){
    
    mapItr->second.deallocAmpVecs();
  }
  
  for( map< string, AmpVecs >::iterator mapItr = 
      m_accMCMap.begin();
      mapItr != m_accMCMap.end();
      ++mapItr ){
    
    mapItr->second.deallocAmpVecs();
  }
  
  for( map< string, AmpVecs >::iterator mapItr = 
      m_genMCMap.begin();
      mapItr != m_genMCMap.end();
      ++mapItr ){
    
    mapItr->second.deallocAmpVecs();
  }    
}


const vector< string >& 
DalitzPlotGenerator::availablePlots() const { 
    
    return m_histTitles; 
}


void 
DalitzPlotGenerator::registerPhysics( AmplitudeManager* ampManager ){
    
  ampManager->registerAmplitudeFactor( BreitWigner() );

}

void
DalitzPlotGenerator::bookHistograms(){
  
  clearHistograms();
  
  m_histTitles[khm12] = "Mass( 1 2 )";
  m_histTitles[khm13] = "Mass( 1 3 )";
  m_histTitles[khm23] = "Mass( 2 3 )";

  m_histVect[khm12]     = Histogram( 60, 0.0, 3.0 ); 
  m_histVect[khm13]     = Histogram( 60, 0.0, 3.0 ); 
  m_histVect[khm23]     = Histogram( 60, 0.0, 3.0 ); 

}


vector< Histogram > 
DalitzPlotGenerator::fillProjections( const string& reactName,
                                     PlotType type ){
    
  bool isData;
  AmpVecs* ampVecsPtr;
  
  switch( type ){
      
    case kAccMC:
      
      assert( m_accMCMap.find( reactName ) != m_accMCMap.end() );
      ampVecsPtr = &(m_accMCMap[reactName]);
      isData = false;
      break;
      
    case kGenMC:
      
      assert( m_genMCMap.find( reactName ) != m_genMCMap.end() );
      ampVecsPtr= &(m_genMCMap[reactName]);
      isData = false;
      break;
      
    case kData:
      
      assert( m_dataMap.find( reactName ) != m_dataMap.end() );
      ampVecsPtr = &(m_dataMap[reactName]);
      isData = true;
      break;
  }
  
  clearHistograms();
  
  // calculate intensities for MC:
  if( !isData )ampManager( reactName ).calcIntensities( *ampVecsPtr, false );
  
  // loop over ampVecs and fill histograms
  for( unsigned int i = 0; i < ampVecsPtr->m_iNTrueEvents; ++i ){
    
    try{
      
      const Kinematics* kinematics = ampVecsPtr->getEvent( i );
      
      float weight = ( isData ? 1.0 : ampVecsPtr->m_pdIntensity[i] );
      weight *= kinematics->weight();
      
      HepLorentzVector P1 = kinematics->particle(0);
      HepLorentzVector P2 = kinematics->particle(1);
      HepLorentzVector P3 = kinematics->particle(2);
      
      m_histVect[khm12].fill((P1+P2).m(),weight);
      m_histVect[khm13].fill((P1+P3).m(),weight);
      m_histVect[khm23].fill((P2+P3).m(),weight);
      
      delete kinematics;
    }
    catch( ... ){
      
      cout << "Exception thrown in projection plotting...skipping event." << endl;
      continue;
    }
  }
  
  return m_histVect;
}

void
DalitzPlotGenerator::clearHistograms(){
  
  for( vector< Histogram >::iterator hist = m_histVect.begin();
      hist != m_histVect.end();
      ++hist ){
    
    hist->clear();
  }
}
