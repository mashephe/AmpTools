#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "gammaKKExe/constants.h"
#include "gammaKKPlot/gammaKKPlotGenerator.h"
#include "gammaKKAmp/gammaKKHelicityAmp.h"
#include "gammaKKAmp/TwoPiAngles.h"
#include "gammaKKAmp/MultipoleAmps.h"

// Include ROOT 2D histograms currently impossible,
// as gammaKKPlotGenerator::fillProjections
// will return vector<Histogram>, and this will
// be used in PlotGenerator???
// #include "TH2F.h"

gammaKKPlotGenerator::gammaKKPlotGenerator( const ConfigurationInfo* cfgInfo,
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
    gammaKKDataReader dataReader( (**react).data().second );
    m_dataMap[name].loadData( &dataReader );
    
    cout << "Caching accepted MC for reaction:  " << name << endl;
    gammaKKDataReader accMCReader( (**react).accMC().second );
    m_accMCMap[name].loadData( &accMCReader );
    m_accMCMap[name].allocateAmps( ampManager( name ), true );
    ampManager( name ).calcAmplitudes( m_accMCMap[name] );
    
    cout << "Caching generated MC for reaction:  " << name << endl;
    gammaKKDataReader genMCReader( (**react).genMC().second );
    m_genMCMap[name].loadData( &genMCReader );
    m_genMCMap[name].allocateAmps( ampManager( name ), true );
    ampManager( name ).calcAmplitudes( m_genMCMap[name] );

  }
  
  bookHistograms();
  
}


gammaKKPlotGenerator::~gammaKKPlotGenerator(){
  
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
gammaKKPlotGenerator::availablePlots() const { 
    
    return m_histTitles; 
}


void 
gammaKKPlotGenerator::registerPhysics( AmplitudeManager* ampManager ){
    
  ampManager->registerAmplitudeFactor( gammaKKHelicityAmp() );
  ampManager->registerAmplitudeFactor( TwoPiAngles() );
  ampManager->registerAmplitudeFactor( MultipoleAmps() );

}

void
gammaKKPlotGenerator::bookHistograms(){
  
  clearHistograms();
  
  m_histTitles[khm12] = "Mass( 1 2 )";
  m_histTitles[khm13] = "Mass( 1 3 )";
  m_histTitles[khm23] = "Mass( 2 3 )";

  // Photon angles in CM frame
  m_histTitles[khPhotonCosTheta] = "cos#theta_{#gamma}";
  m_histTitles[khPhotonPhi] = "#phi_{#gamma}";
  // 1st Kaon angles in KK CM frame
  m_histTitles[khKaonCosTheta] = "cos#theta_{K_{1}}";
  m_histTitles[khKaonPhi] = "#phi_{K_{1}}";

  m_histVect[khm12]     = Histogram( 60, 0.0, 3.0 ); 
  m_histVect[khm13]     = Histogram( 60, 0.0, 3.0 ); 
  m_histVect[khm23]     = Histogram( 60, 0.0, 3.0 ); 

  m_histVect[khPhotonCosTheta]     = Histogram( 100, -1.0, 1.0 ); 
  m_histVect[khPhotonPhi]          = Histogram( 100, -_PI, _PI ); 

  m_histVect[khKaonCosTheta]     = Histogram( 100, -1.0, 1.0 ); 
  m_histVect[khKaonPhi]          = Histogram( 100, -_PI, _PI ); 

}


vector< Histogram > 
gammaKKPlotGenerator::fillProjections( const string& reactName,
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

      // gamma      
      HepLorentzVector P1 = kinematics->particle(0);
      // Kaon 1
      HepLorentzVector P2 = kinematics->particle(1);
      // Kaon 2
      HepLorentzVector P3 = kinematics->particle(2);
      
      m_histVect[khm12].fill((P1+P2).m(),weight);
      m_histVect[khm13].fill((P1+P3).m(),weight);
      m_histVect[khm23].fill((P2+P3).m(),weight);

      float costhetaPhoton = P1.cosTheta();
      float phiPhoton = P1.phi();

      m_histVect[khPhotonCosTheta].fill(costhetaPhoton,weight);
      m_histVect[khPhotonPhi].fill(phiPhoton,weight);

      // Need to boost to KK rest frame to get angles of Kaon      
      HepLorentzVector resonance = P2 + P3;
      HepLorentzRotation resRestBoost( -resonance.boostVector() );
      HepLorentzVector P2_KKrest = resRestBoost * P2;

      float costhetaKaon = P2_KKrest.cosTheta();
      float phiKaon = P2_KKrest.phi();
      
      m_histVect[khKaonCosTheta].fill(costhetaKaon,weight);
      m_histVect[khKaonPhi].fill(phiKaon,weight);

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
gammaKKPlotGenerator::clearHistograms(){
  
  for( vector< Histogram >::iterator hist = m_histVect.begin();
      hist != m_histVect.end();
      ++hist ){
    
    hist->clear();
  }
}
