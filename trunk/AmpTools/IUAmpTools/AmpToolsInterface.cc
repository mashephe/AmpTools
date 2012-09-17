#include "IUAmpTools/AmpToolsInterface.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/AmpToolsInterface.h"


vector<Amplitude*> AmpToolsInterface::m_userAmplitudes;
vector<DataReader*> AmpToolsInterface::m_userDataReaders;


AmpToolsInterface::AmpToolsInterface(ConfigurationInfo* configurationInfo){

  m_ampVecsReactionName = "";

  resetConfigurationInfo(configurationInfo);

}


void
AmpToolsInterface::resetConfigurationInfo(ConfigurationInfo* configurationInfo){

  m_configurationInfo = configurationInfo;

  clear();


    // ************************
    // create a MinuitMinimizationManager
    // ************************

  m_minuitMinimizationManager = new MinuitMinimizationManager(100);
  m_minuitMinimizationManager->setPrecision( 1E-13 );


    // ************************
    // create an AmplitudeManager for each reaction
    // ************************

  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());

    AmplitudeManager* ampMan = new AmplitudeManager(reaction->particleList(),reactionName);
    for (unsigned int i = 0; i < m_userAmplitudes.size(); i++){
      ampMan->registerAmplitudeFactor( *m_userAmplitudes[i] );
    }
    ampMan->setupFromConfigurationInfo( m_configurationInfo );
    m_amplitudeManagers.push_back(ampMan);

  }

    // ************************
    // create a ParameterManager
    // ************************

  m_parameterManager = new ParameterManager ( *m_minuitMinimizationManager, m_amplitudeManagers );
  m_parameterManager->setupFromConfigurationInfo( m_configurationInfo );



    // ************************
    // loop over reactions
    // ************************

  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());
    AmplitudeManager* ampMan = amplitudeManager(reactionName);


      // ************************
      // create DataReaders
      // ************************

    for (unsigned int i = 0; i < m_userDataReaders.size(); i++){
      if (reaction->data().first == m_userDataReaders[i]->name()) 
        m_dataReaderMap[reactionName] 
          = m_userDataReaders[i]->newDataReader(reaction->data().second);
      if (reaction->genMC().first == m_userDataReaders[i]->name()) 
        m_genMCReaderMap[reactionName]
          = m_userDataReaders[i]->newDataReader(reaction->genMC().second);
      if (reaction->accMC().first == m_userDataReaders[i]->name()) 
        m_accMCReaderMap[reactionName]
          = m_userDataReaders[i]->newDataReader(reaction->accMC().second);
    }
    DataReader* dataRdr  =  dataReader(reactionName);
    DataReader* genMCRdr = genMCReader(reactionName);
    DataReader* accMCRdr = accMCReader(reactionName);


      // ************************
      // create a NormIntInterface
      // ************************

    NormIntInterface* normInt = NULL;
    if (genMCRdr && accMCRdr && ampMan){
      normInt = new NormIntInterface(genMCRdr, accMCRdr, *ampMan);
      m_normIntMap[reactionName] = normInt;
      if (reaction->normIntFile() == "")
	cout << "AmpToolsInterface WARNING:  no name given to NormInt file for reaction " 
	 << reactionName << endl;
    }
    else{
      cout << "AmpToolsInterface WARNING:  not creating a NormIntInterface for reaction " 
       << reactionName << endl;
    }



      // ************************
      // create a LikelihoodCalculator
      // ************************

    LikelihoodCalculator* likCalc = NULL;
    if (ampMan && normInt && dataRdr && m_parameterManager){
      likCalc = new LikelihoodCalculator(*ampMan, *normInt, *dataRdr, *m_parameterManager); 
      m_likCalcMap[reactionName] = likCalc;
    }
    else{
      cout << "AmpToolsInterface WARNING:  not creating a LikelihoodCalculator for reaction " 
       << reactionName << endl;
    }

  }

}



double
AmpToolsInterface::likelihood (const string& reactionName){
  LikelihoodCalculator* likCalc = likelihoodCalculator(reactionName);
  if (likCalc) return (*likCalc)();
  return 0.0;
}


double
AmpToolsInterface::likelihood (){
  double L = 0.0;
  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    L += likelihood(reaction->reactionName());
  }
  return L;
}


void
AmpToolsInterface::finalizeFit(){


  MinuitMinimizationManager* fitManager = minuitMinimizationManager();


    // ************************
    // save fit parameters
    // ************************
  
  string outputFile(m_configurationInfo->fitName());  outputFile += ".fit";
  ofstream outFile(outputFile.c_str());
  parameterManager()->writeParameters( outFile );


    // ************************
    // save normalization integrals
    // ************************
  
  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());
    NormIntInterface* normInt = normIntInterface(reactionName);
    normInt->forceCacheUpdate();
    normInt->exportNormIntCache( reaction->normIntFile() );

  }

}



AmplitudeManager*
AmpToolsInterface::amplitudeManager(const string& reactionName){
  for (unsigned int i = 0; i < m_amplitudeManagers.size(); i++){
    if (m_amplitudeManagers[i]->reactionName() == reactionName)
      return m_amplitudeManagers[i];
  }
  cout << "AmpToolsInterface WARNING:  can't find an AmplitudeManager associated with reaction " 
       << reactionName << endl;
  return (AmplitudeManager*) NULL;
}


DataReader*
AmpToolsInterface::dataReader (const string& reactionName){
  if (m_dataReaderMap.find(reactionName) != m_dataReaderMap.end())
    return m_dataReaderMap[reactionName];
  cout << "AmpToolsInterface WARNING:  can't find a DataReader for data associated with reaction " 
       << reactionName << endl;
  return (DataReader*) NULL;
}


DataReader*
AmpToolsInterface::genMCReader (const string& reactionName){
  if (m_genMCReaderMap.find(reactionName) != m_genMCReaderMap.end())
    return m_genMCReaderMap[reactionName];
  cout << "AmpToolsInterface WARNING:  can't find a DataReader for genMC associated with reaction " 
       << reactionName << endl;
  return (DataReader*) NULL;
}


DataReader*
AmpToolsInterface::accMCReader (const string& reactionName){
  if (m_accMCReaderMap.find(reactionName) != m_accMCReaderMap.end())
    return m_accMCReaderMap[reactionName];
  cout << "AmpToolsInterface WARNING:  can't find a DataReader for accMC associated with reaction " 
       << reactionName << endl;
  return (DataReader*) NULL;
}


NormIntInterface*
AmpToolsInterface::normIntInterface (const string& reactionName){
  if (m_normIntMap.find(reactionName) != m_normIntMap.end())
    return m_normIntMap[reactionName];
  cout << "AmpToolsInterface WARNING:  can't find a NormIntInterface associated with reaction " 
       << reactionName << endl;
  return (NormIntInterface*) NULL;
}


LikelihoodCalculator*
AmpToolsInterface::likelihoodCalculator (const string& reactionName){
  if (m_likCalcMap.find(reactionName) != m_likCalcMap.end())
    return m_likCalcMap[reactionName];
  cout << "AmpToolsInterface WARNING:  can't find a LikelihoodCalculator associated with reaction " 
       << reactionName << endl;
  return (LikelihoodCalculator*) NULL;
}



void
AmpToolsInterface::registerAmplitude( const Amplitude& amplitude){

  m_userAmplitudes.push_back(amplitude.clone());

}



void
AmpToolsInterface::registerDataReader( const DataReader& dataReader){

  m_userDataReaders.push_back(dataReader.clone());

}



void
AmpToolsInterface::clear(){

  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());

    if (m_likCalcMap[reactionName]) delete m_likCalcMap[reactionName];
    if (amplitudeManager(reactionName)) delete amplitudeManager(reactionName);
    if (dataReader(reactionName)) delete dataReader(reactionName);
    if (accMCReader(reactionName)) delete accMCReader(reactionName);
    if (genMCReader(reactionName)) delete genMCReader(reactionName);
    if (normIntInterface(reactionName)) delete normIntInterface(reactionName);

  }

  m_amplitudeManagers.clear();
  m_dataReaderMap.clear();
  m_genMCReaderMap.clear();
  m_accMCReaderMap.clear();
  m_normIntMap.clear();
  m_likCalcMap.clear();

  m_ampVecs.deallocAmpVecs();

  if (minuitMinimizationManager()) delete minuitMinimizationManager();
  if (parameterManager()) delete parameterManager();

}


void
AmpToolsInterface::printKinematics(string reactionName, Kinematics* kin){

  ReactionInfo* reaction = m_configurationInfo->reaction(reactionName);
  vector<HepLorentzVector> momenta = kin->particleList();

  if (reaction->particleList().size() != momenta.size()){
    cout << "AmpToolsInterface ERROR:  kinematics incompatible with this reaction" << endl;
    exit(1);
  }

  cout << "  +++++++++++++++++++++++++++++++++" << endl;
  cout << "    EVENT KINEMATICS " << endl;
  for (unsigned int imom = 0; imom < momenta.size(); imom++){
    cout << "      particle " << reaction->particleList()[imom] << endl;
    cout << "          E  = " << momenta[imom].e() << endl;
    cout << "          Px = " << momenta[imom].px() << endl;
    cout << "          Py = " << momenta[imom].py() << endl;
    cout << "          Pz = " << momenta[imom].pz() << endl;
  }
  cout << "  +++++++++++++++++++++++++++++++++" << endl << endl;

}


void
AmpToolsInterface::printAmplitudes(string reactionName, Kinematics* kin){

  AmplitudeManager* ampMan = amplitudeManager(reactionName);
  vector<string> ampNames = ampMan->getAmpNames();
  for (unsigned int iamp = 0; iamp < ampNames.size(); iamp++){
    cout << "    ----------------------------------" << endl;
    cout << "      AMPLITUDE = " << ampNames[iamp] << endl;
    cout << "    ----------------------------------" << endl << endl;
    vector< const Amplitude* > ampFactors = ampMan->getFactors(ampNames[iamp]);
    vector <vector <int> > permutations = ampMan->getPermutations(ampNames[iamp]);
    for (unsigned int iperm = 0; iperm < permutations.size(); iperm++){
      cout << "        PERMUTATION = ";
      for (unsigned int ipar = 0; ipar < permutations[iperm].size(); ipar++){
        cout << permutations[iperm][ipar] << " ";
      }
      cout << endl << endl;
      for (unsigned int ifact = 0; ifact < ampFactors.size(); ifact++){
        cout << "          AMPLITUDE FACTOR = " << ampFactors[ifact]->name() << endl;
        cout << "          RESULT = " 
             << ampFactors[ifact]->calcAmplitude(kin,permutations[iperm]) << endl << endl;
      }
    }
  }

}


void
AmpToolsInterface::printIntensity(string reactionName, Kinematics* kin){

  AmplitudeManager* ampMan = amplitudeManager(reactionName);

  cout << "      ---------------------------------" << endl;
  cout << "        CALCULATING INTENSITY" << endl;
  cout << "      ---------------------------------" << endl << endl;
  double intensity = ampMan->calcIntensity(kin);
  cout << endl << "          INTENSITY = " << intensity << endl << endl << endl;

}


void
AmpToolsInterface::printEventDetails(string reactionName, Kinematics* kin){

  printKinematics(reactionName,kin);
  printAmplitudes(reactionName,kin);
  printIntensity(reactionName,kin);

}


void
AmpToolsInterface::printTestEvents (string reactionName, DataReader* dataReader, int nEvents){


    // ************************
    // find an AmplitudeManager
    // ************************

  AmplitudeManager* ampMan = amplitudeManager(reactionName);


    // ************************
    // print out detailed information for a few events
    // ************************

  cout << "**********************************************" << endl;
  cout << "  AMPLITUDES FOR REACTION " << reactionName << endl;
  cout << "**********************************************" << endl << endl;

  for (unsigned int ievt = 0; ievt < nEvents; ievt++){
    Kinematics* kin = dataReader->getEvent();

    printKinematics(reactionName, kin);
    printAmplitudes(reactionName, kin);
    printIntensity (reactionName, kin);

    delete kin;

  }

}


void
AmpToolsInterface::clearEvents(){

  m_ampVecsReactionName = "";

  m_ampVecs.deallocAmpVecs();

}


void
AmpToolsInterface::loadEvents(DataReader* dataReader){

  clearEvents();
  m_ampVecs.loadData(dataReader);

}


void
AmpToolsInterface::loadEvent(Kinematics* kin, int iEvent, int nEventsTotal){

  m_ampVecs.loadEvent(kin, iEvent, nEventsTotal);

}

double
AmpToolsInterface::processEvents(string reactionName){

  m_ampVecsReactionName = reactionName;

  AmplitudeManager* ampMan = amplitudeManager(m_ampVecsReactionName);

  m_ampVecs.allocateAmps(*ampMan,true);

  return ampMan->calcIntensities(m_ampVecs);

}



double
AmpToolsInterface::intensity(int iEvent){

  if (iEvent >= m_ampVecs.m_iNTrueEvents || iEvent < 0){
    cout << "AmpToolsInterface ERROR:  out of bounds in intensity call" << endl;
    exit(1);
  }

  return m_ampVecs.m_pdIntensity[iEvent];

}


complex<double>
AmpToolsInterface::decayAmplitude (int iEvent, string ampName){

  if (iEvent >= m_ampVecs.m_iNTrueEvents || iEvent < 0){
    cout << "AmpToolsInterface ERROR:  out of bounds in decayAmplitude call" << endl;
    exit(1);
  }

  AmplitudeManager* ampMan = amplitudeManager(m_ampVecsReactionName);

  int iAmp = ampMan->ampIndex(ampName);

  return complex<double>
            (m_ampVecs.m_pdAmps[2*m_ampVecs.m_iNEvents*iAmp+2*iEvent],
             m_ampVecs.m_pdAmps[2*m_ampVecs.m_iNEvents*iAmp+2*iEvent+1]);

}

complex<double>
AmpToolsInterface::productionAmplitude (string ampName){

  return parameterManager()->findParameter(ampName)->value();

}


double
AmpToolsInterface::alternateIntensity(int iEvent){

  double runningIntensity = 0.0;

    // loop over sums

  vector<CoherentSumInfo*> sums = m_configurationInfo->coherentSumList(m_ampVecsReactionName);
  for (unsigned int iSum = 0; iSum < sums.size(); iSum++){

    complex<double> runningAmplitude(0.0,0.0);

    // loop over amps

    vector<AmplitudeInfo*> amps = m_configurationInfo->amplitudeList(m_ampVecsReactionName,sums[iSum]->sumName());
    for (unsigned int iAmp = 0; iAmp < amps.size(); iAmp++){

      complex<double> P = productionAmplitude(amps[iAmp]->fullName());
      complex<double> D = decayAmplitude(iEvent,amps[iAmp]->fullName());

      runningAmplitude += P*D;

    }

    runningIntensity += norm(runningAmplitude);

  }

  return runningIntensity;

}

