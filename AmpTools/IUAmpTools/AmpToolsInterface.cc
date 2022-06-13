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


#include "IUAmpTools/AmpToolsInterface.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/AmpVecs.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "IUAmpTools/report.h"

static const char* kModule = "AmpToolsInterface";

vector<Amplitude*> AmpToolsInterface::m_userAmplitudes;
vector<DataReader*> AmpToolsInterface::m_userDataReaders;
unsigned int AmpToolsInterface::m_randomSeed = 0;

AmpToolsInterface::AmpToolsInterface( FunctionalityFlag flag ) :
m_functionality( flag ),
m_configurationInfo( NULL ),
m_minuitMinimizationManager(NULL),
m_parameterManager(NULL),
m_fitResults(NULL)
{
  srand( m_randomSeed );
}


AmpToolsInterface::AmpToolsInterface(ConfigurationInfo* configurationInfo, FunctionalityFlag flag ):
m_functionality( flag ),
m_configurationInfo(configurationInfo),
m_minuitMinimizationManager(NULL),
m_parameterManager(NULL),
m_fitResults(NULL){
  
  resetConfigurationInfo(configurationInfo);
  srand( m_randomSeed );
}


void
AmpToolsInterface::resetConfigurationInfo(ConfigurationInfo* configurationInfo){
  
  m_configurationInfo = configurationInfo;
  
  clear();
  
  report( DEBUG, kModule ) << "resetting from configuration info" << endl;
  
  if( m_functionality == kFull ){
    
    // ************************
    // create a MinuitMinimizationManager
    // ************************
    
    m_minuitMinimizationManager = new MinuitMinimizationManager(500);
  }
  
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
    
    if( m_functionality == kFull ){
      ampMan->setOptimizeParIteration( true );
      ampMan->setFlushFourVecsIfPossible( true );
    }
    
    if( m_functionality == kMCGeneration ){
      ampMan->setOptimizeParIteration( false );
      ampMan->setFlushFourVecsIfPossible( false );
      ampMan->setForceUserVarRecalculation( true );
    }
    
    m_intensityManagers.push_back(ampMan);
  }
  
  if( m_functionality == kFull ){
    
    // ************************
    // create a ParameterManager
    // ************************
    
    m_parameterManager = new ParameterManager ( m_minuitMinimizationManager, m_intensityManagers );
    m_parameterManager->setupFromConfigurationInfo( m_configurationInfo );
    
  }
  
  // ************************
  // loop over reactions
  // ************************
  
  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
    
    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());
    IntensityManager* intenMan = intensityManager(reactionName);
    
    if (!intenMan)
      report( WARNING, kModule ) << "not creating an AmplitudeManager for reaction " << reactionName << endl;
    
    
    if( m_functionality == kFull || m_functionality == kPlotGeneration ){
      
      // ************************
      // create DataReaders
      // ************************
      
      for (unsigned int i = 0; i < m_userDataReaders.size(); i++){
        if (reaction->data().first == m_userDataReaders[i]->name())
          m_dataReaderMap[reactionName]
	  = m_userDataReaders[i]->newDataReader(reaction->data().second);
        if (reaction->bkgnd().first == m_userDataReaders[i]->name())
          m_bkgndReaderMap[reactionName]
          = m_userDataReaders[i]->newDataReader(reaction->bkgnd().second);
        if (reaction->genMC().first == m_userDataReaders[i]->name())
          m_genMCReaderMap[reactionName]
          = m_userDataReaders[i]->newDataReader(reaction->genMC().second);
        if (reaction->accMC().first == m_userDataReaders[i]->name())
          m_accMCReaderMap[reactionName]
          = m_userDataReaders[i]->newDataReader(reaction->accMC().second);
      }
      DataReader* dataRdr  =  dataReader(reactionName);
      DataReader* bkgndRdr = bkgndReader(reactionName);
      DataReader* genMCRdr = genMCReader(reactionName);
      DataReader* accMCRdr = accMCReader(reactionName);
      
      m_uniqueDataSets.insert( dataRdr );
      m_uniqueDataSets.insert( bkgndRdr );
      m_uniqueDataSets.insert( genMCRdr );
      m_uniqueDataSets.insert( accMCRdr );
      
      if (!dataRdr)
        report( WARNING, kModule ) << "not creating a DataReader for data associated with reaction " << reactionName << endl;
      if (!genMCRdr)
        report( WARNING, kModule ) << "not creating a DataReader for generated MC associated with reaction " << reactionName << endl;
      if (!accMCRdr)
        report( WARNING, kModule ) << "not creating a DataReader for accepted MC associated with reaction " << reactionName << endl;
      
      if( m_functionality == kFull ){
 
        // ************************
        // create a NormIntInterface
        // ************************
        
        // note that in the case that the ATI is being used for plot generation
        // then the NI's should be obtained from the FitResults object as this
        // contains the NI cache at the end of the fit
        
        NormIntInterface* normInt = NULL;
        if (genMCRdr && accMCRdr && intenMan && !(reaction->normIntFileInput())){
          
          normInt = new NormIntInterface(genMCRdr, accMCRdr, *intenMan);
          m_normIntMap[reactionName] = normInt;
          if (reaction->normIntFile() == "")
            report( WARNING, kModule ) << "no name given to NormInt file for reaction " << reactionName << endl;
        }
        else if (reaction->normIntFileInput()){

          normInt = new NormIntInterface(reaction->normIntFile());
          m_normIntMap[reactionName] = normInt;
        }
        else{

          report( WARNING, kModule ) << "not creating a NormIntInterface for reaction " << reactionName << endl;
        }
        
        // ************************
        // create a LikelihoodCalculator
        // ************************
        
        LikelihoodCalculator* likCalc = NULL;
        if (intenMan && normInt && dataRdr && m_parameterManager){
          likCalc = new LikelihoodCalculator(*intenMan, *normInt, dataRdr, bkgndRdr, *m_parameterManager);
          m_likCalcMap[reactionName] = likCalc;
        }
        else{
          report( WARNING, kModule ) << "not creating a LikelihoodCalculator for reaction " << reactionName << endl;
        }
      }
    }
  }
  
  // ************************
  // create FitResults
  // ************************
    
  if( m_functionality == kFull ){
    
    m_fitResults = new FitResults( m_configurationInfo,
                                  m_intensityManagers,
                                  m_likCalcMap,
                                  m_normIntMap,
                                  m_minuitMinimizationManager,
                                  m_parameterManager );    
  }
  
  // if functionality is for PlotGeneration then fitResults
  // is left NULL -- in general m_fitResults is used to keep track
  // of an output of the interface
}



double
AmpToolsInterface::likelihood (const string& reactionName) const {
  LikelihoodCalculator* likCalc = likelihoodCalculator(reactionName);
  if (likCalc) return (*likCalc)();
  return 0.0;
}


double
AmpToolsInterface::likelihood () const {
  double L = 0.0;
  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    L += likelihood(reaction->reactionName());
  }
  return L;
}

void
AmpToolsInterface::reinitializePars(){
  
  // shouldn't be callin' unless you're fittin'
  if( m_functionality != kFull ) return;

  // reinitialize production parameters
  vector< AmplitudeInfo* > amps = m_configurationInfo->amplitudeList();
  for( vector< AmplitudeInfo* >::iterator ampItr = amps.begin();
      ampItr != amps.end();
      ++ampItr ){
   
    if( (**ampItr).fixed() ) continue;
    
    string ampName = (**ampItr).fullName();
    complex< double > prodPar = (**ampItr).value();

    m_parameterManager->setProductionParameter( ampName, prodPar );
  }

  // reinitialize amplitude parameters
  vector<ParameterInfo*> parInfoVec = m_configurationInfo->parameterList();
  for ( vector< ParameterInfo* >::iterator parItr = parInfoVec.begin(); parItr != parInfoVec.end(); ++parItr ){

    if( (**parItr).fixed() ) continue;
 
    string parName = (**parItr).parName();
    double value = (**parItr).value();

    m_parameterManager->setAmpParameter( parName, value );
  }

}

void
AmpToolsInterface::randomizeProductionPars( float maxFitFraction ){
  
  // shouldn't be callin' unless you're fittin'
  if( m_functionality != kFull ) return;
  
  vector< AmplitudeInfo* > amps = m_configurationInfo->amplitudeList();
    
  for( vector< AmplitudeInfo* >::iterator ampItr = amps.begin();
      ampItr != amps.end();
      ++ampItr ){
   
    if( (**ampItr).fixed() ) continue;
    
    string ampName = (**ampItr).fullName();
    string reac = (**ampItr).reactionName();
  
    double numSignalEvents = likelihoodCalculator(reac)->numSignalEvents();
    // for the NI interface to use the cache to avoid a mess with MPI jobs and
    // retrigging of NI calculation -- we don't require a precise number
    // for this anyway!
    double normInt = normIntInterface(reac)->normInt(ampName,ampName,true).real();
    double ampScale = intensityManager(reac)->getScale(ampName);
    
    // fit fraction = ( scale^2 * prodPar^2 * normInt ) / numSignalEvents
    
    double fitFraction = random( maxFitFraction );
    double prodParMag = sqrt( ( numSignalEvents * fitFraction ) /
                              ( ampScale * ampScale * normInt ) );
    double prodParPhase = ( (**ampItr).real() ? 0 : random( 2*PI ) );
    
    complex< double > prodPar( prodParMag*cos( prodParPhase ),
                              prodParMag*sin( prodParPhase ) );
    
    m_parameterManager->setProductionParameter( ampName, prodPar );
  }
}

void
AmpToolsInterface::randomizeParameter( const string& parName, float min, float max ){

  vector<ParameterInfo*> parInfoVec = m_configurationInfo->parameterList();

  auto parItr = parInfoVec.begin();
  for( ; parItr != parInfoVec.end(); ++parItr ){
    
    if( (**parItr).parName() == parName ) break;
  }
  
  if( parItr == parInfoVec.end() ){
    
    report( ERROR, kModule ) << "request to randomize nonexistent parameter:  " << parName << endl;
    return;
  }
  
  if( (**parItr).fixed() ){
    
    report( ERROR, kModule ) << "request to randomize a parameter named " << parName
         << " that is fixed.  Ignoring this request." << endl;
    return;
  }
  
  double value = min + random( max - min );
  m_parameterManager->setAmpParameter( parName, value );
}

void
AmpToolsInterface::finalizeFit( const string& tag ){
  
  // ************************
  // save fit parameters
  // ************************
  
  m_fitResults->saveResults();
  m_fitResults->writeResults( m_configurationInfo->fitOutputFileName( tag ) );
  
  // ************************
  // save normalization integrals
  // ************************
  
  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
    
    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    
    if( !reaction->normIntFileInput() ){
      string reactionName(reaction->reactionName());
      NormIntInterface* normInt = normIntInterface(reactionName);
      // the call to FitResults::writeResults will force a cache update
      // there is no need to do it twice
      //      if (normInt->hasAccessToMC()) normInt->forceCacheUpdate();
      normInt->exportNormIntCache( reaction->normIntFile() );
    }
  }
}

IntensityManager*
AmpToolsInterface::intensityManager(const string& reactionName) const {
  for (unsigned int i = 0; i < m_intensityManagers.size(); i++){
    if (m_intensityManagers[i]->reactionName() == reactionName)
      return m_intensityManagers[i];
  }
  return (IntensityManager*) NULL;
}


DataReader*
AmpToolsInterface::dataReader (const string& reactionName) const {
  if (m_dataReaderMap.find(reactionName) != m_dataReaderMap.end())
    return m_dataReaderMap.find(reactionName)->second;
  return (DataReader*) NULL;
}

DataReader*
AmpToolsInterface::bkgndReader (const string& reactionName) const {
  if (m_bkgndReaderMap.find(reactionName) != m_bkgndReaderMap.end())
    return m_bkgndReaderMap.find(reactionName)->second;
  return (DataReader*) NULL;
}

DataReader*
AmpToolsInterface::genMCReader (const string& reactionName) const {
  if (m_genMCReaderMap.find(reactionName) != m_genMCReaderMap.end())
    return m_genMCReaderMap.find(reactionName)->second;
  return (DataReader*) NULL;
}


DataReader*
AmpToolsInterface::accMCReader (const string& reactionName) const {
  if (m_accMCReaderMap.find(reactionName) != m_accMCReaderMap.end())
    return m_accMCReaderMap.find(reactionName)->second;
  return (DataReader*) NULL;
}


NormIntInterface*
AmpToolsInterface::normIntInterface (const string& reactionName) const {
  if (m_normIntMap.find(reactionName) != m_normIntMap.end())
    return m_normIntMap.find(reactionName)->second;
  return (NormIntInterface*) NULL;
}


LikelihoodCalculator*
AmpToolsInterface::likelihoodCalculator (const string& reactionName) const {
  if (m_likCalcMap.find(reactionName) != m_likCalcMap.end())
    return m_likCalcMap.find(reactionName)->second;
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
  
  if( m_configurationInfo != NULL ){
    
    for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
      
      ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
      string reactionName(reaction->reactionName());
      
      if (likelihoodCalculator(reactionName)) delete m_likCalcMap[reactionName];
      if (normIntInterface(reactionName)) delete normIntInterface(reactionName);
    }

    // logic above won't work for these since the lookup function for each reaction
    // depends on the reaction itself
    for (unsigned int i = 0; i < m_intensityManagers.size(); i++){
      delete m_intensityManagers[i];
    }

    for( auto dataReader = m_uniqueDataSets.begin();
        dataReader != m_uniqueDataSets.end(); ++dataReader ){
      
      if( *dataReader ) delete *dataReader;
    }
  }
  
  
  m_intensityManagers.clear();
  m_dataReaderMap.clear();
  m_genMCReaderMap.clear();
  m_accMCReaderMap.clear();
  m_bkgndReaderMap.clear();
  m_uniqueDataSets.clear();
  m_normIntMap.clear();
  m_likCalcMap.clear();
  
  for (unsigned int i = 0; i < MAXAMPVECS; i++){
    m_ampVecs[i].deallocAmpVecs();
    m_ampVecsReactionName[i] = "";
  }
  
  if (minuitMinimizationManager()) delete minuitMinimizationManager();
  if (parameterManager()) delete parameterManager();
  if (fitResults()) delete fitResults();
}



void
AmpToolsInterface::clearEvents(unsigned int iDataSet){
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  m_ampVecsReactionName[iDataSet] = "";
  m_ampVecs[iDataSet].deallocAmpVecs();
  
}


void
AmpToolsInterface::loadEvents(DataReader* dataReader,
                              unsigned int iDataSet){
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  clearEvents(iDataSet);
  
  // if the data reader is null then this will just leave the AmpVecs
  // object in a state where it contains no data, which is consistent
  // with no data reader available (some DataReaders, like those for
  // background) are optional
  if( dataReader != NULL ) m_ampVecs[iDataSet].loadData(dataReader);
}


void
AmpToolsInterface::loadEvent(Kinematics* kin, int iEvent, int nEventsTotal,
                             unsigned int iDataSet){
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  m_ampVecs[iDataSet].loadEvent(kin, iEvent, nEventsTotal);
  
}


double
AmpToolsInterface::processEvents(string reactionName,
                                 unsigned int iDataSet) {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  bool isFirstPass = (m_ampVecs[iDataSet].m_pdAmps == 0);
  
  m_ampVecsReactionName[iDataSet] = reactionName;
  
  IntensityManager* intenMan = intensityManager(reactionName);
  
  if (isFirstPass) m_ampVecs[iDataSet].allocateTerms(*intenMan,true);
  
  return intenMan->calcIntensities(m_ampVecs[iDataSet]);
  
}


int
AmpToolsInterface::numEvents(unsigned int iDataSet) const {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  return m_ampVecs[iDataSet].m_iNTrueEvents;
  
}

double
AmpToolsInterface::sumWeights(unsigned int iDataSet) const {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  return m_ampVecs[iDataSet].m_dSumWeights;
  
}


Kinematics*
AmpToolsInterface::kinematics(int iEvent,
                              unsigned int iDataSet){
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  return m_ampVecs[iDataSet].getEvent(iEvent);
  
}


double
AmpToolsInterface::intensity(int iEvent,
                             unsigned int iDataSet) const {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  if (iEvent >= m_ampVecs[iDataSet].m_iNTrueEvents || iEvent < 0){
    report( ERROR, kModule ) << "out of bounds in intensity call" << endl;
    assert(false);
  }
  
  return m_ampVecs[iDataSet].m_pdIntensity[iEvent];
  
}


complex<double>
AmpToolsInterface::decayAmplitude (int iEvent, string ampName,
                                   unsigned int iDataSet) const {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  if (iEvent >= m_ampVecs[iDataSet].m_iNTrueEvents || iEvent < 0){
    report( ERROR, kModule ) << "out of bounds in decayAmplitude call" << endl;
    assert(false);
  }
  
  IntensityManager* intenMan = intensityManager(m_ampVecsReactionName[iDataSet]);
  
  int iAmp = intenMan->termIndex(ampName);
  
  // fix!! this experession is not generally correct for all intensity
  // managers -- put as helper function in AmpVecs?
  
  return complex<double>
  (m_ampVecs[iDataSet].m_pdAmps[2*m_ampVecs[iDataSet].m_iNEvents*iAmp+2*iEvent],
   m_ampVecs[iDataSet].m_pdAmps[2*m_ampVecs[iDataSet].m_iNEvents*iAmp+2*iEvent+1]);
  
}

complex<double>
AmpToolsInterface::scaledProductionAmplitude (string ampName, unsigned int iDataSet) const {
  
  const IntensityManager* intenMan = intensityManager(m_ampVecsReactionName[iDataSet]);
  
  double scale = intenMan->getScale( ampName );
  complex< double > prodAmp = intenMan->productionFactor( ampName );
  
  return scale * prodAmp;
}


double
AmpToolsInterface::alternateIntensity(int iEvent,
                                      unsigned int iDataSet) const {
  
  if (iDataSet >= MAXAMPVECS){
    report( ERROR, kModule ) << "data set index out of range" << endl;
    assert(false);
  }
  
  double runningIntensity = 0.0;
  
  // loop over sums
  
  vector<CoherentSumInfo*> sums = m_configurationInfo->coherentSumList(m_ampVecsReactionName[iDataSet]);
  for (unsigned int iSum = 0; iSum < sums.size(); iSum++){
    
    complex<double> runningAmplitude(0.0,0.0);
    
    // loop over amps
    
    vector<AmplitudeInfo*> amps =
    m_configurationInfo->amplitudeList(m_ampVecsReactionName[iDataSet],sums[iSum]->sumName());
    for (unsigned int iAmp = 0; iAmp < amps.size(); iAmp++){
      
      complex<double> P = scaledProductionAmplitude(amps[iAmp]->fullName());
      complex<double> D = decayAmplitude(iEvent,amps[iAmp]->fullName(),iDataSet);
      
      runningAmplitude += P*D;
      
    }
    
    runningIntensity += norm(runningAmplitude);
    
  }
  
  return runningIntensity;
}

void
AmpToolsInterface::printKinematics(string reactionName, Kinematics* kin) const {
  
  ReactionInfo* reaction = m_configurationInfo->reaction(reactionName);
  vector<TLorentzVector> momenta = kin->particleList();
  
  if (reaction->particleList().size() != momenta.size()){
    report( ERROR, kModule ) << "kinematics incompatible with this reaction" << endl;
    assert(false);
  }
  
  report( INFO, kModule ) << "  +++++++++++++++++++++++++++++++++" << endl;
  report( INFO, kModule ) << "    EVENT KINEMATICS " << endl;
  streamsize defaultStreamSize = cout.precision(15);
  for (unsigned int imom = 0; imom < momenta.size(); imom++){
    report( INFO, kModule ) << "      particle " << reaction->particleList()[imom] << endl;
    report( INFO, kModule ) << "          E  = " << momenta[imom].E() << endl;
    report( INFO, kModule ) << "          Px = " << momenta[imom].Px() << endl;
    report( INFO, kModule ) << "          Py = " << momenta[imom].Py() << endl;
    report( INFO, kModule ) << "          Pz = " << momenta[imom].Pz() << endl;
  }
  cout.precision(defaultStreamSize);
  report( INFO, kModule ) << "  +++++++++++++++++++++++++++++++++" << endl << endl;
}


void
AmpToolsInterface::printAmplitudes(string reactionName, Kinematics* kin) const {
  
  IntensityManager* intenMan = intensityManager(reactionName);
  
  if( intenMan->type() != IntensityManager::kAmplitude ){
    
    report( NOTICE, kModule ) << "printAmplitudes is being called for a reaction "
                              << "       that is not setup for an amplitude fit."
                              << "       (Nothing more will be printed.)" << endl;
    
    return;
  }
  
  AmplitudeManager* ampMan = dynamic_cast< AmplitudeManager* >( intenMan );
  
  vector<string> ampNames = ampMan->getTermNames();
  
  // we need to use the AmplitudeManager for this call in order to
  // exercise the GPU code for the amplitude calculation
  
  AmpVecs aVecs;
  aVecs.loadEvent(kin);
  aVecs.allocateTerms(*intenMan,true);
  
  ampMan->calcTerms(aVecs);
  
#ifdef GPU_ACCELERATION
  aVecs.allocateCPUAmpStorage( *intenMan );
  aVecs.m_gpuMan.copyAmpsFromGPU( aVecs );
#endif
  
  int nAmps = ampNames.size();
  for (unsigned int iamp = 0; iamp < nAmps; iamp++){
    
    report( INFO, kModule ) << "    ----------------------------------" << endl;
    report( INFO, kModule ) << "      AMPLITUDE = " << ampNames[iamp] << endl;
    report( INFO, kModule ) << "    ----------------------------------" << endl << endl;
    
    vector< const Amplitude* > ampFactors = ampMan->getFactors(ampNames[iamp]);
    vector <vector <int> > permutations = ampMan->getPermutations(ampNames[iamp]);
    
    report( INFO, kModule ) << "      PRODUCT OF FACTORS" << endl
                            << "      SUMMED OVER PERMUTATIONS = ( "
                            << aVecs.m_pdAmps[iamp*2] << ", "
                            << aVecs.m_pdAmps[iamp*2+1] << " )" << endl << endl;
    
    int nPerm = permutations.size();
    
    if( iamp == nAmps-1 ){
      
      // for the last amplitude, the pdAmpFactors array will still hold
      // the data for all of the factors and permutations of the amplitude
      // so go ahead and print those to the screen for the user
      
      for (unsigned int iperm = 0; iperm < nPerm; iperm++){
        
        report( INFO, kModule ) << "        PERMUTATION = ";
        for (unsigned int ipar = 0; ipar < permutations[iperm].size(); ipar++){
          report( INFO, kModule ) << permutations[iperm][ipar] << " ";
        }
        
        report( INFO, kModule ) << endl << endl;
        
        int nFact = ampFactors.size();
        
        for (unsigned int ifact = 0; ifact < nFact; ifact++){
          
          report( INFO, kModule ) << "          AMPLITUDE FACTOR = " << ampFactors[ifact]->name() << endl;
          report( INFO, kModule ) << "          IDENTIFIER = " << ampFactors[ifact]->identifier() << endl;
          report( INFO, kModule ) << "          RESULT = ( "
                                  << aVecs.m_pdAmpFactors[ifact*nPerm*2+iperm*2] << ", "
                                  << aVecs.m_pdAmpFactors[ifact*nPerm*2+iperm*2+1] << " )"
                                  << endl << endl;
        }
      }
    }
  }
  
  // Deallocate memory and return
  aVecs.deallocAmpVecs();
}


void
AmpToolsInterface::printIntensity(string reactionName, Kinematics* kin) const {
  
  IntensityManager* intenMan = intensityManager(reactionName);
  
  report( INFO, kModule ) << "      ---------------------------------" << endl;
  report( INFO, kModule ) << "        CALCULATING INTENSITY" << endl;
  report( INFO, kModule ) << "      ---------------------------------" << endl << endl;
  double intensity = intenMan->calcIntensity(kin);
  report( INFO, kModule ) << endl << "          INTENSITY = " << intensity << endl << endl << endl;
  
}


void
AmpToolsInterface::printEventDetails(string reactionName, Kinematics* kin) const {
  
  printKinematics(reactionName,kin);
  printAmplitudes(reactionName,kin);
  printIntensity(reactionName,kin);
  
}

float
AmpToolsInterface::random( float randMax ){
  
  return ( (float) rand() / RAND_MAX ) * randMax;
}
