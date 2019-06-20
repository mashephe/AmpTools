#include <mpi.h>
#include <pthread.h>

#include "IUAmpTools/AmpToolsInterface.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"

#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpToolsMPI/LikelihoodCalculatorMPI.h"
#include "IUAmpToolsMPI/NormIntInterfaceMPI.h"
#include "IUAmpToolsMPI/AmpToolsInterfaceMPI.h"

AmpToolsInterfaceMPI::AmpToolsInterfaceMPI(ConfigurationInfo* configurationInfo){


  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &m_numProc );


  m_configurationInfo = configurationInfo;


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
    
    if( m_functionality == kFull ){
      ampMan->setOptimizeParIteration( true );
      ampMan->setFlushFourVecsIfPossible( true );
    }
    
    m_intensityManagers.push_back(ampMan);

  }

    // ************************
    // create a ParameterManager
    // ************************

  ParameterManagerMPI* parameterManagerMPI = NULL;
  if (m_rank != 0){
    parameterManagerMPI = new ParameterManagerMPI( m_intensityManagers );
  }
  else{
    parameterManagerMPI = new ParameterManagerMPI( m_minuitMinimizationManager,
                                                   m_intensityManagers );
  }
  parameterManagerMPI->setupFromConfigurationInfo( m_configurationInfo );
  m_parameterManager = parameterManagerMPI;




    // ************************
    // loop over reactions
    // ************************

  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());
    IntensityManager* intenMan = intensityManager(reactionName);


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

      // ************************
      // create a NormIntInterface
      // ************************

    NormIntInterface* normInt = NULL;
    if (genMCRdr && accMCRdr && intenMan && !(reaction->normIntFileInput())){
      normInt = new NormIntInterfaceMPI(genMCRdr, accMCRdr, *intenMan);
      m_normIntMap[reactionName] = normInt;
      if (reaction->normIntFile() == "")
	cout << "AmpToolsInterface WARNING:  no name given to NormInt file for reaction " 
	 << reactionName << endl;
    }
    else if (reaction->normIntFileInput()){
      normInt = new NormIntInterfaceMPI(reaction->normIntFile());
      m_normIntMap[reactionName] = normInt;
    }
    else{
      cout << "AmpToolsInterface WARNING:  not creating a NormIntInterface for reaction " 
       << reactionName << endl;
    }

      // ************************
      // create a LikelihoodCalculator
      // ************************

    LikelihoodCalculatorMPI* likCalc = NULL;
    if (intenMan && normInt && dataRdr && m_parameterManager){
      likCalc = new LikelihoodCalculatorMPI(*intenMan, *normInt, dataRdr, bkgndRdr, *parameterManagerMPI);
      m_likCalcMap[reactionName] = likCalc;
    }
    else{
      cout << "AmpToolsInterface ERROR:  not creating a LikelihoodCalculator for reaction " 
       << reactionName << endl;
      exit(1);
    }
  }

    // ************************
    // create FitResults
    // ************************

  m_fitResults = new FitResults( m_configurationInfo,
                                 m_intensityManagers,
                                 m_likCalcMap,
                                 m_normIntMap,
                                 m_minuitMinimizationManager,
                                 m_parameterManager );


  if (m_rank != 0) LikelihoodManagerMPI::deliverLikelihood();

}





void
AmpToolsInterfaceMPI::finalizeFit(){

  if (m_rank != 0){

    // this forceCacheUpdate on the workers executes simultaneously with the
    // same call in FitResults::writeResults on the master

    for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
      ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
      string reactionName(reaction->reactionName());
      NormIntInterface* normInt = normIntInterface(reactionName);
      if (normInt->hasAccessToMC()) normInt->forceCacheUpdate();
    }
  }

  if (m_rank == 0){

      // ************************
      // save fit parameters
      // ************************

    string outputFile(m_configurationInfo->fitOutputFileName());
    ofstream outFile(outputFile.c_str());
    m_fitResults->saveResults();
 
    // after saving the results destroy the LikelihoodCalculatorMPI
    // class on the master -- this will break the worker nodes out
    // of their deliverLikelihood loop and let them enter this routine

    for (unsigned int irct = 0;
         irct < m_configurationInfo->reactionList().size(); irct++){
      ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
      string reactionName(reaction->reactionName());
      if (m_likCalcMap.find(reactionName) != m_likCalcMap.end())
        delete m_likCalcMap[reactionName];
      m_likCalcMap.erase(reactionName);
    }

    m_fitResults->writeResults( outputFile );
 }
 
   // ************************
   // save normalization integrals
   // ************************

  // this runs on both the workers and the masters simultaneously

  for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){

    ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
    string reactionName(reaction->reactionName());
    NormIntInterface* normInt = normIntInterface(reactionName);
    // the call to FitResults::writeResults will force a cache update
    // there is no need to do it twice
    //      if (normInt->hasAccessToMC()) normInt->forceCacheUpdate();
    if( m_rank == 0 ) normInt->exportNormIntCache( reaction->normIntFile() );
  }
}

