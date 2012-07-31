#include <mpi.h>
#include <pthread.h>

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
    m_amplitudeManagers.push_back(ampMan);

  }

    // ************************
    // create a ParameterManager
    // ************************

  ParameterManagerMPI* parameterManagerMPI = NULL;
  if (m_rank != 0){
    parameterManagerMPI = new ParameterManagerMPI( m_amplitudeManagers );
  }
  else{
    parameterManagerMPI = new ParameterManagerMPI( *m_minuitMinimizationManager,
                                                   m_amplitudeManagers );
  }
  parameterManagerMPI->setupFromConfigurationInfo( m_configurationInfo );
  m_parameterManager = parameterManagerMPI;




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
      normInt = new NormIntInterfaceMPI(genMCRdr, accMCRdr, *ampMan);
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

    LikelihoodCalculatorMPI* likCalc = NULL;
    if (ampMan && normInt && dataRdr && m_parameterManager){
      likCalc = new LikelihoodCalculatorMPI(*ampMan, *normInt, *dataRdr, *parameterManagerMPI); 
      m_likCalcMap[reactionName] = likCalc;
    }
    else{
      cout << "AmpToolsInterface WARNING:  not creating a LikelihoodCalculator for reaction " 
       << reactionName << endl;
    }


    if (likCalc && m_rank != 0) LikelihoodManagerMPI::deliverLikelihood();

  }

}





void
AmpToolsInterfaceMPI::finalizeFit(){

  if (m_rank != 0){
    for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
      ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
      string reactionName(reaction->reactionName());
      NormIntInterface* normInt = normIntInterface(reactionName);
      normInt->forceCacheUpdate();
    }
    return;
  }

  if (m_rank == 0){
    for (unsigned int irct = 0; irct < m_configurationInfo->reactionList().size(); irct++){
      ReactionInfo* reaction = m_configurationInfo->reactionList()[irct];
      string reactionName(reaction->reactionName());
      if (m_likCalcMap[reactionName]) delete m_likCalcMap[reactionName];
    }
    m_likCalcMap.clear();
  }

  AmpToolsInterface::finalizeFit();

}

