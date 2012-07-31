#if !defined(AMPTOOLSINTERFACE)
#define AMPTOOLSINTERFACE

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <fstream>

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"

using namespace std;


class AmpToolsInterface{

  public:

    AmpToolsInterface() { }

    AmpToolsInterface(ConfigurationInfo* cfgInfo);

    virtual ~AmpToolsInterface();

    static void registerAmplitude( const Amplitude& defaultAmplitude);

    static void registerDataReader( const DataReader& defaultDataReader);


    MinuitMinimizationManager* minuitMinimizationManager() 
                                 { return m_minuitMinimizationManager; }

    ParameterManager*          parameterManager()          
                                 { return m_parameterManager;}


    AmplitudeManager*     amplitudeManager     (const string& reactionName);

    DataReader*           dataReader           (const string& reactionName);
    DataReader*           accMCReader          (const string& reactionName);
    DataReader*           genMCReader          (const string& reactionName);

    NormIntInterface*     normIntInterface     (const string& reactionName);

    LikelihoodCalculator* likelihoodCalculator (const string& reactionName);

    double likelihood();
    double likelihood(const string& reactionName);

    virtual void finalizeFit();

    void printTestAmplitudes (int nEvents = 10, string dataType = "genMC");



/*
  calculateIntensities?
*/


  protected:

    ConfigurationInfo*          m_configurationInfo;
    MinuitMinimizationManager*  m_minuitMinimizationManager;
    ParameterManager*           m_parameterManager;

    vector<AmplitudeManager*>  m_amplitudeManagers;

    map<string,DataReader*> m_dataReaderMap;
    map<string,DataReader*> m_genMCReaderMap;
    map<string,DataReader*> m_accMCReaderMap;

    map<string,NormIntInterface*>     m_normIntMap;
    map<string,LikelihoodCalculator*> m_likCalcMap;

    static vector<Amplitude*>  m_userAmplitudes;
    static vector<DataReader*> m_userDataReaders;


};



#endif
