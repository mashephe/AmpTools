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

    virtual ~AmpToolsInterface() { clear(); }

    static void registerAmplitude( const Amplitude& defaultAmplitude);

    static void registerDataReader( const DataReader& defaultDataReader);

    void resetConfigurationInfo(ConfigurationInfo* cfgInfo);

    ConfigurationInfo* configurationInfo()
                                 { return m_configurationInfo; }

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

    void printKinematics   (string reactionName, Kinematics* kin);
    void printAmplitudes   (string reactionName, Kinematics* kin);
    void printIntensity    (string reactionName, Kinematics* kin);
    void printEventDetails (string reactionName, Kinematics* kin);
    void printTestEvents(string reactionName, DataReader* dataReader, int nEvents = 10);

    double processEvents(string reactionName, DataReader* dataReader);

    double intensity(int iEvent);

    complex<double> decayAmplitude (int iEvent, string ampName);

    complex<double> productionAmplitude (string ampName);

    Kinematics* kinematics(int iEvent) {return m_ampVecs.getEvent(iEvent);}

  protected:

    void clear();

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

    AmpVecs m_ampVecs;

    AmplitudeManager* m_ampVecsAmpManager;

};



#endif
