#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include "IUAmpTools/LHContributionManager.h"

using namespace std;

LHContributionManager::LHContributionManager(){
  cout << endl << "## LHContribution Manager INITIALIZATION ##" << endl;
}

LHContributionManager::~LHContributionManager() {

}

void
LHContributionManager::registerLHContribution( const LHContribution& lhcontribution ){
  m_registeredLHContributions[lhcontribution.name()] = lhcontribution.clone();
}

void LHContributionManager::addLHContribution(const string& lhcontName,
                     const vector< string >& args){
  

  map< string, LHContribution* >::iterator defaultLHContribution =
     m_registeredLHContributions.find( lhcontName );


  cout << lhcontName << " " << m_registeredLHContributions.count( lhcontName ) << " " << __FILE__ << " " << __LINE__ << endl;
  if( defaultLHContribution == m_registeredLHContributions.end() ){
    
    cout << "ERROR: LHContribution with name " << lhcontName
	 << " has not been registered." << endl;
    assert( false );
  }   

  LHContribution* newLHContribution = defaultLHContribution->second->newLHContribution( args );
  m_mapNameToLHContributions[lhcontName] = newLHContribution;

}

void
LHContributionManager::setParValue( const string& name, const string& parName,
                               double val ){
   
  m_mapNameToLHContributions[name]->setParValue( parName, val );
}

void LHContributionManager::setupFromConfigurationInfo( const ConfigurationInfo* configInfo ){


  vector<LHContributionInfo*> lhcontInfoVector = configInfo->LHContributionList();
  for (unsigned int i = 0; i < lhcontInfoVector.size(); i++){
    string lhcontName = lhcontInfoVector[i]->fullName();
    vector< ParameterInfo* > pars = lhcontInfoVector[i]->parameters();
    vector< vector<string> > lhcontFactors = lhcontInfoVector[i]->factors();

    for (unsigned int j = 0; j < lhcontFactors.size(); j++){
      string factorName = lhcontFactors[j][0];
      vector<string> lhcontParameters = lhcontFactors[j];
      lhcontParameters.erase(lhcontParameters.begin());
      addLHContribution(lhcontName,lhcontParameters);
    }

    for( vector< ParameterInfo* >::iterator parItr = pars.begin(); parItr != pars.end(); ++parItr ){
      setParValue( lhcontName, (**parItr).parName(), (**parItr).value() );
    }

/*
    // finally initialize the LHContributions
    vector< const Amplitude* > ampVec = m_mapNameToAmps[ampName];
    for( vector< const Amplitude* >::iterator amp = ampVec.begin();
        amp != ampVec.end();
        ++amp ) {
      
      // init needs to be non-const or else the user has to deal with
      // mutable data -- in reality we really want some aspects of the
      // amplitude like the name, arguments, etc. to never change and
      // other aspects to be mutable, but this seems to put an extra
      // burden on the user
      
      const_cast< Amplitude* >(*amp)->init();
    }
  */
  }
}


