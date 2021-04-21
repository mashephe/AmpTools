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


void LHContributionManager::setupFromConfigurationInfo( const ConfigurationInfo* configInfo ){

  vector<LHContributionInfo*> lhcontInfoVector = configInfo->LHContributionList();
  for (unsigned int i = 0; i < lhcontInfoVector.size(); i++){
    string lhcontName = lhcontInfoVector[i]->fullName();
    vector< ParameterInfo* > pars = lhcontInfoVector[i]->parameters();
    for( vector< ParameterInfo* >::iterator parItr = pars.begin(); parItr != pars.end(); ++parItr ){
      //setParValue( lhcontName, (**parItr).parName(), (**parItr).value() );
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


