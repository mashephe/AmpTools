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


  if( defaultLHContribution == m_registeredLHContributions.end() ){
    
    cout << "ERROR: LHContribution with name " << lhcontName
	 << " has not been registered." << endl;
    assert( false );
  }   

  LHContribution* newLHContribution = defaultLHContribution->second->newLHContribution( args );
  m_mapNameToLHContributions[lhcontName] = newLHContribution;

}

void
LHContributionManager::setParPtr( const string& name, const string& parName,
                            const double* ampParPtr ){

  m_mapNameToLHContributions[name]->setParPtr( parName, ampParPtr );
}

void
LHContributionManager::setParValue( const string& name, const string& parName,
                               double val ){
   
  m_mapNameToLHContributions[name]->setParValue( parName, val );
}
/*
void
LHContributionManager::updatePar( const string& parName ) const {
  
  for( map< string, LHContribution* >::const_iterator mapItr = m_mapNameToLHContributions.begin();
      mapItr != m_mapNameToLHContributions.end();
      ++mapItr ){
    

      
      // if we find an amplitude with the parameter update the iteration
      // counter; this may result in multiple increments over one fuction
      // call but that is OK -- iteration numbers just need to be unique
      if( (*mapItr).second->updatePar( parName ) ){
        
       // ++m_ampIteration[*ampItr];
      
    }
  }
}
*/

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

    // finally initialize the LHContributions
    const_cast<LHContribution*>(m_mapNameToLHContributions[lhcontName])->init();
  }
}


