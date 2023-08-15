#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include "IUAmpTools/Neg2LnLikContribManager.h"
#include "IUAmpTools/report.h"

const char* Neg2LnLikContribManager::kModule = "Neg2LnLikContribManager";

using namespace std;

Neg2LnLikContribManager::Neg2LnLikContribManager(){
  report( DEBUG, kModule ) << "## Neg2LnLikContrib Manager INITIALIZATION ##" << endl;
}

Neg2LnLikContribManager::~Neg2LnLikContribManager() {

}

void
Neg2LnLikContribManager::registerNeg2LnLikContrib( const Neg2LnLikContrib& lhcontribution ){
  m_registeredNeg2LnLikContribs[lhcontribution.name()] = lhcontribution.clone();
}

void Neg2LnLikContribManager::addNeg2LnLikContrib(const string& lhcontName,
                     const vector< string >& args){
  

  map< string, Neg2LnLikContrib* >::iterator defaultNeg2LnLikContrib =
     m_registeredNeg2LnLikContribs.find( lhcontName );


  if( defaultNeg2LnLikContrib == m_registeredNeg2LnLikContribs.end() ){
    
    report( ERROR, kModule ) << "Neg2LnLikContrib with name " << lhcontName
                             << " has not been registered." << endl;
    assert( false );
  }   

  Neg2LnLikContrib* newNeg2LnLikContrib = defaultNeg2LnLikContrib->second->newNeg2LnLikContrib( args );
  m_mapNameToNeg2LnLikContribs[lhcontName] = newNeg2LnLikContrib;

}

void
Neg2LnLikContribManager::setParPtr( const string& name, const string& parName,
                            const double* ampParPtr ){

  m_mapNameToNeg2LnLikContribs[name]->setParPtr( parName, ampParPtr );
}

void
Neg2LnLikContribManager::setParValue( const string& name, const string& parName,
                               double val ){
   
  m_mapNameToNeg2LnLikContribs[name]->setParValue( parName, val );
}

void Neg2LnLikContribManager::setupFromConfigurationInfo( const ConfigurationInfo* configInfo ){


  vector<Neg2LnLikContribInfo*> lhcontInfoVector = configInfo->neg2LnLikContribList();
  for (unsigned int i = 0; i < lhcontInfoVector.size(); i++){
    string lhcontName = lhcontInfoVector[i]->fullName();
    vector<string> lhcontParameters = lhcontInfoVector[i]->args();
    lhcontParameters.erase(lhcontParameters.begin());
    addNeg2LnLikContrib(lhcontName,lhcontParameters);

    vector< ParameterInfo* > pars = lhcontInfoVector[i]->parameters();
    for( vector< ParameterInfo* >::iterator parItr = pars.begin(); parItr != pars.end(); ++parItr ){
      setParValue( lhcontName, (**parItr).parName(), (**parItr).value() );
    }

    // finally initialize the Neg2LnLikContribs
    const_cast<Neg2LnLikContrib*>(m_mapNameToNeg2LnLikContribs[lhcontName])->init();
  }
}


