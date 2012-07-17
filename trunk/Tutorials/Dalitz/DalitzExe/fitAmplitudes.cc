#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"


using std::complex;
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Performing the Fit *** " << endl << endl;

  if (argc <= 1){
    cout << "Usage:" << endl << endl;
    cout << "\tfitAmplitudes <config file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  cout << "Config file name = " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  ReactionInfo* reaction = cfgInfo->reactionList()[0];


    // ************************
    // create an AmplitudeManager
    // ************************

  cout << endl << endl;
  cout << "Creating AmplitudeManager for reaction " << reaction->reactionName() << endl;
  vector<AmplitudeManager*>  ampManagers;
  ampManagers.push_back(new AmplitudeManager(reaction->particleList(),reaction->reactionName()));
  ampManagers[0]->registerAmplitudeFactor( BreitWigner() );
  ampManagers[0]->setupFromConfigurationInfo( cfgInfo );
  cout << "... Finished creating AmplitudeManager" << endl;


    // ************************
    // create a MinuitMinimizationManager
    // ************************

  MinuitMinimizationManager* fitManager = new MinuitMinimizationManager(100);
  fitManager->setPrecision( 1E-13 );


    // ************************
    // create a ParameterManager
    // ************************

  ParameterManager parManager( *fitManager, ampManagers );
  parManager.setupFromConfigurationInfo( cfgInfo );


    // ************************
    // create DataReaders
    // ************************

  DalitzDataReader dataDataReader(reaction->dataFiles()[0]);
  DalitzDataReader genDataReader(reaction->genMCFiles()[0]);
  DalitzDataReader accDataReader(reaction->accMCFiles()[0]);


    // ************************
    // create a NormIntInterface
    // ************************

  NormIntInterface normInt(&genDataReader, &accDataReader, *ampManagers[0]);


    // ************************
    // create a LikelihoodCalculator
    // ************************

  LikelihoodCalculator likeCalc(*ampManagers[0], normInt, dataDataReader, parManager); 
  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << likeCalc() << endl;


    // ************************
    // do the minimization
    // ************************

  cout << "STARTING MINIMIZATION..." << endl;

  fitManager->setStrategy( 1 );
  fitManager->migradMinimization();
  
  if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
    cout << "ERROR: fit failed..." << endl;
    return 1;
  }

  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << likeCalc() << endl;


    // ************************
    // save fit parameters
    // ************************

  cout << "WRITING PARAMETERS..." << endl;
  
  string outputFile(cfgInfo->fitName());  outputFile += ".fit";
  ofstream outFile(outputFile.c_str());
  parManager.writeParameters( outFile );


    // ************************
    // save normalization integrals
    // ************************
  
  cout << "WRITING FINAL NORMALIZATION INTEGRALS.." << endl;
  
  normInt.forceCacheUpdate();
  normInt.exportNormIntCache( reaction->normIntFile() );


  return 0;
}


