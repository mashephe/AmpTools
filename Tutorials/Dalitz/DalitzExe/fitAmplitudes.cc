#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzAmp/constraint.h"

#include "IUAmpTools/report.h"
static const char* kModule = "fitAmplitudes";


using std::complex;
using namespace std;

int main( int argc, char* argv[] ){

  if (argc <= 1){
    report( INFO, kModule ) << "Usage:" << endl << endl;
    report( INFO, kModule ) << "\tfitAmplitudes <config file name>" << endl << endl;
    return 0;
  }

    // ************************
    // usage
    // ************************

  report( INFO, kModule ) << " *** Performing the Fit *** " << endl;


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  report( INFO, kModule ) << "Config file name:  " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();


    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterface::registerAmplitude(BreitWigner());
  AmpToolsInterface::registerLHContribution(constraint());
  AmpToolsInterface::registerDataReader(DalitzDataReader());

  AmpToolsInterface ATI(cfgInfo);

  double neg2LL = ATI.likelihood();
  report( INFO, kModule ) << "-2 ln(L) BEFORE MINIMIZATION:  " << neg2LL << endl;

  MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
  fitManager->setStrategy(1);

  fitManager->migradMinimization();

  if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
    report( WARNING, kModule ) << "Fit failed." << endl;
  }

  report( INFO, kModule ) << "-2 ln(L) AFTER MINIMIZATION:  " << ATI.likelihood() << endl;

  ATI.finalizeFit();

  return 0;

}


