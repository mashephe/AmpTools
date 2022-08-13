#include <iostream>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

using namespace std;

#include "IUAmpTools/report.h"
static const char* kModule = "printAmplitudes";

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  if (argc <= 1){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\tprintAmplitudes <config file name>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Printing Amplitudes *** " << endl << endl;

    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  report( INFO, kModule ) << "Config file name = " << cfgname << endl << endl;


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
  AmpToolsInterface::registerDataReader(DalitzDataReader());

  AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kPlotGeneration);

  DataReader* dataReader = ATI.genMCReader(cfgInfo->reactionList()[0]->reactionName());
  for (int i = 0; i < 10; i++){
    Kinematics* kin = dataReader->getEvent();
    ATI.printEventDetails(cfgInfo->reactionList()[0]->reactionName(),kin);
    delete kin;
  }

}
