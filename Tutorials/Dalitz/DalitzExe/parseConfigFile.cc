#include <iostream>
#include <string>
#include "IUAmpTools/ConfigFileParser.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"

#include "IUAmpTools/report.h"
static const char* kModule = "parseConfigFile";

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  if (argc <= 1){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\tparseConfigFile <config file>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Parse the Config File *** " << endl << endl;

    // ************************
    // parse the command line parameters
    // ************************

  string configfile(argv[1]);

  report( INFO, kModule ) << "config file name = " << configfile << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser::setVerboseParsing(true);

    // method 1
  //ConfigFileParser parser(configfile);

    // method 2
  ConfigFileParser parser;
  ifstream infile(configfile.c_str());
  infile >> parser;
  infile.close();

  parser.getConfigurationInfo()->display();

  parser.getConfigurationInfo()->write("testWrite.cfg");

}
