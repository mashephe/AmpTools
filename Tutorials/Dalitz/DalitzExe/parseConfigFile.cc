#include <iostream>
#include <string>
#include "IUAmpTools/ConfigFileParser.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"


using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Parse the Config File *** " << endl << endl;
  
  if (argc <= 1){
    cout << "Usage:" << endl << endl;
    cout << "\tparseConfigFile <config file>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string configfile(argv[1]);

  cout << "config file name = " << configfile << endl << endl;


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
