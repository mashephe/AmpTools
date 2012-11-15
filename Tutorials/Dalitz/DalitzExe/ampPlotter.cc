#include <iostream>
#include <string>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"

#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "DalitzPlot/DalitzPlotGenerator.h"

typedef DalitzPlotGenerator PlotGen;

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter *** " << endl << endl;

  if (argc <= 1){
    cout << "Usage:" << endl << endl;
    cout << "\tampPlotter <config file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  cout << "Config file name    = " << cfgname << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  string parFile( cfgInfo->fitOutputFileName() );


    // ************************
    // set up the plot generator
    // ************************

  
  cout << "Setting up the plot generator with fit output: " << parFile << endl;

  PlotGen plotGen( cfgInfo, parFile );
  

    // ************************
    // start the GUI
    // ************************

  cout << ">> Plot generator ready, starting GUI..." << endl;

  int dummy_argc = 0;
  char* dummy_argv[] = {};  
  TApplication app( "app", &dummy_argc, dummy_argv );
  
  gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetFillStyle(1001);
  gStyle->SetPalette(1);
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameFillStyle(1001);
  
  PlotFactory factory( plotGen );	
  PlotterMainWindow mainFrame( gClient->GetRoot(), factory );
	
  app.Run();
    
  return 0;

}

