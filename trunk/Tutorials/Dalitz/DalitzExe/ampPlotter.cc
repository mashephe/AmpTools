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

//#define NOMPI

using namespace std;

void printUsage( char* exe ){

  cout << endl << " Usage for: " << exe << endl << endl;
  cout << "\t -c \t Config file" << endl << endl;
}

int main( int argc, char* argv[] ){    
    
  string cfgFile( "" );
 
 // parse command line
  
  if( argc == 1 ){

    printUsage( argv[0] );
    exit( 1 );
  }

  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
    else  cfgFile = argv[++i]; }
    if( (arg == "-h") ){
      printUsage( argv[0] );
      exit(1);
    }
  }
  
  ConfigFileParser parser( cfgFile );
  const ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();  

  string parFile( cfgInfo->fitName() );
  parFile += ".fit";
  
  cout << ">> Using fit output: " << parFile << endl;

  PlotGen plotGen( cfgInfo, parFile );
  
  cout << ">> Plot generator ready, starting GUI.." << endl;

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

