#include <iostream>
#include <string>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "DalitzPlot/DalitzPlotGenerator.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

typedef DalitzPlotGenerator PlotGen;

void atiSetup(){
  
  // the PlotGenerator will create an AmpToolsInterface in order
  // to create plots - this setup must happen before the
  // AmpToolsInterface is created

  AmpToolsInterface::registerAmplitude( BreitWigner() );
  AmpToolsInterface::registerDataReader( DalitzDataReader() );
}


//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Viewing Results Using AmpPlotter *** " << endl << endl;

  if (argc <= 1){
    
    cout << "Usage:" << endl << endl;
    cout << "\tampPlotter <fit results name>" << endl << endl;
    return 0;
  }

    // ************************
    // parse the command line parameters
    // ************************

  string resultsName(argv[1]);
  FitResults results( resultsName );
  if( !results.valid() ){
    
    cout << "Invalid fit results in file:  " << resultsName << endl;
    exit( 1 );
  }

    // ************************
    // set up the plot generator
    // ************************

  atiSetup();  
  PlotGen plotGen( results );

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

