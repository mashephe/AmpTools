#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TFile.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzPlot/DalitzPlotGenerator.h"

using std::complex;
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Plotting Results from the Fit *** " << endl << endl;

  if (argc <= 3){
    cout << "Usage:" << endl << endl;
    cout << "\tplotResults <config file name> <fit results file> <output file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  string parname(argv[2]);
  string outname(argv[3]);

  cout << "Config file name    = " << cfgname << endl;
  cout << "Parameter file name = " << parname << endl << endl;
  cout << "Output file name    = " << outname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  ReactionInfo* reaction = cfgInfo->reactionList()[0];


    // ************************
    // set up an output Root file
    // ************************

  TFile* plotfile = new TFile( outname.c_str(), "recreate");
  TH1::AddDirectory(kFALSE);


    // ************************
    // set up a PlotGenerator and make plots
    // ************************

  AmpToolsInterface::registerDataReader( DalitzDataReader() );
  AmpToolsInterface::registerAmplitude( BreitWigner() );

  AmpToolsInterface ati( cfgInfo );
  
  DalitzPlotGenerator plotGenerator( ati );
  plotGenerator.enableReaction(reaction->reactionName());
  vector<string> amps = plotGenerator.uniqueAmplitudes();


    // loop over amplitude configurations

  for (unsigned int iamp = 0; iamp <= amps.size(); iamp++){

      // turn on all amplitudes

    for (unsigned int i = 0; i < amps.size(); i++){
      plotGenerator.enableAmp(i);
    }

      // turn off all amplitudes but one (or keep all on)

    if (iamp < amps.size()){
      for (unsigned int i = 0; i < amps.size(); i++){
        if (i != iamp) plotGenerator.disableAmp(i);
      }      
    }


    // loop over data, accMC, and genMC

  for (unsigned int iplot = 0; iplot < PlotGenerator::kNumTypes; iplot++){
    if (iamp < amps.size() && iplot == PlotGenerator::kData) continue;


    // loop over different variables

  for (unsigned int ivar  = 0; ivar  < DalitzPlotGenerator::kNumHists; ivar++){

             string histname =  "h";
    if (ivar == DalitzPlotGenerator::khm12)  histname += "m12";
    if (ivar == DalitzPlotGenerator::khm13)  histname += "m13";
    if (ivar == DalitzPlotGenerator::khm23)  histname += "m23";
    if (iplot == PlotGenerator::kData) histname += "dat";
    if (iplot == PlotGenerator::kAccMC) histname += "acc";
    if (iplot == PlotGenerator::kGenMC) histname += "gen";
    if (iamp < amps.size()){
      ostringstream sdig;  sdig << (iamp + 1);
      histname += sdig.str();
    }

    Histogram hist = plotGenerator.projection(ivar,
                         reaction->reactionName(), iplot);

            string xtitle = "Mass(P_{1}P_{2}) (GeV/c^{2})";
    if (ivar == 1) xtitle = "Mass(P_{1}P_{3}) (GeV/c^{2})";
    if (ivar == 2) xtitle = "Mass(P_{2}P_{3}) (GeV/c^{2})";

    TH1F thist = hist.toRoot();
    thist.SetName(histname.c_str());
    thist.SetStats(0);
    thist.SetTitleOffset(1.9,"Y");
    thist.SetTitleOffset(1.9,"X");
    thist.SetTitle("IUAmpTools Dalitz Tutorial");
    thist.SetXTitle(xtitle.c_str());
    thist.SetYTitle("Events / 50 MeV/c^{2}");
    plotfile->cd();
    thist.Write();

  }}}

  plotfile->Close();



    // ************************
    // print results to the screen
    // ************************


  cout << "TOTAL EVENTS = " << plotGenerator.intensity().first << " +- "
                            << plotGenerator.intensity().second << endl;
  vector<string> fullamps = plotGenerator.fullAmplitudes();
  for (unsigned int i = 0; i < fullamps.size(); i++){
    vector<string> useamp;  useamp.push_back(fullamps[i]);
    cout << "FIT FRACTION " << i+1 << " = "
         << plotGenerator.intensity(useamp).first /
            plotGenerator.intensity().first <<  " +- "
         << plotGenerator.intensity(useamp).second /
            plotGenerator.intensity().first <<  endl;
  }

  return 0;

}

