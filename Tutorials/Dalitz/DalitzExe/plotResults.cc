#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TFile.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/Histogram.h"
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

  DalitzPlotGenerator plotGenerator(cfgInfo,parname);
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

  for (unsigned int iplot = 0; iplot < 3; iplot++){
    if (iamp < amps.size() && iplot == 0) continue;


    // loop over different variables

  for (unsigned int ivar  = 0; ivar  < 3; ivar++){

             string histname =  "h";
    if (ivar == 0)  histname += "m12";
    if (ivar == 1)  histname += "m13";
    if (ivar == 2)  histname += "m23";
    if (iplot == 0) histname += "dat";
    if (iplot == 1) histname += "acc";
    if (iplot == 2) histname += "gen";
    if (iamp < amps.size()){
      ostringstream sdig;  sdig << (iamp + 1);
      histname += sdig.str();
    }

    PlotGenerator::PlotType kplot = PlotGenerator::kData;
            if (iplot == 1) kplot = PlotGenerator::kAccMC;
            if (iplot == 2) kplot = PlotGenerator::kGenMC;

    DalitzPlotGenerator::HistIndex kvar = DalitzPlotGenerator::khm12; 
                    if (ivar == 1) kvar = DalitzPlotGenerator::khm13;
                    if (ivar == 2) kvar = DalitzPlotGenerator::khm23;

    Histogram hist = plotGenerator.projection(kvar,
                         reaction->reactionName(), kplot);

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
    thist.SetYTitle("Events / 5 MeV/c^{2}");
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

