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
#include "gammaKKDataIO/gammaKKDataReader.h"
#include "gammaKKPlot/gammaKKPlotGenerator.h"

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

      //xxxxx once an amplitude is disabled it is always disabled????
      //        can't turn an amplitude off then back on?????

  gammaKKPlotGenerator plotGenerator(cfgInfo,parname);
  plotGenerator.enableReaction(reaction->reactionName());
  vector<string> amps = plotGenerator.uniqueAmplitudes();

  for(int i=0;i<amps.size();i++){
    cout << "amps[" << i << "] = " << amps[i] << endl;
  }


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

    for (unsigned int iplot = 0; iplot < PlotGenerator::kNumPlotTypes; iplot++){
    if (iamp < amps.size() && iplot == 0) continue;


    // loop over different variables

    for (unsigned int ivar  = 0; ivar  < gammaKKPlotGenerator::kNumHist; ivar++){

             string histname =  "h";
    if (ivar == 0)  histname += "m12";
    if (ivar == 1)  histname += "m13";
    if (ivar == 2)  histname += "m23";
    if (ivar == 3)  histname += "PhotonCosTheta";
    if (ivar == 4)  histname += "PhotonPhi";
    if (ivar == 5)  histname += "KaonCosTheta";
    if (ivar == 6)  histname += "KaonPhi";

    if (iplot == 0) histname += "dat";
    if (iplot == 1) histname += "acc";
    if (iplot == 2) histname += "gen";
    if (iamp < amps.size()){
      ostringstream sdig;  sdig << (iamp + 1);
      histname += sdig.str();
    }

    // data type
    PlotGenerator::PlotType kplot = PlotGenerator::kData;
    if (iplot == 1) kplot = PlotGenerator::kAccMC;
    if (iplot == 2) kplot = PlotGenerator::kGenMC;

    // histogram type    
    gammaKKPlotGenerator::HistIndex kvar = gammaKKPlotGenerator::khm12; 
    if (ivar == 1) kvar = gammaKKPlotGenerator::khm13;
    if (ivar == 2) kvar = gammaKKPlotGenerator::khm23;
    if (ivar == 3) kvar = gammaKKPlotGenerator::khPhotonCosTheta;
    if (ivar == 4) kvar = gammaKKPlotGenerator::khPhotonPhi;
    if (ivar == 5) kvar = gammaKKPlotGenerator::khKaonCosTheta;
    if (ivar == 6) kvar = gammaKKPlotGenerator::khKaonPhi;
    
    Histogram hist = plotGenerator.projection(kvar,
					      reaction->reactionName(), kplot);
    
    string xtitle = "Mass(P_{1}P_{2}) (GeV/c^{2})";
    if (ivar == 1) xtitle = "Mass(P_{1}P_{3}) (GeV/c^{2})";
    if (ivar == 2) xtitle = "Mass(P_{2}P_{3}) (GeV/c^{2})";
    if (ivar == 3) xtitle = "cos#theta_{#gamma}";
    if (ivar == 4) xtitle = "#phi_{#gamma}";
    if (ivar == 5) xtitle = "cos#theta_{K_{1}}";
    if (ivar == 6) xtitle = "#phi_{K_{1}}";

    TH1F thist = hist.toRoot();
    thist.SetName(histname.c_str());
    thist.SetStats(0);
    thist.SetTitleOffset(1.9,"Y");
    thist.SetTitleOffset(1.9,"X");
    thist.SetTitle("IUAmpTools gammaKK Tutorial");
    thist.SetXTitle(xtitle.c_str());
    thist.SetYTitle("counts");
    plotfile->cd();
    thist.Write();

  }}}

  plotfile->Close();



    // ************************
    // print results to the screen
    // ************************

       // xxxxx error on the total intensity < sqrt(intensity) 
       //         when using 100% acceptance????

  cout << "TOTAL EVENTS = " << plotGenerator.intensity().first << " +- "
                            << plotGenerator.intensity().second << endl;
  vector<string> fullamps = plotGenerator.fullAmplitudes();

  for(int i=0;i<fullamps.size();i++){
    cout << "fullamps[" << i << "] = " << fullamps[i] << endl;
  }

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

