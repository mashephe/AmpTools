#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"

#include "gammaKKPlot/gammaKKPlotGenerator.h"
#include "gammaKKDataIO/gammaKKDataReader.h"
#include "gammaKKAmp/gammaKKHelicityAmp.h"
#include "gammaKKAmp/TwoPiAngles.h"
#include "gammaKKAmp/MultipoleAmps.h"
#include "gammaKKExe/constants.h"

using std::complex;
using namespace std;
using namespace CLHEP;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Performing the Fit *** " << endl << endl;

  if (argc <= 1){
    cout << "Usage:" << endl << endl;
    cout << "\tfitAmplitudes <config file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  cout << "Config file name = " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();
  string fitname = cfgInfo->fitName();

  int nbins;
  double min, max;

  vector< vector<string> > args_binning = cfgInfo->userKeywordArguments("binning");
  // The keyword "binning" takes in 3 parameters,
  // the number of bins,
  // the min mass of KK,
  // the max mass of KK.
  assert(args_binning.size()==1);
  nbins = atoi(args_binning[0][0].c_str());
  min = atof(args_binning[0][1].c_str());
  max = atof(args_binning[0][2].c_str());

  // Now that we have created the AmpToolsInterface,
  // we can use it over and over in different bins by
  // resetting the parameters.
  const int NBINS = nbins;
  const double MIN = min;
  const double MAX = max;
  const double BINNING = (MAX - MIN) / NBINS;

    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterface::registerAmplitude(gammaKKHelicityAmp());
  AmpToolsInterface::registerAmplitude(TwoPiAngles());
  AmpToolsInterface::registerAmplitude(MultipoleAmps());

  AmpToolsInterface::registerDataReader(gammaKKDataReader());


  vector<string> args;
  string value;
  ostringstream ss;

  // A map that ties together the amplitude name and
  // the previous fit result (a complex double)
  //
  // For each fit done, this will be filled with the fitted
  // prodcution amplitude
  map<string, complex<double> > previous_fit;
  
  for(int i=0;i<NBINS;i++){

    cout << "******** Starting bin " << i << " **********" << endl;

    // Get the reaction name
    // (Do not change reaction name)
    assert(cfgInfo->reactionList().size()==1);
    string reactionName = (cfgInfo->reactionList())[0]->reactionName();

    // Set the fit name, so that the .fit file is set
    ss.clear();
    ss.str("");
    ss << fitname << "_" << setw(2) << setfill('0') << i;
    string fitname_new = ss.str();
    cfgInfo->setFitName(fitname_new);

    // Set the norm int file
    string normIntFileName = reactionName;
    ss.clear();
    ss.str("");
    ss << normIntFileName << "_" << setw(2) << setfill('0') << i << ".ni";
    normIntFileName = ss.str();
    cout << "normIntFileName = " << normIntFileName << endl;
    (cfgInfo->reactionList())[0]->setNormIntFile(normIntFileName);

    //////////////////////////////////////////////////////////////////    
    //    1. We will need to set the input data file                //
    //       1st argument is file name for data file                //
    //////////////////////////////////////////////////////////////////    
    args.clear();
    // ReactionInfo::data() returns pair<string, vector<string> >,
    // where the 1st string is the DataReader,
    // and the second <vector>string is the arguments with
    // which it was created.
    // The file name is the 1st element of the vector<string>.
    value = cfgInfo->reaction(reactionName)->data().second[0];
    cout << "data file = " << value << endl;
    args.push_back(value);
    // 2nd argument is minimum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * i;
    value = ss.str();
    args.push_back(value);
    // 3rd argument is maximum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * (i+1.);
    value = ss.str();
    args.push_back(value);
    cfgInfo->reaction(reactionName)->setData("gammaKKDataReader",args);

    //////////////////////////////////////////////////////////////////    
    //    2. Next, set up the MC file                               //
    //////////////////////////////////////////////////////////////////    
    args.clear();
    // ReactionInfo::accMC() returns pair<string, vector<string> >,
    // where the 1st string is the DataReader,
    // and the second <vector>string is the arguments with
    // which it was created.
    // The file name is the 1st element of the vector<string>.
    value = cfgInfo->reaction(reactionName)->accMC().second[0];
    args.push_back(value);
    // 2nd argument is minimum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * i;
    value = ss.str();
    args.push_back(value);
    // 3rd argument is maximum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * (i+1.);
    value = ss.str();
    args.push_back(value);
    cfgInfo->reaction(reactionName)->setAccMC("gammaKKDataReader",args);

    //////////////////////////////////////////////////////////////////    
    //    3. Set up the raw MC file                                 //
    //////////////////////////////////////////////////////////////////    
    args.clear();
    // ReactionInfo::genMC() returns pair<string, vector<string> >,
    // where the 1st string is the DataReader,
    // and the second <vector>string is the arguments with
    // which it was created.
    // The file name is the 1st element of the vector<string>.
    value = cfgInfo->reaction(reactionName)->genMC().second[0];
    args.push_back(value);
    // 2nd argument is minimum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * i;
    value = ss.str();
    args.push_back(value);
    // 3rd argument is maximum for mass range (as string)
    ss.clear();
    ss.str("");
    ss << MIN + BINNING * (i+1.);
    value = ss.str();
    args.push_back(value);
    cfgInfo->reaction(reactionName)->setGenMC("gammaKKDataReader",args);

    //////////////////////////////////////////////////////////////////    
    //    4. Finally, initialize the fit parameters according to    //
    //       the previous fit result.                               //
    //////////////////////////////////////////////////////////////////
    // For the 1st bin, the values are set in the configuration file.
    // For the other bins, use the previous fit result.
    if(i!=0){
      vector<AmplitudeInfo*> ampInfos = cfgInfo->amplitudeList();
      for (unsigned int iampInfo = 0; iampInfo < ampInfos.size(); iampInfo++){
        ampInfos[iampInfo]->setValue(previous_fit[ampInfos[iampInfo]->fullName()]); 
        ampInfos[iampInfo]->setReal(true); 
      }
    }

    // Show what config info is
    cout << "--------------------------------------------------------------------------" << endl;
    cfgInfo->display();
    cout << "--------------------------------------------------------------------------" << endl;

    //////////////////////////////////////////////////////////////////    
    //    5. Create the ATI manager using the cfgInfo               //
    //////////////////////////////////////////////////////////////////    
    AmpToolsInterface ATI(cfgInfo);

    // Set up MinuitMinimizationManager
    MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
    fitManager->setPrecision(1E-13);
    fitManager->setStrategy(1);
    
    cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ATI.likelihood() << endl;

    // Minimize!!
    fitManager->migradMinimization();

    bool fitSuccess = true;

    if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
      cout << "ERROR: fit failed..." << endl;
      fitSuccess = false;
    }

    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ATI.likelihood() << endl;
    
    ATI.finalizeFit();

    if(fitSuccess == true){
      // If the fit succeeded, save the current fit value
      // so we can use it for the next fit

      ss.clear();
      ss.str("");
      ss << fitname_new << ".fit";
      string fitfilename = ss.str();
      gammaKKPlotGenerator plotGenerator(cfgInfo,fitfilename);
      plotGenerator.enableReaction(reactionName);
      vector<string> amps = plotGenerator.fullAmplitudes();

      for(int i=0;i<amps.size();i++){
	previous_fit[amps[i]] = ATI.parameterManager()->findParameter(amps[i])->value();
	cout << "amplitude number " << i
	     << ", ampname =  " << amps[i]
	     << ", value = " << previous_fit[amps[i]] << endl;
      }
    }
    
  } // end of loop over different mass bins
  
  return 0;

}


