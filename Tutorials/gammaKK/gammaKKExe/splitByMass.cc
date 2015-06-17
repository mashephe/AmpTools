#include <iostream>
#include <iomanip>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "gammaKKDataIO/gammaKKDataReader.h"
#include "gammaKKDataIO/gammaKKDataWriter.h"
#include "gammaKKExe/constants.h"

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Split ROOT file into mass bins *** " << endl << endl;

  if (argc != 4){
    cout << "Usage:" << endl << endl;
    cout << "\tsplitByMass <config file name> <infile name> <output file name base>" << endl << endl;
    cout << "\t\t (The binning should be given by the config file)" << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  const string infilename(argv[2]);
  const string outfilenamebase(argv[3]);


  cout << "Config file name      = " << cfgname << endl << endl;
  cout << "Input file name       = " << infilename << endl;
  cout << "Output file name base = " << outfilenamebase << endl;

    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  // Set up the binning from the config file
  // (search for keyword 'binning')
  int nbins = -999;
  double min, max;
  vector<ConfigFileLine> lines = parser.getConfigFileLines();
  for(unsigned int i=0;i<lines.size();i++){
    if(lines[i].keyword() == "binning"){
      vector<string>args = lines[i].arguments();
      nbins = atoi(args[0].c_str());
      min = atof(args[1].c_str());
      max = atof(args[2].c_str());
    }
  }
  assert(nbins != -999);

  // Now that we have created the AmpToolsInterface,
  // we can use it over and over in different bins by
  // resetting the parameters.
  const int NBINS = nbins;
  const double MIN = min;
  const double MAX = max;
  const double BINNING = (MAX - MIN) / NBINS;
  cout << "- - - - - - - - - - - - Mass binning - - - - - - - - - - - - " << endl;
  cout << " NBINS    = " << NBINS << endl;
  cout << " MIN MASS = " << MIN   << endl;
  cout << " MAX MASS = " << MAX   << endl;
  cout << "- - - - - - - - - - - - Mass binning - - - - - - - - - - - - " << endl;

  ReactionInfo* reaction = cfgInfo->reactionList()[0];

    // ************************
    // create multiple DataWriters
    // ************************

  cout << "Creating Data Writers..." << endl;
  gammaKKDataWriter *dataWriter[NBINS];
  char filename[200];
  for(int i=0;i<NBINS;i++){
    sprintf(filename,"%s_%2.2d.root",outfilenamebase.c_str(),i);
    dataWriter[i] = new gammaKKDataWriter(filename);
  }
  cout << "... Finished creating Data Writers" << endl << endl;

    // ************************
    // read in ROOT file
    // ************************
  vector<string> args;
  args.push_back(infilename);
  gammaKKDataReader dataReader(args);
  Kinematics* kin;
  while( (kin = dataReader.getEvent()) ){

    if(dataReader.eventCounter() % 50000==0)
      cout << "processing event " << setw(12) << dataReader.eventCounter() << endl;

    // Calculate bin for this event
    TLorentzVector p2 = kin->particle(1);
    TLorentzVector p3 = kin->particle(2);
    double m23 = (p2 + p3).M();
    int bin = static_cast<int>( floor((m23 - MIN) / BINNING));
    if(!(0<=bin && bin<NBINS)){
      // If the bin is out of range, show what the mass was
      cout << "m23 = " << m23 << " bin = " << bin << endl;
      continue;
    }
    
    dataWriter[bin]->writeEvent(*kin);
    delete kin;
  }


  // For having multiple dataWriter objects,
  // it seems it is necessary to call the write of the ROOT file
  // and TTree explicitly, instead of having the destructor
  // take care of it.
  int totalEvents = 0;
  for(int i=0;i<NBINS;i++){
    cout << "KEPT " << setw(12) << dataWriter[i]->eventCounter() << " events in bin " << i << endl;
    dataWriter[i]->write();
    totalEvents += dataWriter[i]->eventCounter();
  }
  cout << "KEPT " << setw(12) << totalEvents << " TOTAL events" << endl;



}
