#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzDataIO/DalitzDataWriter.h"

#include "IUAmpTools/report.h"
static const char* kModule = "generatePhaseSpace";

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************


  if (argc <= 2){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\tgeneratePhaseSpace <output file name> <number of events>" << endl << endl;
    return 0;
  }

  cout << endl << " *** Generate Phase Space *** " << endl << endl;


    // ************************
    // parse the command line parameters
    // ************************

  string filename(argv[1]);
  int nevents(atoi(argv[2]));

  report( INFO, kModule ) << "Output file name = " << filename << endl;
  report( INFO, kModule ) << "Number of events = " << nevents << endl << endl;


    // ************************
    // set up a DalitzDataWriter object
    // ************************

  report( DEBUG, kModule ) << "Creating a Data Writer..." << endl;

  DalitzDataWriter dataWriter(filename);

  report( DEBUG, kModule ) << "... Finished creating a Data Writer" << endl << endl;


    // ************************
    // use ROOT to generate phase space
    // ************************

  TGenPhaseSpace generator;
  TLorentzVector parent(0.0, 0.0, 0.0, 3.0);
  double daughterMasses[3] = {0.2, 0.2, 0.2};
  generator.SetDecay(parent, 3, daughterMasses);
  double maxWeight = generator.GetWtMax();


    // ************************
    // use the DalitzDataWriter object to write events to a file
    // ************************

  for (int i = 0; i < nevents; ++i){


      // generate the decay

    double weight = generator.Generate();
    if (weight < drand48() * maxWeight){
      i--;  continue;
    }


      // pack the decay products into a Kinematics object

    vector<TLorentzVector> fourvectors;
    fourvectors.push_back( TLorentzVector( *generator.GetDecay(0) ) );
    fourvectors.push_back( TLorentzVector( *generator.GetDecay(1) ) );
    fourvectors.push_back( TLorentzVector( *generator.GetDecay(2) ) );
    Kinematics kin(fourvectors);


      // write to a file

    dataWriter.writeEvent(kin);

    if (dataWriter.eventCounter() % 1000 == 0)
      report( INFO, kModule ) << "Event counter = " << dataWriter.eventCounter() << endl;

  }


}
