#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzDataIO/DalitzDataWriter.h"

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Generate Phase Space *** " << endl << endl;

  if (argc <= 2){
    cout << "Usage:" << endl << endl;
    cout << "\tgeneratePhaseSpace <output file name> <number of events>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string filename(argv[1]);
  int nevents(atoi(argv[2]));

  cout << "Output file name = " << filename << endl;
  cout << "Number of events = " << nevents << endl << endl;


    // ************************
    // set up a DalitzDataWriter object
    // ************************

  cout << "Creating a Data Writer..." << endl;

  DalitzDataWriter dataWriter(filename);

  cout << "... Finished creating a Data Writer" << endl << endl;


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

    vector<HepLorentzVector> fourvectors;
    TLorentzVector* p1 = generator.GetDecay(0);
    TLorentzVector* p2 = generator.GetDecay(1);
    TLorentzVector* p3 = generator.GetDecay(2);
    fourvectors.push_back(HepLorentzVector(p1->Px(),p1->Py(),p1->Pz(),p1->E()));
    fourvectors.push_back(HepLorentzVector(p2->Px(),p2->Py(),p2->Pz(),p2->E()));
    fourvectors.push_back(HepLorentzVector(p3->Px(),p3->Py(),p3->Pz(),p3->E()));
    Kinematics kin(fourvectors);


      // write to a file

    dataWriter.writeEvent(kin);

    if (dataWriter.eventCounter() % 1000 == 0)
      cout << "Event counter = " << dataWriter.eventCounter() << endl;

  }


}
