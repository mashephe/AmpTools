#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"

// We don't want to rely on MCToolkit for now...
// The ROOT class TGenPhaseSpace is equivalent.
// #include "MCToolkit/NBodyPhaseSpaceFactory.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom2.h"

#include "gammaKKDataIO/gammaKKDataWriter.h"
#include "gammaKKExe/constants.h"

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Generating Phase Space *** " << endl << endl;

  if (argc <= 2){
    cout << "Usage:" << endl << endl;
    cout << "\tgeneratePhaseSpace <output file name> <number of events>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string filename(argv[1]);
  const int NEVENTS(atoi(argv[2]));

  cout << "Output file name = " << filename << endl;
  cout << "Number of events = " << NEVENTS << endl << endl;


    // ************************
    // set up a gammaKKDataWriter object
    // ************************

  cout << "Creating Data Writer..." << endl;

  gammaKKDataWriter dataWriter(filename);

  cout << "... Finished creating Data Writer" << endl << endl;


    // ************************
    // set up an TGenPhaseSpace object (from ROOT)
    // ************************

  TGenPhaseSpace *generator = new TGenPhaseSpace();
  // This takes in the TLorentzVector of the parent particle,
  // the number of daughters, and an array of masses for the daughters.
  // The option "fermi" can be specified to get fermi motion...
  TLorentzVector p4_Jpsi(0,0,0,m_Jpsi);
  double daughterMasses[3] = {m_photon, m_KZero, m_KZero};
  generator->SetDecay(p4_Jpsi, 3, daughterMasses);
  double maxWeight = generator->GetWtMax();

  // Now generate each event by using TGenPhaseSpace::Generate(),
  // and use function GetDecay(n),
  // which returns the TLorentzVector* for the nth particle.

  // We will need a random number to weight the events properly.
  // The way TGenPhaseSpace works, each event is assigned a weight
  // at generation time. To have a completely flat phase space,
  // the event needs to be filled with that weight.
  //
  // Since we want each weight to have equal weight (=1), we
  // circumvent this by throwing a random number between
  // 0 and the maximum allowed weight (= maxWeight) and
  // keep the event only if the weight is higher.
  // This essentially gives each event a probability
  // of being kept proportional to the weight.

  TRandom2* rand = new TRandom2();
  int nRejected = 0;

  for(int i=0;i<NEVENTS;i++){

    vector<HepLorentzVector> fourvectors(3);

    double weight = generator->Generate();

    if(weight < maxWeight * rand->Rndm()){
      // If the event is not kept, decrement the event counter
      // and go on to the next event.
      i--;
      nRejected++;
      continue;
    }

    TLorentzVector p4gamma(*(generator->GetDecay(0)));
    TLorentzVector p4Ks1(*(generator->GetDecay(1)));
    TLorentzVector p4Ks2(*(generator->GetDecay(2)));

    // Set up vector of 4-vectors
    fourvectors[0].setX(p4gamma.X());
    fourvectors[0].setY(p4gamma.Y());
    fourvectors[0].setZ(p4gamma.Z());
    fourvectors[0].setE(p4gamma.E());

    fourvectors[1].setX(p4Ks1.X());
    fourvectors[1].setY(p4Ks1.Y());
    fourvectors[1].setZ(p4Ks1.Z());
    fourvectors[1].setE(p4Ks1.E());

    fourvectors[2].setX(p4Ks2.X());
    fourvectors[2].setY(p4Ks2.Y());
    fourvectors[2].setZ(p4Ks2.Z());
    fourvectors[2].setE(p4Ks2.E());

    // Create Kinematics object, which holds the vector of 4-vectors
    Kinematics kin(fourvectors);

    // Write out the event to file
    dataWriter.writeEvent(kin);
    if(dataWriter.eventCounter() % 100000==0)
      cout << "Event counter = " << dataWriter.eventCounter() << endl;
  }

  cout << "Total of " << nRejected << " events rejected (" << NEVENTS << " events generated)" << endl;
  cout << "Rejection rate: " << 1. * nRejected / NEVENTS * 100. << " %" << endl;

  // This taken out since it relies on the MCToolkit/NBodyPhaseSpace class
  // double parentMass = 3.0;
  // vector<double> daughterMasses;
  // daughterMasses.push_back(0.2);
  // daughterMasses.push_back(0.2);
  // daughterMasses.push_back(0.2);
  // 
  // NBodyPhaseSpaceFactory generator(parentMass, daughterMasses);
  // // ************************
  // // use the gammaKKDataWriter object to write events to a file
  // // ************************
  // 
  // for (int i = 0; i < NEVENTS; i++){
  // 
  //   vector<HepLorentzVector> fourvectors = generator.generateDecay();
  // 
  //   Kinematics kin(fourvectors);
  // 
  //   dataWriter.writeEvent(kin);
  // 
  //   if (dataWriter.eventCounter() % 1000 == 0)
  //     cout << "Event counter = " << dataWriter.eventCounter() << endl;
  // 
  // }
}
