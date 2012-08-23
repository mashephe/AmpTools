#include <iostream>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "gammaKKDataIO/gammaKKDataWriter.h"
#include "gammaKKExe/constants.h"

// amplitudes
#include "gammaKKAmp/gammaKKHelicityAmp.h"
#include "gammaKKAmp/MultipoleAmps.h"
#include "gammaKKAmp/NBodyPhaseSpaceFactory.h"

#include "TRandom2.h"

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Generating Events According to Amplitudes *** " << endl << endl;

  if (argc <= 4){
    cout << "Usage:" << endl << endl;
    cout << "\tgeneratePhysics <config file name> <output file name> <number of events> <random seed>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  string outfilename(argv[2]);
  int nevents(atoi(argv[3]));
  int seed(atoi(argv[4]));

  cout << "Config file name = " << cfgname << endl << endl;
  cout << "Output file name = " << outfilename << endl;
  cout << "Number of events = " << nevents << endl << endl;
  cout << "Random seed      = " << seed << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  ReactionInfo* reaction = cfgInfo->reactionList()[0];


    // ************************
    // create an AmplitudeManager
    // ************************

  cout << endl << endl;
  cout << "Creating AmplitudeManager for reaction " << reaction->reactionName() << endl;
  AmplitudeManager ampMan(reaction->particleList(),reaction->reactionName());

  // Register all amplitudes that we could use
  ampMan.registerAmplitudeFactor( gammaKKHelicityAmp() );
  ampMan.registerAmplitudeFactor( MultipoleAmps() );
  ampMan.setupFromConfigurationInfo( cfgInfo );
  cout << "... Finished creating AmplitudeManager" << endl;


    // ************************
    // create a DataWriter
    // ************************

  cout << "Creating Data Writer..." << endl;
  gammaKKDataWriter dataWriter(outfilename);
  cout << "... Finished creating Data Writer" << endl << endl;


    // ************************
    // set up an NBodyPhaseSpaceFactory object
    // ************************

  double parentMass = m_Jpsi;
  vector<double> daughterMasses;
  daughterMasses.push_back(m_photon);
  daughterMasses.push_back(m_KZero);
  daughterMasses.push_back(m_KZero);

  NBodyPhaseSpaceFactory generator(parentMass, daughterMasses);

    // ************************
    // Need random generator for accept/reject
    // ************************
  TRandom2 *rand = new TRandom2(seed);

    // ************************
    // first generate a set of phase space events
    // ************************

  cout << "generating phase space..." << endl;

  AmpVecs packedSummary;

  for (int i = 0; i < nevents; i++){

    vector<HepLorentzVector> fourmomenta = generator.generateDecay();

    Kinematics* kin = new Kinematics(fourmomenta);

    packedSummary.loadEvent(kin,i,nevents);

    delete kin;

  }

  packedSummary.allocateAmps(ampMan,true);

  cout << "... finished generating phase space" << endl;

    
    // ************************
    // calculate intensities for all events
    // ************************

  cout << "calculating intensities..." << endl;
  double maxIntensity = 1.2 * ampMan.calcIntensities(packedSummary);
  cout << "... finished calculating intensities" << endl;


    // ************************
    // loop over all events again and do accept/reject
    // ************************

  cout << "doing accept/reject..." << endl;
  for (int i = 0; i < nevents; i++){
    double intensity = packedSummary.m_pdIntensity[i]; 
    double rndm = rand->Rndm() * maxIntensity;
    if (intensity > rndm){
      Kinematics* kin = packedSummary.getEvent(i);
      dataWriter.writeEvent(*kin);
      delete kin;
    }
  }
  cout << "... finished doing accept/reject" << endl;

  cout << "KEPT " << dataWriter.eventCounter() << " events" << endl;



}
