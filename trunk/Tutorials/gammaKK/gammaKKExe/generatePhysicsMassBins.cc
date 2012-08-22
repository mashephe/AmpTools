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
#include "gammaKKAmp/TwoPiAngles.h"
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

  if (argc != 6){
    cout << "Usage:" << endl << endl;
    cout << "\tgeneratePhysics <config file name> <output file name base> <number of events> <random seed> <number of mass bins>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  const string outfilenamebase(argv[2]);
  const int nevents(atoi(argv[3]));
  const int seed(atoi(argv[4]));
  const int NBINS(atoi(argv[5]));


  cout << "Config file name      = " << cfgname << endl << endl;
  cout << "Output file name base = " << outfilenamebase << endl;
  cout << "Number of events      = " << nevents << endl << endl;
  cout << "Random seed           = " << seed << endl << endl;
  cout << "Number of mass bins   = " << NBINS << endl;

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
  ampMan.registerAmplitudeFactor( TwoPiAngles() );
  ampMan.registerAmplitudeFactor( gammaKKHelicityAmp() );
  ampMan.registerAmplitudeFactor( MultipoleAmps() );
  ampMan.setupFromConfigurationInfo( cfgInfo );
  cout << "... Finished creating AmplitudeManager" << endl;


    // ************************
    // create multiple DataWriters
    // ************************

  cout << "Creating Data Writers..." << endl;
  gammaKKDataWriter *dataWriter[NBINS];
  char filename[200];
  for(int i=0;i<NBINS;i++){
    sprintf(filename,"%s_%d.root",outfilenamebase.c_str(),i);
    dataWriter[i] = new gammaKKDataWriter(filename);
  }
  cout << "... Finished creating Data Writers" << endl << endl;


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

  // Calculate the minimum and maximum values possible for KK mass
  const double MIN = 2. * m_KZero;
  const double MAX = m_Jpsi;
  const double binning = (MAX - MIN) / NBINS;

  for (int i = 0; i < nevents; i++){
    double intensity = packedSummary.m_pdIntensity[i]; 
    double rndm = rand->Rndm() * maxIntensity;
    if (intensity > rndm){
      Kinematics* kin = packedSummary.getEvent(i);

      // Calculate the dataWriter bin for this event based on
      // M23
      HepLorentzVector p2 = kin->particle(1);
      HepLorentzVector p3 = kin->particle(2);
      double m23 = (p2 + p3).m();
      int bin = static_cast<int>( floor((m23 - MIN) / binning));

      dataWriter[bin]->writeEvent(*kin);
      delete kin;
    }
  }
  cout << "... finished doing accept/reject" << endl;

  // For having multiple dataWriter objects,
  // it seems it is necessary to call the write of the ROOT file
  // and TTree explicitly, instead of having the destructor
  // take care of it.
  int totalEvents = 0;
  for(int i=0;i<NBINS;i++){
    cout << "KEPT " << dataWriter[i]->eventCounter() << " events in bin " << i << endl;
    dataWriter[i]->write();
    totalEvents += dataWriter[i]->eventCounter();
  }
  cout << "KEPT " << totalEvents << " TOTAL events" << endl;



}
