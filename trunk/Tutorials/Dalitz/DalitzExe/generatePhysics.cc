#include <iostream>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzDataIO/DalitzDataWriter.h"
#include "DalitzAmp/BreitWigner.h"


using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Generating Events According to Amplitudes *** " << endl << endl;

  if (argc <= 3){
    cout << "Usage:" << endl << endl;
    cout << "\tgeneratePhysics <config file name> <output file name> <number of events>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  string outfilename(argv[2]);
  int nevents(atoi(argv[3]));

  cout << "Config file name = " << cfgname << endl << endl;
  cout << "Output file name = " << outfilename << endl;
  cout << "Number of events = " << nevents << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  ReactionInfo* reaction = cfgInfo->reactionList()[0];


    // ************************
    // create an AmpToolsInterface
    // ************************

  cout << endl << endl;
  cout << "Creating an AmpToolsInterface..." << endl;
  AmpToolsInterface::registerAmplitude(BreitWigner());
  AmpToolsInterface::registerDataReader(DalitzDataReader());
  AmpToolsInterface ATI(cfgInfo,AmpToolsInterface::kMCGeneration);
  cout << "... Finished creating AmpToolsInterface" << endl;


    // ************************
    // create a DataWriter
    // ************************

  cout << "Creating a Data Writer..." << endl;
  DalitzDataWriter dataWriter(outfilename);
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
    // first generate a set of phase space events
    // ************************

  cout << "generating phase space..." << endl;

  for (int i = 0; i < nevents; i++){

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
    Kinematics* kin = new Kinematics(fourvectors);


      // load events into the AmpToolsInterface

    ATI.loadEvent(kin,i,nevents);

    delete kin;

  }

  cout << "... finished generating phase space" << endl;

    
    // ************************
    // calculate intensities for all events
    // ************************

  cout << "calculating intensities..." << endl;
  double maxIntensity = 1.2 * ATI.processEvents(reaction->reactionName());
  cout << "... finished calculating intensities" << endl;


    // ************************
    // loop over all events again and do accept/reject
    // ************************

  cout << "doing accept/reject..." << endl;
  for (int i = 0; i < nevents; i++){
    double intensity = ATI.intensity(i); 
    double rndm = drand48() * maxIntensity;
    if (intensity > rndm){
      Kinematics* kin = ATI.kinematics(i);
      dataWriter.writeEvent(*kin);
      delete kin;
    }
  }
  cout << "... finished doing accept/reject" << endl;

  cout << "KEPT " << dataWriter.eventCounter() << " events" << endl;



}
