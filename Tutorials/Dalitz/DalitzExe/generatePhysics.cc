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

#include "IUAmpTools/report.h"
static const char* kModule = "generatePhysics";


using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************


  if (argc <= 3){
    report( NOTICE, kModule ) << "Usage:" << endl;
    report( NOTICE, kModule ) << "\tgeneratePhysics <config file name> <output file name> <number of events>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Generating Events According to Amplitudes *** " << endl << endl;


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  string outfilename(argv[2]);
  int nevents(atoi(argv[3]));

  report( INFO, kModule ) << "Config file name = " << cfgname << endl << endl;
  report( INFO, kModule ) << "Output file name = " << outfilename << endl;
  report( INFO, kModule ) << "Number of events = " << nevents << endl << endl;


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
  report( DEBUG, kModule ) << "Creating an AmpToolsInterface..." << endl;
  AmpToolsInterface::registerAmplitude(BreitWigner());
  AmpToolsInterface::registerDataReader(DalitzDataReader());
  AmpToolsInterface ATI(cfgInfo,AmpToolsInterface::kMCGeneration);
  report( DEBUG, kModule ) << "... Finished creating AmpToolsInterface" << endl;


    // ************************
    // create a DataWriter
    // ************************

  report( DEBUG, kModule ) << "Creating a Data Writer..." << endl;
  DalitzDataWriter dataWriter(outfilename);
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
    // first generate a set of phase space events
    // ************************

  report( DEBUG, kModule ) << "generating phase space..." << endl;

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

  report( DEBUG, kModule ) << "... finished generating phase space" << endl;

    
    // ************************
    // calculate intensities for all events
    // ************************

  report( DEBUG, kModule ) << "calculating intensities..." << endl;
  double maxIntensity = 1.2 * ATI.processEvents(reaction->reactionName());
  report( DEBUG, kModule ) << "... finished calculating intensities" << endl;


    // ************************
    // loop over all events again and do accept/reject
    // ************************

  report( DEBUG, kModule ) << "doing accept/reject..." << endl;
  for (int i = 0; i < nevents; i++){
    double intensity = ATI.intensity(i); 
    double rndm = drand48() * maxIntensity;
    if (intensity > rndm){
      Kinematics* kin = ATI.kinematics(i);
      dataWriter.writeEvent(*kin);
      delete kin;
    }
  }
  report( DEBUG, kModule ) << "... finished doing accept/reject" << endl;

  report( INFO, kModule ) << "KEPT " << dataWriter.eventCounter() << " events" << endl;

}
