#include <iostream>
#include <string>
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

using namespace std;

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Printing Amplitudes *** " << endl << endl;

  if (argc <= 2){
    cout << "Usage:" << endl << endl;
    cout << "\tprintAmplitudes <config file name> <input file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);
  string infilename(argv[2]);

  cout << "Config file name = " << cfgname << endl;
  cout << "Input file name  = " << infilename << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();


    // ************************
    // loop over reactions
    // ************************

  vector<ReactionInfo*> reactions = cfgInfo->reactionList();

  for (unsigned int irct = 0; irct < reactions.size(); irct++){


      // ************************
      // create an AmplitudeManager
      // ************************

    cout << endl << endl;
    cout << "Creating AmplitudeManager for reaction " << reactions[irct]->reactionName() << endl;
    AmplitudeManager ampMan(reactions[irct]->particleList(),reactions[irct]->reactionName());
    ampMan.registerAmplitudeFactor( BreitWigner() );
    ampMan.setupFromConfigurationInfo( cfgInfo );
    cout << "... Finished creating AmplitudeManager" << endl;


      // ************************
      // create a DataReader
      // ************************

    cout << "Creating Data Reader..." << endl;
    //DalitzDataReader dataReader(reactions[irct]->dataFiles());
    vector<string> args;
    args.push_back(infilename);
    DalitzDataReader dataReader(args);
    cout << "... Finished creating Data Reader" << endl << endl << endl;


      // ************************
      // print out detailed information for ten events
      // ************************

    cout << "**********************************************" << endl;
    cout << "  AMPLITUDES FOR REACTION " << reactions[irct]->reactionName() << endl;
    cout << "**********************************************" << endl << endl;

    for (unsigned int ievt = 0; ievt < 10; ievt++){
      Kinematics* kin = dataReader.getEvent();
      vector<HepLorentzVector> momenta = kin->particleList();

      cout << "  +++++++++++++++++++++++++++++++++" << endl;
      cout << "    EVENT NUMBER " << dataReader.eventCounter() << endl;
      //cout << "      from file " << reactions[irct]->dataFiles()[0] << endl;
      cout << "      from file " << infilename << endl;
      for (unsigned int imom = 0; imom < momenta.size(); imom++){
        cout << "        " << reactions[irct]->particleList()[imom] << endl;
        cout << "          E  = " << momenta[imom].e() << endl;
        cout << "          Px = " << momenta[imom].px() << endl;
        cout << "          Py = " << momenta[imom].py() << endl;
        cout << "          Pz = " << momenta[imom].pz() << endl;
      }
      cout << "  +++++++++++++++++++++++++++++++++" << endl << endl;

      vector<string> ampNames = ampMan.getAmpNames();
      for (unsigned int iamp = 0; iamp < ampNames.size(); iamp++){
        cout << "    ----------------------------------" << endl;
        cout << "      AMPLITUDE = " << ampNames[iamp] << endl;
        cout << "    ----------------------------------" << endl << endl;
        vector< const Amplitude* > ampFactors = ampMan.getFactors(ampNames[iamp]);
        vector <vector <int> > permutations = ampMan.getPermutations(ampNames[iamp]);
        for (unsigned int iperm = 0; iperm < permutations.size(); iperm++){
          cout << "        PERMUTATION = ";
          for (unsigned int ipar = 0; ipar < permutations[iperm].size(); ipar++){
            cout << permutations[iperm][ipar] << " ";
          }
          cout << endl << endl;
          for (unsigned int ifact = 0; ifact < ampFactors.size(); ifact++){
            cout << "          AMPLITUDE FACTOR = " << ampFactors[ifact]->name() << endl;
            cout << "          RESULT = " 
                 << ampFactors[ifact]->calcAmplitude(kin,permutations[iperm]) << endl << endl;
          }
        }
      }
      cout << "      ---------------------------------" << endl;
      cout << "        CALCULATING INTENSITY" << endl;
      cout << "      ---------------------------------" << endl << endl;
      double intensity = ampMan.calcIntensity(kin);
      cout << endl << "          INTENSITY = " << intensity << endl << endl << endl;

    }

  }

}
