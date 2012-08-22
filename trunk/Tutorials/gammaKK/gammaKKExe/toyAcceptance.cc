#include <iostream>
#include <string>
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "gammaKKDataIO/gammaKKDataReader.h"
#include "gammaKKDataIO/gammaKKDataWriter.h"

#include "TRandom2.h"

using namespace std;

/**
 *
 * This class will calculate the pseudo-acceptance of each
 * event based on the event 4-vectors.
 *
 * The function calcEfficiency takes in a Kinematics object
 * from the DataReader, and calculates the efficiency.
 * Modify this function to whatever acceptance is required.
 *
 * A random number is thrown using the ROOT TRandom2 class,
 * and the event is accepted proportinately to the acceptance
 * calculated by calcEfficiency.
 *
 * The accepted events are written out with a DataWriter
 * to a ROOT file.
 *
 */

double calcEfficiency(Kinematics *kin);

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Simulating Detector Effects *** " << endl << endl;

  if (argc <= 2){
    cout << "Usage:" << endl << endl;
    cout << "\ttoyAcceptance <input file name> <output file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  vector<string> infilenames;
  infilenames.push_back(argv[1]);

  string outfilename(argv[2]);

  cout << "Input file name  = " << infilenames[0] << endl;
  cout << "Output file name = " << outfilename << endl;


    // ************************
    // create a DataReader
    // ************************

  cout << "Creating Data Reader..." << endl;
  gammaKKDataReader dataReader(infilenames);
  cout << "... Finished creating Data Reader" << endl << endl;


    // ************************
    // create a DataWriter
    // ************************

  cout << "Creating Data Writer..." << endl;
  gammaKKDataWriter dataWriter(outfilename);
  cout << "... Finished creating Data Writer" << endl << endl;


    // ************************
    // simulate acceptance
    // ************************

  Kinematics* kin;

  // Use ROOT TRandom2 class for random generator
  TRandom2* rand = new TRandom2();

  while (kin = dataReader.getEvent()){

    double efficiency = calcEfficiency(kin);

    double rndm = rand->Rndm();

    if (rndm < efficiency) dataWriter.writeEvent(*kin);

    if (dataReader.eventCounter() % 100000 == 0)
      cout << "Event counter = " << dataReader.eventCounter() << endl;

    delete kin;

  }

  return 0;
}

double calcEfficiency(Kinematics *kin){

  // Calculate the efficiency of this event based on the
  // Kinematics object, which contains the 4-vectors of
  // all the particles.

  double acceptance = 1.;

  vector<HepLorentzVector> pList = kin->particleList();
  for(int n=0;n<pList.size();n++){
    acceptance *= (1. - 0.05 * pow(cos(pList[n].theta()),2.));
    acceptance *= (1. - 0.03 * pow(sin(pList[n].phi()),2.));
  }
  
  return acceptance;
}

