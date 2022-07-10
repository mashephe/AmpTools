#include <cstdlib>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzDataIO/DalitzDataWriter.h"

using namespace std;

#include "IUAmpTools/report.h"
static const char* kModule = "toyAcceptance";

int main(int argc, char** argv){


    // ************************
    // usage
    // ************************

  if (argc <= 2){
    report( NOTICE, kModule ) << "Usage:" << endl << endl;
    report( NOTICE, kModule ) << "\ttoyAcceptance <input file name> <output file name>" << endl << endl;
    return 0;
  }

  report( INFO, kModule ) << endl << " *** Simulating Detector Effects *** " << endl << endl;

    // ************************
    // parse the command line parameters
    // ************************

  string infilename(argv[1]);
  string outfilename(argv[2]);

  report( INFO, kModule ) << "Input file name  = " << infilename << endl;
  report( INFO, kModule ) << "Output file name = " << outfilename << endl;


    // ************************
    // create a DataReader
    // ************************

  report( DEBUG, kModule ) << "Creating a Data Reader..." << endl;
  vector<string> args;
  args.push_back(infilename);
  DalitzDataReader dataReader(args);
  report( DEBUG, kModule ) << "... Finished creating a Data Reader" << endl << endl;


    // ************************
    // create a DataWriter
    // ************************

  report( DEBUG, kModule ) << "Creating a Data Writer..." << endl;
  DalitzDataWriter dataWriter(outfilename);
  report( DEBUG, kModule ) << "... Finished creating a Data Writer" << endl << endl;


    // ************************
    // simulate acceptance
    // ************************

  Kinematics* kin;

  while( (kin = dataReader.getEvent()) ){

    vector<TLorentzVector> pList = kin->particleList();

    double m12 = (pList[0]+pList[1]).M();

    double efficiency = 0.1 + (0.9/3.0)*m12;

    double rndm = drand48();

    if (rndm < efficiency) dataWriter.writeEvent(*kin);

    if (dataReader.eventCounter() % 1000 == 0)
      report( INFO, kModule ) << "Event counter = " << dataReader.eventCounter() << endl;

    delete kin;

  }


}


