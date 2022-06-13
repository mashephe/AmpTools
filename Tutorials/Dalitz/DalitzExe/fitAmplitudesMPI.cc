#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <mpi.h>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpToolsMPI/AmpToolsInterfaceMPI.h"
#include "IUAmpToolsMPI/DataReaderMPI.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

#include "IUAmpTools/report.h"
static const char* kModule = "fitAmplitudesMPI";

using std::complex;
using namespace std;

int main( int argc, char* argv[] ){

  MPI_Init( &argc, &argv );
  
  int rank;
  int size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );


    // ************************
    // usage
    // ************************


  if (argc <= 1){
    report( INFO, kModule ) << "Usage:" << endl << endl;
    report( INFO, kModule ) << "\tfitAmplitudesMPI <config file name>" << endl << endl;
    MPI_Finalize();
    return 0;
  }

  report( INFO, kModule ) << " *** Performing the Fit *** " << endl;


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  report( INFO, kModule ) << "Config file name:  " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();


    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterfaceMPI::registerAmplitude(BreitWigner());
  AmpToolsInterfaceMPI::registerDataReader(DataReaderMPI<DalitzDataReader>());

  AmpToolsInterfaceMPI ATI(cfgInfo);

  if (rank == 0){
    
    double neg2LL = ATI.likelihood();
    report( INFO, kModule ) << "-2 ln(L) BEFORE MINIMIZATION:  " << neg2LL << endl;

    MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
    fitManager->setStrategy(1);

    fitManager->migradMinimization();

    if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
      report( WARNING, kModule ) << "Fit failed." << endl;
    }

    report( INFO, kModule ) << "-2 ln(L) AFTER MINIMIZATION:  " << ATI.likelihood() << endl;

    ATI.finalizeFit();
  }

  ATI.exitMPI();
  MPI_Finalize();

  return 0;

}


