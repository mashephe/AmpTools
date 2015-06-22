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

  if (rank == 0) cout << endl << " *** Performing the Fit *** " << endl << endl;

  if (argc <= 1){
    if (rank == 0) cout << "Usage:" << endl << endl;
    if (rank == 0) cout << "\tfitAmplitudesMPI <config file name>" << endl << endl;
    MPI_Finalize();
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  if (rank == 0) cout << "Config file name = " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  if (rank == 0) cfgInfo->display();


    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterfaceMPI::registerAmplitude(BreitWigner());
  AmpToolsInterfaceMPI::registerDataReader(DataReaderMPI<DalitzDataReader>());

  AmpToolsInterfaceMPI ATI(cfgInfo);

  if (rank == 0){
    cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ATI.likelihood() << endl;

    MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
    fitManager->setPrecision(1E-13);
    fitManager->setStrategy(1);

    fitManager->migradMinimization();

    if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
      cout << "ERROR: fit failed..." << endl;
    }

    cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ATI.likelihood() << endl;
  }

  ATI.finalizeFit();

  MPI_Finalize();

  return 0;

}


