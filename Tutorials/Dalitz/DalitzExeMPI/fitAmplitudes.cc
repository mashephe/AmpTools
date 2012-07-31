
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <sys/time.h>

#include <mpi.h>

#include "MinuitInterface/MinuitMinimizationManager.h"

#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

// MPI enabled versions of the standard classes

#include "IUAmpToolsMPI/DataReaderMPI.h"
#include "IUAmpToolsMPI/NormIntInterfaceMPI.h"
#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpToolsMPI/LikelihoodCalculatorMPI.h"
#include "IUAmpToolsMPI/LikelihoodManagerMPI.h"

using std::complex;
using namespace std;
using namespace CLHEP;

typedef DataReaderMPI<DalitzDataReader> Reader;

int main( int argc, char* argv[] ){
	
  MPI_Init( &argc, &argv );
  
  int rank;
  int size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  // set default parameters
  
  string  configfile( "" );
  string  treeName( "nt" );
  
  // parse command line
  
  for (int i = 1; i < argc; i++){
    
    string arg(argv[i]);
    
    if (arg == "-c"){  
      if ((i+1 == argc) || (argv[i+1][0] == '-')) arg = "-h";
      else  configfile = argv[++i]; }
    if (arg == "-h"){
      cout << endl << " Usage for: " << argv[0] << endl << endl;
      cout << "\t -c <file>\t Config file" << endl;
      exit(1);}
  }
  
  if (configfile.size() == 0){
    cout << "No config file specified" << endl;
    exit(1);
  }
  
  ConfigFileParser parser(configfile);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();
  
  string fitname = cfgInfo->fitName();
  
  vector< ReactionInfo* > reactions = cfgInfo->reactionList();
  
  vector<AmplitudeManager*>        ampManagers;
  vector<NormIntInterfaceMPI*>     normInts;
  vector<string>                   normIntFiles;
  vector<Reader*>       dataReaders;
  vector<Reader*>       accMCReaders;
  vector<Reader*>       genMCReaders;
  vector<LikelihoodCalculatorMPI*> likCalcs;
  
  for( unsigned int i = 0; i < reactions.size(); i++){
    
    string fsName = reactions[i]->reactionName();

    vector<string> fsParticles  = reactions[i]->particleList();
    vector<string> dataFiles    = reactions[i]->data().second;
    vector<string> accMCFiles   = reactions[i]->accMC().second;
    vector<string> genMCFiles   = reactions[i]->genMC().second;

    ampManagers.push_back( new AmplitudeManager(fsParticles,fsName) );
    ampManagers[i]->registerAmplitudeFactor( BreitWigner() );
    ampManagers[i]->setupFromConfigurationInfo( cfgInfo );

    dataReaders.push_back( new Reader( dataFiles ) );
    accMCReaders.push_back( new Reader( accMCFiles ) );
    genMCReaders.push_back( new Reader( genMCFiles ) );
    
    normInts.push_back( new NormIntInterfaceMPI( genMCReaders[i],
                                                 accMCReaders[i],
                                                 *(ampManagers[i]) ) );
    normIntFiles.push_back( reactions[i]->normIntFile() );
  }
  
  // all of the above objects get setup on every node whether master or
  // worker -- although master/worker functionality is a little different
  // below here the execution begins to diverge
  
  MinuitMinimizationManager* fitManager;
  
  if( rank != 0 ){
    
    // this block of code runs on the worker nodes
  
    ParameterManagerMPI parManager( ampManagers );
    parManager.setupFromConfigurationInfo( cfgInfo );
  
    for ( unsigned int i = 0; i < ampManagers.size(); i++){
    
      likCalcs.push_back( new LikelihoodCalculatorMPI( *ampManagers[i], 
                                                       *normInts[i], 
                                                       *dataReaders[i],
                                                        parManager ) );
    }
    
    // put the worker nodes in the deliver likelihood loop
    // execution will hang here until LikelihoodCalculators are destroyed
    // on the master
    LikelihoodManagerMPI::deliverLikelihood();

    // now we are done with the likelihood calculators -- free the memory
    for (vector<LikelihoodCalculatorMPI*>::iterator itr = likCalcs.begin();
         itr != likCalcs.end(); itr++) { delete *itr; }

    // the workers need to help in final computation of the NI's
    for( unsigned int i = 0; i < normInts.size(); ++i ){
      
      normInts[i]->forceCacheUpdate();
    }      
    
  }
  else{
    
    // this block of code runs on the master node
    
    fitManager = new MinuitMinimizationManager(100);
    fitManager->setPrecision( 1E-13 );
  
    ParameterManagerMPI parManager( *fitManager, ampManagers );
    parManager.setupFromConfigurationInfo( cfgInfo );
    
    for ( unsigned int i = 0; i < ampManagers.size(); i++){
      
      likCalcs.push_back( new LikelihoodCalculatorMPI( *ampManagers[i], 
                                                       *normInts[i], 
                                                       *dataReaders[i],
                                                       parManager ) );
    }
    
    cout << "LIKELIHOODS BEFORE MINIMIZATION:  " << endl;
    for ( unsigned int i = 0; i < likCalcs.size(); i++){
      cout << "   -2 ln L:  " << (*(likCalcs[i]))() << endl;}
  
    fitManager->setStrategy( 1 );

    cout << "STARTING MINIMIZATION..." << endl;
  
    timeval tStart, tStop,tSpan;
    double dTime;
    gettimeofday( &(tStart), NULL );
    fitManager->migradMinimization();  
    gettimeofday( &(tStop), NULL );
    timersub( &(tStop), &(tStart), &tSpan );
    dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
    cout << "--->> Time for minimization after first iteration: " << dTime << " s" << endl; 

    cout << "LIKELIHOODS AFTER MINIMIZATION:  " << endl;
    for ( unsigned int i = 0; i < likCalcs.size(); i++){
    
      cout << "   -2 ln L:  " << (*(likCalcs[i]))() << endl;
    
      // go ahead and destroy these objects which will release the workers
      delete likCalcs[i];
    }
    
    cout << "WRITING FINAL NORMALIZATION INTEGRALS.." << endl;    
    for( unsigned int i = 0; i < normInts.size(); ++i ){
      
      normInts[i]->forceCacheUpdate();
      normInts[i]->exportNormIntCache( normIntFiles[i] );
    }
    
    if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
    
      cout << "ERROR: fit failed..." << endl;
      MPI_Finalize();
      return 1;
    }
  
    cout << "WRITING PARAMETERS..." << endl;
  
    string outputFile("fit.");  outputFile += fitname;  outputFile += ".txt";
    ofstream outFile(outputFile.c_str());
    parManager.writeParameters( outFile );
  
   }
  
  // cleanup allocated memory common on both master and worker nodes
  for (vector<AmplitudeManager*>::iterator itr = ampManagers.begin();
       itr != ampManagers.end(); itr++){ delete *itr; }
  for (vector<NormIntInterfaceMPI*>::iterator itr = normInts.begin();
       itr != normInts.end(); itr++) { delete *itr; }
  for (vector<Reader*>::iterator itr = dataReaders.begin();
       itr != dataReaders.end(); itr++) { delete *itr; }
  for (vector<Reader*>::iterator itr = accMCReaders.begin();
       itr != accMCReaders.end(); itr++) { delete *itr; }
  for (vector<Reader*>::iterator itr = genMCReaders.begin();
       itr != genMCReaders.end(); itr++) { delete *itr; }
  
  if( rank == 0 ) delete fitManager;
  
  MPI_Finalize();
  
  return 0;
}


