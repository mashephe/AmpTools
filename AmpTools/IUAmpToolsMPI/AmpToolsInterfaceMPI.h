#if !defined(AMPTOOLSINTERFACEMPI)
#define AMPTOOLSINTERFACEMPI

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <fstream>

#include <mpi.h>
#include <pthread.h>

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/LikelihoodCalculator.h"

#include "IUAmpTools/AmpToolsInterface.h"

#include "IUAmpToolsMPI/DataReaderMPI.h"
#include "IUAmpToolsMPI/ParameterManagerMPI.h"
#include "IUAmpToolsMPI/LikelihoodCalculatorMPI.h"
#include "IUAmpToolsMPI/NormIntInterfaceMPI.h"


using namespace std;


class AmpToolsInterfaceMPI : public AmpToolsInterface{

  public:

    AmpToolsInterfaceMPI(ConfigurationInfo* cfgInfo);

    void finalizeFit();

  private:

    int m_rank;
    int m_numProc;

};



#endif
