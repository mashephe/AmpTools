#if !defined(LIKELIHOODMANAGERMPI)
#define LIKELIHOODMANAGERMPI

#include <map>

class LikelihoodCalculatorMPI;

using namespace std;

class LikelihoodManagerMPI
{

 public:

  enum FitCommand { kComputeLikelihood,
    kUpdateParameters,
    kExit };

  LikelihoodManagerMPI(){};

  static void registerCalculator( int id, LikelihoodCalculatorMPI* calc );
  static void deliverLikelihood();

 private:

  static void setupMPI();

  static bool m_mpiSetup;
  static bool m_isMaster;

  static map< int, LikelihoodCalculatorMPI* > m_calcMap;
};

#endif
