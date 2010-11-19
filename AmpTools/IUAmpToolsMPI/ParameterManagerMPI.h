#if !defined(PARAMETERMANAGERMPI)
#define PARAMETERMANAGERMPI

#include <vector>
#include <map>
#include <complex>

#include "IUAmpTools/ParameterManager.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "MinuitInterface/MIObserver.h"

using namespace std;

class ParameterManagerMPI : public ParameterManager
{

 public:

  enum { kMaxNameLength = 200 };

  // since there is only one MinuitMinimizationManager we must add
  // an additional constructor for the worker nodes

  // this constructor is called on the head node
  ParameterManagerMPI( MinuitMinimizationManager& minuitManager,
		       AmplitudeManager* ampManager );
  ParameterManagerMPI( MinuitMinimizationManager& minuitManager,
		       const vector< AmplitudeManager* >& ampManagers );

  // this constructor should be called on the worker nodes
  ParameterManagerMPI( AmplitudeManager* ampManager );
  ParameterManagerMPI( const vector< AmplitudeManager* >& ampManagers );

  ~ParameterManagerMPI();

  // override this function to do the appropriate thing depending
  // on whether this instance is on the head or worker node
  void addProductionParameter( const string& ampName );

  // the likelihood calculator will need to call this routine in order
  // to update the parameters in advance of the likelihood calculator
  void updateParameters();

 private:

  void setupMPI();

  vector<AmplitudeManager*> m_ampManagers;

  int m_rank;
  int m_numProc;
  bool m_isMaster;

  map< string, complex< double >* > m_parMap;

};

#endif
