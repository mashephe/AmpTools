#if !defined(AMPTOOLSINTERFACEMPI)
#define AMPTOOLSINTERFACEMPI

#include "IUAmpTools/AmpToolsInterface.h"

using namespace std;

class AmpToolsInterfaceMPI : public AmpToolsInterface{
  
public:
  
  AmpToolsInterfaceMPI(ConfigurationInfo* cfgInfo);
  ~AmpToolsInterfaceMPI(){}
  
  void finalizeFit( const string& tag = "" );

  // exit MPI should be called on the leader process before
  // MPI_Finalize() or variables go out of scope
  void exitMPI();
  
private:
  
  int m_rank;
  int m_numProc;
  
};



#endif
