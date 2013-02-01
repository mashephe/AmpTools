#if !defined(AMPTOOLSINTERFACEMPI)
#define AMPTOOLSINTERFACEMPI

#include "IUAmpTools/AmpToolsInterface.h"

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
