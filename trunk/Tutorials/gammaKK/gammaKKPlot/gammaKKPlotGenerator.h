#if !(defined GAMMAKKPLOTGENERATOR)
#define GAMMAKKPLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"

class FitResults;
class Kinematics;

using namespace std;

class gammaKKPlotGenerator : public PlotGenerator
{
    
public:
    
  gammaKKPlotGenerator( const FitResults& results );

  enum { 
    khm12 = 0, khm13, khm23,
    // angles of photon in CM frame
    khPhotonCosTheta, khPhotonPhi,
    // angles of 1st K in KK CM frame
    khKaonCosTheta, khKaonPhi,
    kNumHists
  };
    
private:
        
  void projectEvent( Kinematics* kin );

};

#endif
