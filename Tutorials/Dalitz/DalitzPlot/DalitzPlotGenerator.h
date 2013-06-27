#if !(defined DALITZPLOTGENERATOR)
#define DALITZPLOTGENERATOR

using namespace std;

#include "IUAmpTools/PlotGenerator.h"

class AmpToolsInterface;
class Kinematics;

class DalitzPlotGenerator : public PlotGenerator
{
    
public:
    
  DalitzPlotGenerator( const FitResults& results );

  enum { 
    khm12 = 0, khm13, khm23, 
    kNumHists 
  };
    
private:
        
  void projectEvent( Kinematics* kin );

};

#endif
