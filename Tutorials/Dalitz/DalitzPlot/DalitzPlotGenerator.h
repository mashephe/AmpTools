#if !(defined DALITZPLOTGENERATOR)
#define DALITZPLOTGENERATOR

#include "IUAmpTools/PlotGenerator.h"

class FitResults;
class Kinematics;

class DalitzPlotGenerator : public PlotGenerator
{
    
public:
    
  DalitzPlotGenerator( const FitResults& results );

  enum { 
    khm12 = 0, khm13, khm23, kdltz,
    kNumHists 
  };
    
private:
        
  void projectEvent( Kinematics* kin );

};

#endif
