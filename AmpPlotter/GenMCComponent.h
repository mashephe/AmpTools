#if !(defined GENMCCOMPONENT)
#define GENMCCOMPONENT

#include <string>

#include "TFile.h"

#include "AmpPlotter/PlotComponent.h"
#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class GenMCComponent : public PlotComponent
{
	
public:
    
    GenMCComponent(  const string& title, 
                                     const string& reaction,
                                     PlotGenerator& pltGen );
    
	TH1F deliverPlot( unsigned int plotIndex,
                      double scale = 1.0 );
    
private:

};

#endif
