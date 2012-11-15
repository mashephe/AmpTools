#if !(defined ACCMCCOMPONENT)
#define ACCMCCOMPONENT

#include <string>

#include "TFile.h"

#include "AmpPlotter/PlotComponent.h"
#include "IUAmpTools/PlotGenerator.h"

class AccMCComponent : public PlotComponent
{
	
public:
    
	AccMCComponent( const string& title, 
                    const string& reaction,
                    PlotGenerator& pltGen );
    
	TH1F deliverPlot( unsigned int plotIndex,
                      double scale = 1.0 );

private:

};

#endif
