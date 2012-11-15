
#include <iostream>
#include <string>

#include "TH1.h"

#include "AmpPlotter/GenMCComponent.h"
#include "IUAmpTools/Histogram.h"

using namespace std;

GenMCComponent::GenMCComponent( const string& title, 
                                const string& reaction,
                                PlotGenerator& pltGen ) :
PlotComponent( title, reaction, pltGen )
{
    setDataOK();
}

TH1F
GenMCComponent::deliverPlot( unsigned int plotIndex, 
                             double scale ) {
    
    Histogram hist = generator().projection( plotIndex, reaction(), 
                                             PlotGenerator::kGenMC );
    
    TH1F plot = hist.toRoot();	
    
	plot.SetFillStyle( fillStyle() );
    plot.SetFillColor( fillColor() );
    plot.Scale( scale );
    plot.SetLineColor( fillColor() );

	return plot;
}
