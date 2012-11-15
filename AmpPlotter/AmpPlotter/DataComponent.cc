
#include <iostream>
#include <string>

#include "TH1.h"

#include "AmpPlotter/DataComponent.h"
#include "IUAmpTools/Histogram.h"

using namespace std;

DataComponent::DataComponent( const string& title, 
                              const string& reaction,
                              PlotGenerator& pltGen ) :
PlotComponent( title, reaction, pltGen )
{
    setDataOK();
}

TH1F
DataComponent::deliverPlot( unsigned int plotIndex, 
                            double scale ) {
			
    Histogram hist = generator().projection( plotIndex, reaction(), 
                                             PlotGenerator::kData );
	
    TH1F plot = hist.toRoot();
    	
	plot.SetMarkerStyle( markerStyle() );
	plot.SetMarkerSize( markerSize() );
    plot.SetFillStyle( fillStyle() );
    plot.Scale( scale );
    
	return plot;
}
