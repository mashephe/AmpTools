
#include <iostream>
#include <string>

#include "TH1.h"

#include "AmpPlotter/AccMCComponent.h"
#include "IUAmpTools/Histogram.h"

using namespace std;

AccMCComponent::AccMCComponent( const string& title, 
                                const string& reaction,
                                PlotGenerator& pltGen ) :
PlotComponent( title, reaction, pltGen )
{
  setDataOK();
}

TH1F
AccMCComponent::deliverPlot( unsigned int plotIndex, 
                             double scale ) {
  
  Histogram hist = generator().projection( plotIndex, reaction(), 
                                           PlotGenerator::kAccMC );
	
  TH1F plot = hist.toRoot();
  
	plot.SetFillStyle( fillStyle() );
  plot.SetFillColor( fillColor() );
  plot.Scale( scale );
  plot.SetLineColor( fillColor() );
  
	return plot;
}
