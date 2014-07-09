
#include "DalitzPlot/DalitzPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

DalitzPlotGenerator::DalitzPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  bookHistogram( khm12, "Mass( 1 2 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( khm13, "Mass( 1 3 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( khm23, "Mass( 2 3 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( kdltz, "Dalitz Plot", new Histogram2D( 60, 0.0, 9.0, 60, 0.0, 9.0 ) );
}

void
DalitzPlotGenerator::projectEvent( Kinematics* kin ){
          
  HepLorentzVector P1 = kin->particle(0);
  HepLorentzVector P2 = kin->particle(1);
  HepLorentzVector P3 = kin->particle(2);
      
  fillHistogram( khm12, (P1+P2).m() );
  fillHistogram( khm13, (P1+P3).m() );  
  fillHistogram( khm23, (P2+P3).m() );
  fillHistogram( kdltz, (P1+P2).m2(), (P2+P3).m2() );
}