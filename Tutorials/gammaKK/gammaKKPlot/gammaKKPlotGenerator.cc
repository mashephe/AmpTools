#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "IUAmpTools/Histogram1D.h"

#include "gammaKKExe/constants.h"
#include "gammaKKPlot/gammaKKPlotGenerator.h"

// Include ROOT 2D histograms currently impossible,
// as gammaKKPlotGenerator::fillProjections
// will return vector<Histogram>, and this will
// be used in PlotGenerator???
// #include "TH2F.h"

gammaKKPlotGenerator::gammaKKPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
  bookHistogram( khm12, "Mass( 1 2 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( khm13, "Mass( 1 3 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( khm23, "Mass( 2 3 )", new Histogram1D( 60, 0.0, 3.0 ) );
  bookHistogram( khPhotonCosTheta, "cos#theta_{#gamma}", new Histogram1D( 100, -1.0, 1.0 ) );
  bookHistogram( khPhotonPhi, "#phi_{#gamma}", new Histogram1D( 100, -_PI, _PI ) );
  bookHistogram( khKaonCosTheta, "cos#theta_{K_{1}}", new Histogram1D( 100, -1.0, 1.0 ) );
  bookHistogram( khKaonPhi, "#phi_{K_{1}}", new Histogram1D( 100, -_PI, _PI ) );
}
                
void 
gammaKKPlotGenerator::projectEvent( Kinematics* kinematics ){
      
  // gamma      
  TLorentzVector P1 = kinematics->particle(0);
  // Kaon 1
  TLorentzVector P2 = kinematics->particle(1);
  // Kaon 2
  TLorentzVector P3 = kinematics->particle(2);
  
  fillHistogram( khm12, (P1+P2).M() );
  fillHistogram( khm13, (P1+P2).M() );
  fillHistogram( khm23, (P2+P3).M() );

  float costhetaPhoton = P1.CosTheta();
  float phiPhoton = P1.Phi();
    
  fillHistogram( khPhotonCosTheta, costhetaPhoton );
  fillHistogram( khPhotonPhi, phiPhoton );
  
  // Need to boost to KK rest frame to get angles of Kaon      
  TLorentzVector resonance = P2 + P3;
  TLorentzRotation resRestBoost( -resonance.BoostVector() );
  TLorentzVector P2_KKrest = resRestBoost * P2;
  
  float costhetaKaon = P2_KKrest.CosTheta();
  float phiKaon = P2_KKrest.Phi();
  
  fillHistogram( khKaonCosTheta, costhetaKaon );
  fillHistogram( khKaonPhi, phiKaon );
}
