#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "gammaKKExe/constants.h"
#include "gammaKKPlot/gammaKKPlotGenerator.h"
#include "gammaKKAmp/gammaKKHelicityAmp.h"
#include "gammaKKAmp/MultipoleAmps.h"

// Include ROOT 2D histograms currently impossible,
// as gammaKKPlotGenerator::fillProjections
// will return vector<Histogram>, and this will
// be used in PlotGenerator???
// #include "TH2F.h"

gammaKKPlotGenerator::gammaKKPlotGenerator( AmpToolsInterface& ati ) :
PlotGenerator( ati )
{
  bookHistogram( khm12, "Mass( 1 2 )", Histogram( 60, 0.0, 3.0 ) );
  bookHistogram( khm13, "Mass( 1 3 )", Histogram( 60, 0.0, 3.0 ) );
  bookHistogram( khm23, "Mass( 2 3 )", Histogram( 60, 0.0, 3.0 ) );
  bookHistogram( khPhotonCosTheta, "cos#theta_{#gamma}", Histogram( 100, -1.0, 1.0 ) );
  bookHistogram( khPhotonPhi, "#phi_{#gamma}", Histogram( 100, -_PI, _PI ) );
  bookHistogram( khKaonCosTheta, "cos#theta_{K_{1}}", Histogram( 100, -1.0, 1.0 ) );
  bookHistogram( khKaonPhi, "#phi_{K_{1}}", Histogram( 100, -_PI, _PI ) );
}
                
void 
gammaKKPlotGenerator::projectEvent( Kinematics* kinematics ){
      
  // gamma      
  HepLorentzVector P1 = kinematics->particle(0);
  // Kaon 1
  HepLorentzVector P2 = kinematics->particle(1);
  // Kaon 2
  HepLorentzVector P3 = kinematics->particle(2);
  
  fillHistogram( khm12, (P1+P2).m() );
  fillHistogram( khm13, (P1+P2).m() );
  fillHistogram( khm23, (P2+P3).m() );

  float costhetaPhoton = P1.cosTheta();
  float phiPhoton = P1.phi();
    
  fillHistogram( khPhotonCosTheta, costhetaPhoton );
  fillHistogram( khPhotonPhi, phiPhoton );
  
  // Need to boost to KK rest frame to get angles of Kaon      
  HepLorentzVector resonance = P2 + P3;
  HepLorentzRotation resRestBoost( -resonance.boostVector() );
  HepLorentzVector P2_KKrest = resRestBoost * P2;
  
  float costhetaKaon = P2_KKrest.cosTheta();
  float phiKaon = P2_KKrest.phi();
  
  fillHistogram( khKaonCosTheta, costhetaKaon );
  fillHistogram( khKaonPhi, phiKaon );
}
