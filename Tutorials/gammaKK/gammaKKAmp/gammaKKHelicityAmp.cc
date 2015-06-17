
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "TMath.h"

#include "IUAmpTools/Kinematics.h"
#include "gammaKKHelicityAmp.h"
#include "clebschGordan.h"
#include "wignerD.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

// returns normalization for spin J
GDouble norm(int J);

gammaKKHelicityAmp::gammaKKHelicityAmp(const vector<string> &args) :
  UserAmplitude<gammaKKHelicityAmp>(){

  m_j       = atoi(args[0].c_str()); // total J of resonance
  m_mu      = atoi(args[1].c_str()); // helicity of ressonance J
  m_M       = atoi(args[2].c_str()); // polarization of parent (J/psi)
  m_lambda  = atoi(args[3].c_str()); // helicity of photon

  assert(abs(m_mu) <= m_j);
  assert( m_M == 1 || m_M == -1 );
  assert( m_lambda == 1 || m_lambda == -1 );
}

complex< GDouble >
gammaKKHelicityAmp::calcAmplitude( GDouble** pKin ) const {
  
  TLorentzVector p4gamma ( pKin[0][1], pKin [0][2], pKin[0][3], pKin[0][0] );
  TLorentzVector p4K1    ( pKin[1][1], pKin [1][2], pKin[1][3], pKin[1][0] );
  TLorentzVector p4K2    ( pKin[2][1], pKin [2][2], pKin[2][3], pKin[2][0] );

  // In helicity amplitudes, the decay amplitude for a given resonance
  // with quantum numbers (j,mu) will be
  // N_J D(J,M,mu-lambda)(dir. of KK) N_j D(j,mu,0)(dir. of K)

  // All 4-vectors are in J/psi rest frame.
  // First we need to get the angles of the KK system.
  TLorentzVector p4KK = p4K1 + p4K2;
  // wignerD_1 takes in theta variable in degrees
  GDouble thetaKK = p4KK.Theta() * 180. / TMath::Pi();
  GDouble phiKK   = p4KK.Phi();

  // Next we calculate the angles of K1 in the KK rest frame.
  // Boost into the KK rest frame.
  TLorentzRotation resRestBoost(-p4KK.BoostVector());
  TLorentzVector p4K1_res = resRestBoost * p4K1;
  GDouble thetaK1 = p4K1_res.Theta() * 180. / TMath::Pi();
  GDouble phiK1   = p4K1_res.Phi();

  // Normalization factors for J/psi (J=1) and KK system (spin j)
  GDouble coeff = norm(1) * norm(m_j);

  complex<GDouble > amp( 0, 0 );

  amp = coeff * wignerD_1(1,m_M,m_mu-m_lambda,thetaKK,phiKK) * wignerD_1(m_j,m_mu,0,thetaK1,phiK1);

  return amp;
}

// returns normalization for spin J
GDouble norm(int J){
  return sqrt((2. * J + 1.) / (4. * TMath::Pi()));
}
