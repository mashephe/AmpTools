
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "TMath.h"

#include "IUAmpTools/Kinematics.h"
#include "MultipoleAmps.h"
#include "./clebschGordan.h"
#include "./wignerD.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// returns normalization for spin J
// (need to rename this as having multiple definitions
//  in different amplitude classes will make conflicting libraries)
GDouble norm_multipoleAmps(int J);

MultipoleAmps::MultipoleAmps(const vector<string> &args) :
  UserAmplitude<MultipoleAmps>()
{
  m_j       = atoi(args[0].c_str()); // total J of resonance
  m_J_gamma = atoi(args[1].c_str()); // J from radiative multipole transition
  m_M       = atoi(args[2].c_str()); // polarization of parent (J/psi)
  m_lambda  = atoi(args[3].c_str()); // helicity of photon

  assert( m_M == 1 || m_M == -1 );
  assert( m_lambda == 1 || m_lambda == -1 );

  if(m_J_gamma==0){
    cout << "Total angular momentum of photon CANNOT BE 0!!!" << endl;
    cout << "aborting..." << endl;
    assert(false);
  }

}

complex< GDouble >
MultipoleAmps::calcAmplitude( GDouble** pKin ) const {

  HepLorentzVector p4gamma ( pKin[0][1], pKin [0][2], pKin[0][3], pKin[0][0] );
  HepLorentzVector p4K1    ( pKin[1][1], pKin [1][2], pKin[1][3], pKin[1][0] );
  HepLorentzVector p4K2    ( pKin[2][1], pKin [2][2], pKin[2][3], pKin[2][0] );

  // In helicity amplitudes, the decay amplitude for a given resonance
  // with quantum numbers (j,mu) will be
  // N_J D(J,M,mu-lambda)(dir. of KK) N_j D(j,mu,0)(dir. of K)

  // All 4-vectors are in J/psi rest frame.
  // First we need to get the angles of the KK system.
  HepLorentzVector p4KK = p4K1 + p4K2;
  // wignerD_1 takes in theta variable in degrees
  GDouble thetaKK = p4KK.theta() * 180. / TMath::Pi();
  GDouble phiKK   = p4KK.phi();

  // Next we calculate the angles of K1 in the KK rest frame.
  // Boost into the KK rest frame.
  HepLorentzRotation resRestBoost(-p4KK.boostVector());
  HepLorentzVector p4K1_res = resRestBoost * p4K1;
  GDouble thetaK1 = p4K1_res.theta() * 180. / TMath::Pi();
  GDouble phiK1   = p4K1_res.phi();

  // Normalization factors for J/psi (J=1) and KK system (spin j)
  GDouble coeff = norm_multipoleAmps(1) * norm_multipoleAmps(m_J_gamma);

  // Parity relates the lambda = +1 and -1 states by +1 or -1.
  // If we are thinking of a radiative transition that
  // changes the parity, and furthermore, if the final
  // state has positive parity, the negative helicity
  // state acquires a phase of (-1)^(J_gamma - 1)
  if(m_lambda == -1 && (m_J_gamma - 1)%2==1) coeff *= -1.;

  complex<GDouble > amp( 0, 0 );

  // We will need to sum over the different helicities with
  // the correct Clebsch-Gordan to get the multipole amplitude.
  for(int mu=-m_j;mu<=m_j;mu++){
    amp += coeff
      * wignerD_1(1,m_M,mu-m_lambda,thetaKK,phiKK)             // Wigner D for decay of KK system of spin j
      * wignerD_1(m_j,mu,0,thetaK1,phiK1)                      // Wigner D for decay of ccbar of state |JM>
      * clebschGordan( m_J_gamma, m_j, -m_lambda, mu, 1, mu - m_lambda);   // CGC for this helicity amplitude that goes into multipole
  }
   
  return amp;
}

// returns normalization for spin J
GDouble norm_multipoleAmps(int J){
  return sqrt((2. * J + 1.) / (4. * TMath::Pi()));
}
