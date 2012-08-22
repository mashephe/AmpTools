
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "IUAmpTools/Kinematics.h"
#include "TwoPiAngles.h"
#include "./clebschGordan.h"
#include "./wignerD.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/ThreeVector.h"

TwoPiAngles::TwoPiAngles(const vector<string> &args) :
  UserAmplitude<TwoPiAngles>()
{
  m_j       = atoi(args[0].c_str()); // total J of resonance
  m_J_gamma = atoi(args[1].c_str()); // J from radiative multipole transition
  m_M       = atoi(args[2].c_str()); // polarization of parent (J/psi)
  m_l       = atoi(args[3].c_str()); // helicity of photon

  assert( m_M == 1 || m_M == -1 );
  assert( m_l == 1 || m_l == -1 );

  if(m_J_gamma==0){
    cout << "Total angular momentum of photon CANNOT BE 0!!!" << endl;
    cout << "aborting..." << endl;
    assert(false);
  }

}

complex< GDouble >
TwoPiAngles::calcAmplitude( GDouble** pKin ) const {

  HepLorentzVector gamma ( pKin[0][1], pKin [0][2], pKin[0][3], pKin[0][0] );
  HepLorentzVector p1    ( pKin[1][1], pKin [1][2], pKin[1][3], pKin[1][0] );
  HepLorentzVector p2    ( pKin[2][1], pKin [2][2], pKin[2][3], pKin[2][0] );

  GDouble cosTheta_g = -gamma.vect().cosTheta();
  GDouble phi_g = gamma.vect().phi();

  HepLorentzVector resonance = p1 + p2;

  HepLorentzRotation resRestBoost( -resonance.boostVector() );

  HepLorentzVector p1_res = resRestBoost * p1;

  // phi is generated randomly, we just need to pick a direction?
  Hep3Vector z = -gamma.vect().unit();
  //  Hep3Vector y = -p1.vect().cross(z).unit();
  Hep3Vector bdir ( 0.0, 0.0, 1.0 );
  Hep3Vector y = bdir.cross(z).unit();
  Hep3Vector x = y.cross(z).unit();

  Hep3Vector angles( (p1_res.vect()).dot(x),
		     (p1_res.vect()).dot(y),
		     (p1_res.vect()).dot(z) );

  GDouble cosTheta = angles.cosTheta();
  GDouble phi = angles.phi();

  // we have to pass in theta_1 instead of cos(theta_1) because 
  // cos(-theta) = cos(theta) - we will lose the minus sign difference on theta
  
  GDouble theta_1 = acos( cosTheta ) * 180.0 / 3.14159265;

  // normalization factor

  GDouble coeff = ( sqrt( 2. * m_J_gamma + 1 ) * sqrt ( 2. * m_j + 1 ) / 
		    ( sqrt(2) * 4 * 3.14159265 ) );

  if( m_l == -1 ){
    if( m_J_gamma == 2 || m_J_gamma == 4 ){
      coeff = coeff*(-1);
    }
  }

  // sum over helicities of the resonance

  complex< GDouble > term( 0, 0 );

  for( int mu = -m_j; mu <= m_j; ++mu ){

    term += ( coeff * wignerD( 1, m_M, ( mu - m_l ), cosTheta_g, (phi_g+3.14159265), phi ) *
			     wignerD_1( m_j, mu, 0, theta_1, 0 ) *
	      clebschGordan( m_J_gamma, m_j, -m_l, mu, 1, ( mu - m_l ) ) );
    /*
    cout << "Event: " << endl;
    cout << "coef: " << coeff << endl;
    cout << "gamma: (" << pKin[0][1] << ", " << pKin[0][2] << ", " << pKin[0][3] << ", " << pKin[0][0] << ")" << endl;
    cout << "pi0_1: (" << pKin[1][1] << ", " << pKin[1][2] << ", " << pKin[1][3] << ", " << pKin[1][0] << ")" << endl;
    cout << "pi0_2: (" << pKin[2][1] << ", " << pKin[2][2] << ", " << pKin[2][3] << ", " << pKin[2][0] << ")" << endl;
    cout << "(j,mu,j_gamma,M,l): " <<  m_j << ", " << mu << ", " << m_J_gamma << ", " << m_M << ", " << m_l << endl;
    cout << "CG: " << clebschGordan( m_J_gamma, m_j, -m_l, mu, 1, ( mu - m_l ) ) << endl;
    cout << "D1: " << wignerD( 1, m_M, ( mu - m_l ), cosTheta_g, (phi_g+3.14159265), phi ) << endl;
    cout << "D2: " << wignerD( m_j, mu, 0, cosTheta, 0 ) << endl;
    */
  }

  //  cout << "Angles (cos(th_g), phi_g, cos(th_1), phi_1): " << cosTheta_g << ",  " << phi_g << ", " << cosTheta << ", " << phi << endl;
  //  cout << "TwoPiAngles: " << term << endl << endl;

  return term;

}
