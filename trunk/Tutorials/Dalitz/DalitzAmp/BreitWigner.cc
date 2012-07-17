#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include "CLHEP/Vector/LorentzVector.h"
#include "IUAmpTools/Kinematics.h"
#include "DalitzAmp/BreitWigner.h"

BreitWigner::BreitWigner( const AmpParameter& mass0, 
                          const AmpParameter& width0,
                          int daughter1, int daughter2) :
Amplitude(),
m_mass0( mass0 ),
m_width0( width0 ),
m_daughter1(daughter1),
m_daughter2(daughter2)
{
  
  // this is not the default constructor
  setDefaultStatus( false );
  
  // need to register any free parameters so the framework knows about them
  registerParameter( m_mass0 );
  registerParameter( m_width0 );
  
}


complex< GDouble >
BreitWigner::calcAmplitude( GDouble** pKin ) const {

  HepLorentzVector P1(pKin[m_daughter1-1][1], pKin[m_daughter1-1][2],
                      pKin[m_daughter1-1][3], pKin[m_daughter1-1][0]);

  HepLorentzVector P2(pKin[m_daughter2-1][1], pKin[m_daughter2-1][2],
                      pKin[m_daughter2-1][3], pKin[m_daughter2-1][0]);

  GDouble m = (P1+P2).m();

  complex<GDouble> bwdenominator(m*m - m_mass0*m_mass0, m_mass0*m_width0);

  return  complex<GDouble>(1.0,0.0) / bwdenominator;

}


BreitWigner*
BreitWigner::newAmplitude( const vector< string >& args ) const {

  assert( args.size() == 4 );

  AmpParameter mass(args[0]);
  AmpParameter width(args[1]);
  int daughter1(atoi(args[2].c_str()));
  int daughter2(atoi(args[3].c_str()));

  return new BreitWigner( mass, width, daughter1, daughter2 );

}


BreitWigner*
BreitWigner::clone() const {

  return ( isDefault() ? new BreitWigner() : 
    new BreitWigner( m_mass0, m_width0, m_daughter1, m_daughter2 ) );

}


#ifdef GPU_ACCELERATION
void
BreitWigner::launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
  
  // use integers to endcode the string of daughters -- one index in each
  // decimal place
  
  GPUBreitWigner_exec( dimGrid,  dimBlock, GPU_AMP_ARGS, 
                       m_mass0, m_width0, m_daughter1, m_daughter2 );

}
#endif //GPU_ACCELERATION

