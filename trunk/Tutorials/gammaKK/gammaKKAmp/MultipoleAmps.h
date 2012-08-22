#if !defined(MULTIPOLEAMPS)
#define MULTIPOLEAMPS

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

// j,mu are the total and z projection of the spin of X
// J_gamma is the angular momentum related to the radiative transition
// M is the polarization of the parent (J/psi)
// l is the helicity of the photon

class Kinematics;

class MultipoleAmps : public UserAmplitude<MultipoleAmps>
{
    
public:
	
 MultipoleAmps() : UserAmplitude<MultipoleAmps>() {}
  MultipoleAmps(const vector<string> &args); 
  ~MultipoleAmps(){}

  string name() const { return "MultipoleAmps"; }
    
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  // void printEvent( GDouble** pKin ) const;
	
private:
        
  int m_j;
  int m_J_gamma;
  int m_M;
  int m_lambda;

};

#endif
