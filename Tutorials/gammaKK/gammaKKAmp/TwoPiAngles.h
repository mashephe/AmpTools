#if !defined(TWOPIANGLES)
#define TWOPIANGLES

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

class TwoPiAngles : public UserAmplitude<TwoPiAngles>
{
    
public:
	
 TwoPiAngles() : UserAmplitude<TwoPiAngles>() {}
  TwoPiAngles(const vector<string> &args); 
  ~TwoPiAngles(){}

  string name() const { return "TwoPiAngles"; }
    
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  // void printEvent( GDouble** pKin ) const;
	
private:
        
  int m_j;
  int m_J_gamma;
  int m_M;
  int m_l;

};

#endif
