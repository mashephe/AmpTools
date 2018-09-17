#if !defined(BREITWIGNER)
#define BREITWIGNER

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

#ifdef GPU_ACCELERATION

void BreitWigner_exec( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                          GDouble mass, GDouble width,
                          int daughter1, int daughter2 );

#endif // GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

class BreitWigner : public UserAmplitude< BreitWigner >{

public:
  
  BreitWigner() : UserAmplitude< BreitWigner >() { }

  BreitWigner( const vector< string >& args );

  ~BreitWigner(){}

  string name() const { return "BreitWigner"; }

  //  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userData ) const;

  enum UserVars { kMass2 = 0, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  void calcUserData( GDouble** pKin, GDouble* userData ) const;
  
#ifdef GPU_ACCELERATION

  void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const;

  bool isGPUEnabled() const { return true; }

#endif // GPU_ACCELERATION
  
private:
	
  AmpParameter m_mass;
  AmpParameter m_width;
  int m_daughter1;
  int m_daughter2;  

};

#endif
