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

  complex< GDouble > calcAmplitude( GDouble** pKin, GDouble* userVars ) const;
  
  // **********************
  // The following lines are optional and can be used to precalcualte
  // user-defined data that the amplitudes depend on.
  
  // Use this for indexing a user-defined data array and notifying
  // the framework of the number of user-defined variables.
  enum UserVars { kMass2 = 0, kNumUserVars };
  unsigned int numUserVars() const { return kNumUserVars; }
  
  // This function needs to be defined -- see comments and discussion
  // in the .cc file.
  void calcUserVars( GDouble** pKin, GDouble* userVars ) const;
  
  // This is an optional addition if the calcAmplitude routine
  // can run with only the user-defined data and not the original
  // four-vectors.  It is used to optimize memory usage in GPU
  // based fits.
  //bool needsUserVarsOnly() const { return true; }
  bool needsUserVarsOnly() const { return false; }
  // **  end of optional lines **
  
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
