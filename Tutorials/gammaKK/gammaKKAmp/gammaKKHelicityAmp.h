#if !defined(GAMMAKKHELICITYAMP)
#define GAMMAKKHELICITYAMP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"

#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// 2012/07/27 Kei Moriya                                                  //
//                                                                        //
// Main amplitude to be used in J/psi -> gamma KK analysis.               //
// The amplitude will take in 4 numbers, which are                        //
// 1. J      : resonance spin of KK                                       //
// 2. mu     : helicity of KK resonance                                   //
// 3. M      : z-projection of J/psi                                      //
// 4. lambda : photon helicity                                            //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class Kinematics;

class gammaKKHelicityAmp : public UserAmplitude<gammaKKHelicityAmp>
{
    
public:
 gammaKKHelicityAmp() : UserAmplitude<gammaKKHelicityAmp>() {}
  gammaKKHelicityAmp(const vector<string> &args); 
  
  string name() const { return "gammaKKHelicityAmp"; }
  
  complex< GDouble > calcAmplitude( GDouble** pKin ) const;
  void printEvent( GDouble** pKin ) const;
  
 private:
  int m_j;
  int m_mu;
  int m_M;
  int m_lambda;

};

#endif
