#if !defined(GAMMAKKAMP)
#define GAMMAKKAMP

#include "IUAmpTools/Amplitude.h"

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
// 2. J_gamma: total ang. mom. of photon                                  //
// 3. M      : z-projection of J/psi                                      //
// 4. lambda : photon helicity                                            //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class Kinematics;

class gammaKKAmp : public Amplitude
{
    
public:
	
        gammaKKAmp() : Amplitude() { setDefaultStatus( true ); }
	gammaKKAmp( int j, int J_gamma, int M, int l ); 
	
	string name() const { return "gammaKKAmp"; }
    
	complex< GDouble > calcAmplitude( GDouble** pKin ) const;
	void printEvent( GDouble** pKin ) const;
	
	gammaKKAmp* newAmplitude( const vector< string >& args ) const;
	gammaKKAmp* clone() const;

private:
        
  int m_j;
  int m_J_gamma;
  int m_M;
  int m_l;

};

#endif
