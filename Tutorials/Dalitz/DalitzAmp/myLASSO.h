#if !defined(MYLASSO)
#define MYLASSO

#include "IUAmpTools/LHContribution.h"
#include "IUAmpTools/UserLHContribution.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>

using std::complex;
using namespace std;

class myLASSO : public UserLHContribution< myLASSO >{

public:
  
  myLASSO() : UserLHContribution< myLASSO >() { }

  myLASSO( const vector< string >& args );

  ~myLASSO(){}

  string name() const { return "myLASSO"; }

  double calcLHContribution( double x ) const;

  double neg2LnLikelihood();

  static double drawThis(double *x, double *par);
  
private:

  double lambda;
  int npars;
  AmpParameter amp_re[10], amp_im[10];
	
};

#endif
