#if !defined(CONSTRAINT)
#define CONSTRAINT

#include "IUAmpTools/LHContribution.h"
#include "IUAmpTools/UserLHContribution.h"
#include "IUAmpTools/AmpParameter.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>

#include "TH1F.h"
#include "TF1.h"

using std::complex;
using namespace std;

class constraint : public UserLHContribution< constraint >{

public:
  
  constraint() : UserLHContribution< constraint >() { }

  constraint( const vector< string >& args );

  ~constraint(){}

  string name() const { return "constraint"; }

  double calcLHContribution( double x ) const;

  double neg2LnLikelihood();

  static double drawThis(double *x, double *par);
  
private:

  AmpParameter a, m0, g, bkg;
  TH1F *hData;
  TF1 *func;
	
};

#endif
