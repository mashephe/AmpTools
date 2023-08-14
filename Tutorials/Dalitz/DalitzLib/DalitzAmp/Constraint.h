#if !defined(CONSTRAINT)
#define CONSTRAINT

#include "IUAmpTools/Neg2LnLikContrib.h"
#include "IUAmpTools/UserNeg2LnLikContrib.h"
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

class Constraint : public UserNeg2LnLikContrib< Constraint >{

public:
  
  Constraint() : UserNeg2LnLikContrib< Constraint >() { }

  Constraint( const vector< string >& args );

  ~Constraint(){}

  string name() const { return "Constraint"; }

  double calcLHContribution( double x ) const;

  double neg2LnLikelihood();

  static double drawThis(double *x, double *par);
  
private:

  AmpParameter a, m0, g, bkg;
  TH1F *hData;
  TF1 *func;
	
};

#endif
