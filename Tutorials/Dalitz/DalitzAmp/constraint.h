#if !defined(CONSTRAINT)
#define CONSTRAINT

#include "IUAmpTools/LHContribution.h"
#include "IUAmpTools/UserLHContribution.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/BinnedData.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>

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

  //AmpParameter par;

  
private:

  AmpParameter a, m0, g;
  BinnedData local_data;
	
};

#endif
