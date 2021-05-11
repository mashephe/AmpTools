#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <fstream>

#include "DalitzAmp/myLASSO.h"

myLASSO::myLASSO( const vector< string >& args ) :
UserLHContribution< myLASSO >(args)
{
  npars = 0;
  lambda = atof(args[0].c_str());
  for(unsigned int i=1;i<args.size()/2+1;i++){
	  amp_re[i-1] = AmpParameter(args[2*i-1]);
	  amp_im[i-1] = AmpParameter(args[2*i]);
	  registerParameter(amp_re[i-1]);
	  registerParameter(amp_im[i-1]);
	  npars++;
  }

}

double
myLASSO::calcLHContribution( double x ) const {
	return 0;
}


double
myLASSO::neg2LnLikelihood(){
	double val=0;
	for(int i=0;i<npars;i++){
		val += lambda*(abs(amp_re[i])+abs(amp_im[i]));
	}
	return val;
}
