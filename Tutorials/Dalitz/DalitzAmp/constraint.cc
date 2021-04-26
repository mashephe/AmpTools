#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <fstream>

#include "DalitzAmp/constraint.h"

constraint::constraint( const vector< string >& args ) :
UserLHContribution< constraint >(args)
{

  ifstream datafile(args[0].c_str());
  local_data = BinnedData();
  string line;

  while(datafile.is_open() && getline(datafile,line)){
    istringstream ss(line);
    string s1, s2, s3;
    ss >> s1 >> s2 >> s3;
    local_data.addPoint(atof(s1.c_str()),atof(s2.c_str()),atof(s3.c_str()));
  }

  a = AmpParameter(args[1]);
  m0 = AmpParameter(args[2]);
  g = AmpParameter(args[3]);
  registerParameter(a);
  registerParameter(m0);
  registerParameter(g);

}

double
constraint::calcLHContribution( double x ) const {
  double s = pow(x,2);
  complex<double> amp = complex<double>(a,0.0) / complex<double>(s-m0*m0,m0*g);
  return norm(amp);
}

double
constraint::neg2LnLikelihood(){
  double chi2=0;
  for(unsigned int i=0;i<local_data.getN();i++){
    datapoint pnt = local_data.getPoint(i);
    chi2 += pow(calcLHContribution(pnt.x)-pnt.y,2)/pow(pnt.eyl,2);
  }
  return chi2;
}