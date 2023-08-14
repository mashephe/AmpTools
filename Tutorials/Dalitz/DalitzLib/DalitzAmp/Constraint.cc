#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <fstream>

#include "DalitzAmp/Constraint.h"

Constraint::Constraint( const vector< string >& args ) :
UserNeg2LnLikContrib< Constraint >(args)
{
  assert(args.size()==4);
  a = AmpParameter(args[0]);
  m0 = AmpParameter(args[1]);
  g = AmpParameter(args[2]);
  bkg = AmpParameter(args[3]);
  registerParameter(a);
  registerParameter(m0);
  registerParameter(g);
  registerParameter(bkg);
  
  /** in principle, the idea is to read in some external (binned) data here
  *   here, we don't have any and thus just generate some toy data
  */
  func = new TF1("func",Constraint::drawThis,0,2,4);
  func->SetParameters(1,1.0,0.2,0.5);
  hData = new TH1F("hData","",100,0,2);
  for(int ev=0;ev<5000;ev++)
  	hData->Fill(func->GetRandom());

  hData->SaveAs("external_data.root");
}

double
Constraint::calcLHContribution( double x ) const {
  double s = pow(x,2);
  complex<double> amp = complex<double>(a,0.0) / complex<double>(s-m0*m0,m0*g);
  return norm(amp) + bkg;
}

double Constraint::drawThis(double *x, double *par){
  double s = pow(x[0],2);
  complex<double> amp = complex<double>(par[0],0.0) / complex<double>(s-par[1]*par[1],par[1]*par[2]);
  return norm(amp) + par[3];
}

double
Constraint::neg2LnLikelihood(){
  double chi2=0;
  for(int i=1;i<=hData->GetNbinsX();i++){
    double _x = hData->GetBinCenter(i);
    double _y = hData->GetBinContent(i);
    double _ey = hData->GetBinError(i);
    if(_ey)
    		chi2 += pow(calcLHContribution(_x)-_y,2)/pow(_ey,2);
  }
  return chi2;
}
