#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include <fstream>
#include <complex>
#include "TCanvas.h"

double bw(double *x, double *par){

  double s = pow(x[0],2);
  double m0 = par[1];
  double g = par[2];
  double a = par[0];
  complex<double> amp = complex<double>(a,0.0) / complex<double>(s-m0*m0,m0*g);
  return norm(amp);
}

void genDist(){

  TF1 *f = new TF1("f",bw,0,2,3);
  f->SetParameters(1,1.000,0.200);

  TH1F *h = new TH1F("h","",50,0,2);

  for(int i=0;i<5000;i++){
    h->Fill(f->GetRandom());
  }

  ofstream out("extradata.dat");
  for(int i=1;i<=h->GetNbinsX();i++){
    out << h->GetBinCenter(i) << " " << h->GetBinContent(i) << " " << (h->GetBinError(i)<1?1:h->GetBinError(i)) << endl;
  }
  out.close();
  exit(0);
}
