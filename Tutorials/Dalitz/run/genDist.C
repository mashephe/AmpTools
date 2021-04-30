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
  
  TCanvas *c = new TCanvas();
  h->Draw("e1");
  h->SetMarkerStyle(23);
  h->SetMarkerSize(2);
  TF1 *f2 = new TF1("f2",bw,0,2,3);
  f2->SetParameters(5,1.0,0.2);
  f2->FixParameter(0,5.09092*sqrt(sqrt(2)));
  f2->FixParameter(1,1.0033);
  f2->FixParameter(2,0.200786);
  h->Fit(f2,"R");
  f2->Draw("same");
  c->SaveAs("dummyfit.png");
  
  exit(0);
}
