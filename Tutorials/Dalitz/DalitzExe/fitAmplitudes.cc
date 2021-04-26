#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzAmp/constraint.h"
#include "IUAmpTools/FitResults.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"


using std::complex;
using namespace std;

int main( int argc, char* argv[] ){


    // ************************
    // usage
    // ************************

  cout << endl << " *** Performing the Fit *** " << endl << endl;

  if (argc <= 1){
    cout << "Usage:" << endl << endl;
    cout << "\tfitAmplitudes <config file name>" << endl << endl;
    return 0;
  }


    // ************************
    // parse the command line parameters
    // ************************

  string cfgname(argv[1]);

  cout << "Config file name = " << cfgname << endl << endl;


    // ************************
    // parse the config file
    // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();


    // ************************
    // AmpToolsInterface
    // ************************

  AmpToolsInterface::registerAmplitude(BreitWigner());
  AmpToolsInterface::registerLHContribution(constraint());
  AmpToolsInterface::registerDataReader(DalitzDataReader());

  AmpToolsInterface ATI(cfgInfo);

  cout << "LIKELIHOOD BEFORE MINIMIZATION:  " << ATI.likelihood() << endl;

  MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
  fitManager->setPrecision(1E-13);
  fitManager->setStrategy(1);

  fitManager->migradMinimization();

  if( fitManager->status() != 0 && fitManager->eMatrixStatus() != 3 ){
    cout << "ERROR: fit failed..." << endl;
  }

  cout << "LIKELIHOOD AFTER MINIMIZATION:  " << ATI.likelihood() << endl;


  ifstream data("extradata.dat");
  string line;
  TGraphErrors *gr = new TGraphErrors();
  gr->SetMarkerStyle(23);
  gr->SetMarkerSize(2);
  int c=0;
  while(data.is_open() && getline(data,line)){
    stringstream ss(line);
    string s1, s2, s3;
    ss >> s1 >> s2 >> s3;
    gr->SetPoint(c,atof(s1.c_str()),atof(s2.c_str()));
    gr->SetPointError(c,0,atof(s3.c_str()));
    c++;
  }
  ATI.finalizeFit();

  FitResults results("dalitz1.fit");

  TF1 *f = new TF1("f",constraint::drawThis,0,2,3);
  f->SetParameters(results.parValue("am"),results.parValue("m1"),results.parValue("g1"));

  TCanvas *can = new TCanvas();
  gr->Draw("ap");
  f->Draw("same");
  can->SaveAs("test.pdf");


  return 0;

}


