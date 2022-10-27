#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "IUAmpTools/PlotInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"

#include "IUAmpTools/report.h"
static const char* kModule = "myPlotter";

int main(int argc, char** argv){

  string cfgfile = argv[1];
  PlotInterface *plotter = new PlotInterface("out_plotter.root");
  PlotInterface::registerDataReader(DalitzDataReader());
  PlotInterface::registerAmplitude(BreitWigner());

  plotter->setupFromCFG(cfgfile.c_str(),1);
  plotter->addToDraw("dalitz","R12");
  plotter->addToDraw("dalitz","R13");

  plotter->plot("dalitz");
  //plotter->addAmplitude("name");
  plotter->writeOutput();
  return 1;
}