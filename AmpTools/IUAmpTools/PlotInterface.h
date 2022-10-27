#ifndef PLOTINTERFACE
#define PLOTINTERFACE

//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
// 
// Copyright Trustees of Indiana University 2010, all rights reserved
// 
// This software written by Matthew Shepherd, Ryan Mitchell, and 
//                  Hrayr Matevosyan at Indiana University, Bloomington
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// 3. Neither the name of the University nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// 
// Creation of derivative forms of this software for commercial
// utilization may be subject to restriction; written permission may be
// obtained from the Trustees of Indiana University.
// 
// INDIANA UNIVERSITY AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES, 
// EXPRESS OR IMPLIED.  By way of example, but not limitation, INDIANA 
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR 
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR 
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS, 
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be 
// held liable for any liability with respect to any claim by the user or 
// any other party arising from use of the program.
//******************************************************************************

#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/FitResults.h"
#include "IUAmpTools/DataReader.h"

using namespace std;

class PlotInterface{

  public:
    PlotInterface(const char *out = "output.root");
    virtual ~PlotInterface() {
      clear();
    }
    static void registerAmplitude( const Amplitude& defaultAmplitude);
    static void registerDataReader( const DataReader& defaultDataReader);
    AmpToolsInterface *ATI(){return m_ati;}
    void setupFromCFG(const char* cfgfile, bool useMC);
    void setupFromFitResult(const char* fitresult);
    void plot(string reaction);
    void writeOutput();
    void addToDraw(string reaction, string name);

  protected:
    void clear();
    pair<float,complex<float>> calcAmp(int iEvent, string reaction, string name);

  private:
    AmpToolsInterface* m_ati;
    FitResults* m_fitResults;
    ConfigurationInfo* m_configurationInfo;
    ConfigFileParser* m_parser;
    map<string,DataReader*> m_reader;
    map<string,TGenPhaseSpace*> m_gen;
    map<int,double> m_daughters;
    double m_mother;
    bool m_useMC;
    map<string,TTree*> m_tree;
    float m_PxP[10];
    float m_PyP[10];
    float m_PzP[10];
    float m_EnP[10];
    float m_weight;
    float m_weightOrig;
    TFile *fOut;
    map<string,string> m_amps_to_draw;
    map<string,float> m_amp_weights_re;
    map<string,float> m_amp_weights_im;
    map<string,float> m_amp_weights_INT;
};

#endif