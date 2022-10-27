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


#include "IUAmpTools/PlotInterface.h"
#include "IUAmpTools/report.h"

static const char* kModule = "PlotInterface";

PlotInterface::PlotInterface(const char *out) :
m_ati(NULL),
m_configurationInfo(NULL),
m_fitResults(NULL),
m_parser(NULL)
{
  fOut = new TFile(out,"recreate");
  m_useMC = true;
}

void PlotInterface::registerAmplitude(const Amplitude& amplitude){
  AmpToolsInterface::registerAmplitude(amplitude);
}

void PlotInterface::registerDataReader(const DataReader& reader){
  AmpToolsInterface::registerDataReader(reader);
}

void PlotInterface::addToDraw(string reaction, string name){
  vector<AmplitudeInfo*> ampList = m_configurationInfo->amplitudeList(reaction,"","");
  string s_index = reaction+name;
  m_amps_to_draw[s_index] = name;
  m_amp_weights_re[s_index] = 1.;
  m_amp_weights_im[s_index] = 0.;
  m_amp_weights_INT[s_index] = 0.;

  fOut->cd();
  for(auto &p : m_amps_to_draw){
    string ampName = p.second;
    m_tree[reaction]->Branch(Form("weightRe_%s",ampName.c_str()),&m_amp_weights_re[ampName],Form("weightRe_%s/F",ampName.c_str()));
    m_tree[reaction]->Branch(Form("weightIm_%s",ampName.c_str()),&m_amp_weights_im[ampName],Form("weightIm_%s/F",ampName.c_str()));
    m_tree[reaction]->Branch(Form("weightINT_%s",ampName.c_str()),&m_amp_weights_INT[ampName],Form("weightINT_%s/F",ampName.c_str()));
  }
}

void PlotInterface::setupFromCFG(const char* cfgfile, bool useMC){
  m_parser = new ConfigFileParser(cfgfile);
  m_configurationInfo = m_parser->getConfigurationInfo();
  m_ati = new AmpToolsInterface(m_configurationInfo);

  for(unsigned int r=0; r<m_configurationInfo->reactionList().size();r++){
    string reacName = m_configurationInfo->reactionList().at(r)->reactionName();
    if(useMC)
      m_reader[reacName] = m_ati->accMCReader(reacName);
    else{
      m_useMC = false;
      m_gen[reacName] = new TGenPhaseSpace();
    }

    fOut->cd();
    m_tree[reacName] = new TTree(reacName.c_str(),reacName.c_str());
    m_tree[reacName]->Branch("weightAmp",&m_weight,"weightAmp/F");
    m_tree[reacName]->Branch("weight",&m_weightOrig,"weight/F");
    for(unsigned int p=0; p<m_configurationInfo->reactionList().at(r)->particleList().size(); p++){
      m_tree[reacName]->Branch(Form("PxP%i",p),&m_PxP[p],Form("PxP%i/F",p));
      m_tree[reacName]->Branch(Form("PyP%i",p),&m_PyP[p],Form("PyP%i/F",p));
      m_tree[reacName]->Branch(Form("PzP%i",p),&m_PzP[p],Form("PzP%i/F",p));
      m_tree[reacName]->Branch(Form("EnP%i",p),&m_EnP[p],Form("EnP%i/F",p));
    }
  }
}

void PlotInterface::setupFromFitResult(const char* fitresult){
  m_fitResults = new FitResults(fitresult);
  m_fitResults->configInfo()->write("xxxxxTEMPxxxxx.cfg");
  m_parser = new ConfigFileParser("xxxxxTEMPxxxxx.cfg");
  m_configurationInfo = m_parser->getConfigurationInfo();
  m_ati = new AmpToolsInterface(m_configurationInfo);

  for(unsigned int r=0; r<m_configurationInfo->reactionList().size();r++){
    string reacName = m_configurationInfo->reactionList().at(r)->reactionName();
    m_reader[reacName] = m_ati->accMCReader(reacName);

    fOut->cd();
    m_tree[reacName] = new TTree(reacName.c_str(),reacName.c_str());
    m_tree[reacName]->Branch("weightAmp",&m_weight,"weightAmp/F");
    m_tree[reacName]->Branch("weight",&m_weightOrig,"weight/F");
    for(unsigned int p=0; p<m_configurationInfo->reactionList().at(r)->particleList().size(); p++){
      m_tree[reacName]->Branch(Form("PxP%i",p),&m_PxP[p],Form("PxP%i/F",p));
      m_tree[reacName]->Branch(Form("PyP%i",p),&m_PyP[p],Form("PyP%i/F",p));
      m_tree[reacName]->Branch(Form("PzP%i",p),&m_PzP[p],Form("PzP%i/F",p));
      m_tree[reacName]->Branch(Form("EnP%i",p),&m_EnP[p],Form("EnP%i/F",p));
    }
  }
}

void PlotInterface::plot(string reaction){
  if(m_useMC){
    m_reader[reaction]->resetSource();
    m_ati->loadEvents(m_reader[reaction]);
    m_ati->processEvents(reaction);
    int numEvents = m_ati->numEvents();
    float totIntensity = 0.;

    for(int ev=0;ev<numEvents;ev++){

      Kinematics *kin = m_ati->kinematics(ev);
      
      if(kin == NULL){
        break;
      }

      m_weight = m_ati->intensity(ev);
      const vector<TLorentzVector> vecs = kin->particleList();
      m_weightOrig = kin->weight();
      for(unsigned int p=0;p<vecs.size();p++){
        m_PxP[p] = vecs.at(p).Px();
        m_PyP[p] = vecs.at(p).Py();
        m_PzP[p] = vecs.at(p).Pz();
        m_EnP[p] = vecs.at(p).E();
      }


      for(auto &p : m_amps_to_draw){
        pair<float,complex<float>> _vals = calcAmp(ev,reaction,p.second);
        m_amp_weights_INT[p.second] = _vals.first/numEvents;
        m_amp_weights_re[p.second] = _vals.second.real()/sqrt(numEvents);
        m_amp_weights_im[p.second] = _vals.second.imag()/sqrt(numEvents);
      }

      totIntensity += m_weight;
      
      m_tree[reaction]->Fill();
      delete kin;
    }
  }
  else{

  }
  fOut->cd();
  m_tree[reaction]->Write();
}


// this is from AmpToolsInterface::alternateIntensity
pair<float,complex<float>> PlotInterface::calcAmp(int iEvent, string reaction, string name){
  float intensity = 0.;
  complex<float> amp(0.,0.);
  vector<CoherentSumInfo*> sums = m_configurationInfo->coherentSumList(reaction);
  for (unsigned int iSum = 0; iSum < sums.size(); iSum++){
    complex<double> runningAmplitude(0.0,0.0);
    vector<AmplitudeInfo*> amps = m_configurationInfo->amplitudeList(reaction,sums[iSum]->sumName());
    for(unsigned int a=0; a < amps.size(); a++){
      if(amps.at(a)->ampName().find(name.c_str())!=string::npos){
        complex<double> P = m_ati->scaledProductionAmplitude(amps[a]->fullName());
        complex<double> D = m_ati->decayAmplitude(iEvent,amps[a]->fullName(),0);
        runningAmplitude += P*D;
      }
    }
    amp = runningAmplitude;
    intensity += norm(runningAmplitude);
  }

  return make_pair(intensity,amp);
}

void PlotInterface::writeOutput(){
  fOut->Write();
  fOut->Close();
}

void PlotInterface::clear(){

}



