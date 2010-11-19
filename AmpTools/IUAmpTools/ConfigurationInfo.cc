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
#include <cstdlib>

#include "IUAmpTools/ConfigurationInfo.h"


vector<ReactionInfo*>
ConfigurationInfo::reactionList  (const string& reactionName) const {

  vector<ReactionInfo*> reactions;
  for (unsigned int i = 0; i < m_reactions.size(); i++){
    if ((reactionName == "") || (m_reactions[i]->reactionName() == reactionName)){
      reactions.push_back(m_reactions[i]);
    }
  }

  return reactions;
}


vector<CoherentSumInfo*>
ConfigurationInfo::coherentSumList  (const string& reactionName,
                                     const string& sumName) const {
  vector<CoherentSumInfo*> sums;
  for (unsigned int i = 0; i < m_sums.size(); i++){
    if (((reactionName == "") || (m_sums[i]->reactionName() == reactionName)) &&
        ((sumName      == "") || (m_sums[i]->sumName()      == sumName))){
      sums.push_back(m_sums[i]);
    }
  }
  return sums;
}


vector<AmplitudeInfo*>
ConfigurationInfo::amplitudeList  (const string& reactionName,
                                   const string& sumName,
                                   const string& ampName) const {
  vector<AmplitudeInfo*> amplitudes;
  for (unsigned int i = 0; i < m_amplitudes.size(); i++){
    if (((reactionName == "") || (m_amplitudes[i]->reactionName() == reactionName)) &&
        ((sumName      == "") || (m_amplitudes[i]->sumName()      == sumName)) &&
        ((ampName      == "") || (m_amplitudes[i]->ampName()      == ampName))){
      amplitudes.push_back(m_amplitudes[i]);
    }
  }
  return amplitudes;
}

vector<ParameterInfo*>
ConfigurationInfo::parameterList  (const string& reactionName,
                                   const string& sumName,
                                   const string& ampName,
                                   const string& parName) const {
  if ((reactionName == "") &&
      (sumName      == "") &&
      (ampName      == "") &&
      (parName      == "")) return m_parameters;
  vector<ParameterInfo*> parameters;
  for (unsigned int i = 0; i < m_amplitudes.size(); i++){
    if (((reactionName == "") || (m_amplitudes[i]->reactionName() == reactionName)) &&
        ((sumName      == "") || (m_amplitudes[i]->sumName()      == sumName)) &&
        ((ampName      == "") || (m_amplitudes[i]->ampName()      == ampName))){
      vector<ParameterInfo*> ampParameters = m_amplitudes[i]->parameters();
      for (unsigned int j = 0; j < ampParameters.size(); j++){
        bool foundParameter = false;
        for (unsigned int k = 0; k < parameters.size(); k++){
          if (parameters[k]->parName() == ampParameters[j]->parName()) foundParameter = true;
        }
        if (!foundParameter){
          if ((parName == "") || (ampParameters[j]->parName() == parName)){
            parameters.push_back(ampParameters[j]);
          }
        }
      }
    }
  }
  return parameters;
}



ReactionInfo*
ConfigurationInfo::reaction  (const string& reactionName) const {
  for (unsigned int i = 0; i < m_reactions.size(); i++){
    if (m_reactions[i]->reactionName() == reactionName){
      return m_reactions[i];
    }
  }
  ReactionInfo* reaction = NULL;
  return reaction;
}


CoherentSumInfo*
ConfigurationInfo::coherentSum  (const string& reactionName,
                                 const string& sumName) const {
  for (unsigned int i = 0; i < m_sums.size(); i++){
    if ((m_sums[i]->reactionName() == reactionName) && 
        (m_sums[i]->sumName()      == sumName)){
      return m_sums[i];
    }
  }
  CoherentSumInfo* sum = NULL;
  return sum;
}


AmplitudeInfo*
ConfigurationInfo::amplitude  (const string& reactionName,
                               const string& sumName,
                               const string& ampName) const {
  for (unsigned int i = 0; i < m_amplitudes.size(); i++){
    if ((m_amplitudes[i]->reactionName() == reactionName) && 
        (m_amplitudes[i]->sumName()      == sumName) &&
        (m_amplitudes[i]->ampName()      == ampName)){
      return m_amplitudes[i];
    }
  }
  AmplitudeInfo* amplitude = NULL;
  return amplitude;
}


ParameterInfo*
ConfigurationInfo::parameter  (const string& parName) const {
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    if (m_parameters[i]->parName() == parName){
      return m_parameters[i];
    }
  }
  ParameterInfo* parameter = NULL;
  return parameter;
}



ReactionInfo*
ConfigurationInfo::createReaction  (const string&         reactionName, 
                                    const vector<string>& particleList){
  ReactionInfo* rctn = reaction(reactionName);
  if (rctn == NULL){
    rctn = new ReactionInfo(reactionName,particleList);
    m_reactions.push_back(rctn);
  }
  else{
    rctn->setParticleList(particleList);
    rctn->clear();
  }
  return rctn;
}


CoherentSumInfo*
ConfigurationInfo::createCoherentSum  (const string& reactionName, 
                                       const string& sumName){
  CoherentSumInfo* sum = coherentSum(reactionName,sumName);
  if (sum == NULL){
    ReactionInfo* rctn = reaction(reactionName);
    if (rctn == NULL){
      cout << "ConfigurationInfo:  Trying to create coherentSum for unknown reaction" << endl;
      cout << "\tcreateCoherentSum(" << reactionName << "," << sumName << ")" << endl;
      exit(0);
    }
    sum = new CoherentSumInfo(reactionName,sumName);
    m_sums.push_back(sum);
  }
  else{
    sum->clear();
  }
  return sum;
}



AmplitudeInfo*
ConfigurationInfo::createAmplitude  (const string& reactionName, 
                                     const string& sumName, 
                                     const string& ampName){
  AmplitudeInfo* amp = amplitude(reactionName,sumName,ampName);
  if (amp == NULL){
    ReactionInfo* rctn = reaction(reactionName);
    if (rctn == NULL){
      cout << "ConfigurationInfo:  Trying to create amplitude for unknown reaction" << endl;
      cout << "\tcreateAmplitude(" << reactionName << "," << sumName << "," << ampName << ")" << endl;
      exit(0);
    }
    CoherentSumInfo* sum = coherentSum(reactionName,sumName);
    if (sum == NULL){
      cout << "ConfigurationInfo:  Trying to create amplitude for unknown coherent sum" << endl;
      cout << "\tcreateAmplitude(" << reactionName << "," << sumName << "," << ampName << ")" << endl;
      exit(0);
    }
    amp = new AmplitudeInfo(reactionName,sumName,ampName);
    m_amplitudes.push_back(amp);
  }
  else{
    amp->clear();
  }
  return amp;
}



ParameterInfo*
ConfigurationInfo::createParameter  (const string&  parName,
                                     double         value){
  ParameterInfo* par = parameter(parName);
  if (par == NULL){
    par = new ParameterInfo(parName,value);
    m_parameters.push_back(par);
  }
  else{
    par->setValue(value);
    par->clear();
  }
  return par;
}



void
ConfigurationInfo::removeReaction (const string& reactionName){
  removeCoherentSum(reactionName);
  vector<ReactionInfo*> removeList = reactionList(reactionName);
  unsigned int removeListSize = removeList.size();
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_reactions.size(); j++){
      if (removeList[i]->reactionName() == m_reactions[j]->reactionName()){
        m_reactions[j]->clear();
        delete m_reactions[j];
        m_reactions.erase(m_reactions.begin()+j,m_reactions.begin()+j+1);
        j = m_reactions.size();
      }
    }
  }
}

void
ConfigurationInfo::removeCoherentSum (const string& reactionName,
                                      const string& sumName){
  removeAmplitude(reactionName,sumName);
  vector<CoherentSumInfo*> removeList = coherentSumList(reactionName, sumName);
  unsigned int removeListSize = removeList.size();
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_sums.size(); j++){
      if ((removeList[i]->reactionName() == m_sums[j]->reactionName()) &&
          (removeList[i]->sumName()      == m_sums[j]->sumName())){
        m_sums[j]->clear();
        delete m_sums[j];
        m_sums.erase(m_sums.begin()+j,m_sums.begin()+j+1);
        j = m_sums.size();
      }
    }
  }
}

void
ConfigurationInfo::removeAmplitude   (const string& reactionName,
                                      const string& sumName,
                                      const string& ampName){
  vector<AmplitudeInfo*> removeList = amplitudeList(reactionName, sumName, ampName);
  unsigned int removeListSize = removeList.size();
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_amplitudes.size(); j++){
      m_amplitudes[j]->removeConstraint(removeList[i]);
    }
  }
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_amplitudes.size(); j++){
      if ((removeList[i]->reactionName() == m_amplitudes[j]->reactionName()) &&
          (removeList[i]->sumName()      == m_amplitudes[j]->sumName()) &&
          (removeList[i]->ampName()      == m_amplitudes[j]->ampName())){
        m_amplitudes[j]->clear();
        delete m_amplitudes[j];
        m_amplitudes.erase(m_amplitudes.begin()+j,m_amplitudes.begin()+j+1);
        j = m_amplitudes.size();
      }
    }
  }
}


void
ConfigurationInfo::removeParameter   (const string& parName){
  vector<ParameterInfo*> removeList;
  if (parName == ""){
    removeList = m_parameters;
  }
  else{
    ParameterInfo* parinfo = parameter(parName);
    if (parinfo) removeList.push_back(parinfo);
  }
  unsigned int removeListSize = removeList.size();
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_amplitudes.size(); j++){
      m_amplitudes[j]->removeParameter(removeList[i]);
    }
  }
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_parameters.size(); j++){
      if (removeList[i]->parName() == m_parameters[j]->parName()){
        m_parameters[j]->clear();
        delete m_parameters[j];
        m_parameters.erase(m_parameters.begin()+j,m_parameters.begin()+j+1);
        j = m_parameters.size();
      }
    }
  }
}

map< string, vector< string > >
ConfigurationInfo::constraintMap() const {
  
  map< string, vector< string > > cMap;
  
  for( vector< AmplitudeInfo* >::const_iterator aItr = m_amplitudes.begin();
      aItr != m_amplitudes.end();
      ++aItr ){
    
    cMap[(**aItr).fullName()] = vector< string >( 0 );
    
    vector< AmplitudeInfo* > constraints = (**aItr).constraints();
    
    for( vector< AmplitudeInfo* >::const_iterator bItr = constraints.begin();
        bItr != constraints.end();
        ++bItr ){
      
      cMap[(**aItr).fullName()].push_back( (**bItr).fullName() );
    }
  }
  
  return cMap;
}


ConfigurationInfo::~ConfigurationInfo(){
  removeReaction();
  removeParameter();
}


void
ConfigurationInfo::display(string fileName, bool append){

  ofstream outfile;
  streambuf* cout_sbuf = cout.rdbuf();
  if (fileName != ""){
    if (append){
      outfile.open(fileName.c_str(), ios::out | ios::app);
    }
    else{
      outfile.open(fileName.c_str());
    }
    cout.rdbuf(outfile.rdbuf());
  }

  cout << endl;
  cout << "## CONFIGURATION INFO DISPLAY ##" << endl;
  cout << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

  for (unsigned int i = 0; i < m_reactions.size(); i++){
    m_reactions[i]->display(fileName,true);
    vector<CoherentSumInfo*> sums = coherentSumList(m_reactions[i]->reactionName());
    for (unsigned int j = 0; j < sums.size(); j++){
      sums[j]->display(fileName,true);
      vector<AmplitudeInfo*> amps = amplitudeList(m_reactions[i]->reactionName(),sums[j]->sumName());
      for (unsigned int k = 0; k < amps.size(); k++){
        amps[k]->display(fileName,true);
      }
    }
  }
  for (unsigned int l = 0; l < m_parameters.size(); l++){
    m_parameters[l]->display(fileName,true);
  }
}



void
ReactionInfo::display(string fileName, bool append){

  ofstream outfile;
  streambuf* cout_sbuf = cout.rdbuf();
  if (fileName != ""){
    if (append){
      outfile.open(fileName.c_str(), ios::out | ios::app);
    }
    else{
      outfile.open(fileName.c_str());
    }
    cout.rdbuf(outfile.rdbuf());
  }

  cout << "############################################" << endl;
  cout << "#############   REACTION INFO  #############" << endl;
  cout << "############################################" << endl;
  cout << "      REACTION NAME:  " << m_reactionName << endl;
  cout << "      PARTICLE LIST:  " << m_particleList.size() << endl;
  for (unsigned int i = 0; i < m_particleList.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_particleList[i] << endl;
  }
  cout << "         DATA FILES:  " << m_dataFiles.size() << endl;
  for (unsigned int i = 0; i < m_dataFiles.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_dataFiles[i] << endl;
  }
  cout << " GENERATED MC FILES:  " << m_genMCFiles.size() << endl;
  for (unsigned int i = 0; i < m_genMCFiles.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_genMCFiles[i] << endl;
  }
  cout << "  ACCEPTED MC FILES:  " << m_accMCFiles.size() << endl;
  for (unsigned int i = 0; i < m_accMCFiles.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_accMCFiles[i] << endl;
  }
  cout << "  NORMALIZATION INTEGRAL FILE: " << endl;
  if (m_normIntFile != "")  cout << "\t\t    " << m_normIntFile << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

}


void
CoherentSumInfo::display(string fileName, bool append){

  ofstream outfile;
  streambuf* cout_sbuf = cout.rdbuf();
  if (fileName != ""){
    if (append){
      outfile.open(fileName.c_str(), ios::out | ios::app);
    }
    else{
      outfile.open(fileName.c_str());
    }
    cout.rdbuf(outfile.rdbuf());
  }

  cout << "********************************************" << endl;
  cout << "***********  COHERENT SUM INFO  ************" << endl;
  cout << "********************************************" << endl;
  cout << "      REACTION NAME:  " << m_reactionName << endl;
  cout << "  COHERENT SUM NAME:  " << m_sumName << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

}


void
AmplitudeInfo::display(string fileName, bool append){

  ofstream outfile;
  streambuf* cout_sbuf = cout.rdbuf();
  if (fileName != ""){
    if (append){
      outfile.open(fileName.c_str(), ios::out | ios::app);
    }
    else{
      outfile.open(fileName.c_str());
    }
    cout.rdbuf(outfile.rdbuf());
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++  AMPLITUDE INFO  +++++++++++++" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "      REACTION NAME:  " << m_reactionName << endl;
  cout << "  COHERENT SUM NAME:  " << m_sumName << endl;
  cout << "     AMPLITUDE NAME:  " << m_ampName << endl;
  cout << "            FACTORS:  " << m_factors.size() << endl;
  for (unsigned int i = 0; i < m_factors.size(); i++){
    vector<string> factor = m_factors[i];
    cout << "\t\t" << i+1 << ".  ";
    for (unsigned int j = 0; j < factor.size(); j++){
      cout << " " << factor[j];
    }
    cout << endl;
  }
  cout << "       PERMUTATIONS:  " << m_permutations.size() << endl;
  for (unsigned int i = 0; i < m_permutations.size(); i++){
    vector<int> permutation = m_permutations[i];
    cout << "\t\t" << i+1 << ".  ";
    for (unsigned int j = 0; j < permutation.size(); j++){
      cout << " " << permutation[j];
    }
    cout << endl;
  }
  cout << "        CONSTRAINTS:  " << m_constraints.size() << endl;
  for (unsigned int i = 0; i < m_constraints.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_constraints[i]->reactionName() 
                          << "  "  << m_constraints[i]->sumName()
                          << "  "  << m_constraints[i]->ampName() << endl;
  }
  cout << "         PARAMETERS:  " << m_parameters.size() << endl;
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    cout << "\t\t" << i+1 << ".  " << m_parameters[i]->parName() << endl;
  }
  cout << "      INITIAL VALUE:  " << m_value << endl;
  cout << "               REAL?  " << m_real << endl;
  cout << "              FIXED?  " << m_fixed << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

}


void
ParameterInfo::display(string fileName, bool append){

  ofstream outfile;
  streambuf* cout_sbuf = cout.rdbuf();
  if (fileName != ""){
    if (append){
      outfile.open(fileName.c_str(), ios::out | ios::app);
    }
    else{
      outfile.open(fileName.c_str());
    }
    cout.rdbuf(outfile.rdbuf());
  }

  cout << "--------------------------------------------" << endl;
  cout << "-------------  PARAMETER INFO  -------------" << endl;
  cout << "--------------------------------------------" << endl;
  cout << "     PARAMETER NAME:  " << m_parName << endl;
  cout << "      INITIAL VALUE:  " << m_value << endl;
  cout << "              FIXED?  " << m_fixed << endl;
  cout << "            BOUNDED?  " << m_bounded << endl;
  cout << "        LOWER BOUND:  " << m_lowerBound << endl;
  cout << "        UPPER BOUND:  " << m_upperBound << endl;
  cout << "   GAUSSIAN BOUNDED?  " << m_gaussianBounded << endl;
  cout << "      CENTRAL VALUE:  " << m_centralValue << endl;
  cout << "     GAUSSIAN ERROR:  " << m_gaussianError << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

}


/*
void
ConfigurationInfo::write(string fileName){
  ofstream ff(fileName.c_str());
  for (unsigned int i = 0; i < m_reactions.size(); i++){
    ff << "REACTION" << endl;
    ff << m_reactions[i]->reactionName() << endl;
    vector<string> pList = m_reactions[i]->particleList();
    ff << pList.size() << endl;
    for (unsigned int j = 0; j < pList.size(); j++){ ff << pList[j] << endl; }
    vector<string> dataFiles = m_reactions[i]->dataFiles();
    ff << dataFiles.size() << endl;
    for (unsigned int j = 0; j < dataFiles.size(); j++){ ff << dataFiles[j] << endl; }
    vector<string> genMCFiles = m_reactions[i]->genMCFiles();
    ff << genMCFiles.size() << endl;
    for (unsigned int j = 0; j < genMCFiles.size(); j++){ ff << genMCFiles[j] << endl; }
    vector<string> accMCFiles = m_reactions[i]->accMCFiles();
    ff << accMCFiles.size() << endl;
    for (unsigned int j = 0; j < accMCFiles.size(); j++){ ff << accMCFiles[j] << endl; }
    ff << m_reactions[i]->normIntFile() << endl;
    vector<AmplitudeInfo*> amps = amplitudeList(m_reactions[i]->reactionName());
    for (unsigned int j = 0; j < amps.size(); j++){
      ff << "AMPLITUDE" << endl;
      ff << amps[j]->reactionName() << endl;
      ff << amps[j]->sumName() << endl;
      ff << amps[j]->ampName() << endl;
      vector< vector<string> > factors = amps[j]->factors();
      ff << factors.size() << endl;
      for (unsigned int k = 0; k < factors.size(); k++){
        vector<string> factor = factors[k];
        ff << factor.size() << endl;
        for (unsigned int l = 0; l < factor.size(); l++){
          ff << factor[l] << endl;
        }
      }
      vector< vector<int> > permutations = amps[j]->permutations();
      ff << permutations.size() << endl;
      for (unsigned int k = 0; k < permutations.size(); k++){
        vector<int> permutation = permutations[k];
        ff << permutation.size() << endl;
        for (unsigned int l = 0; l < permutation.size(); l++){
          ff << permutation[l] << endl;
        }
      }
      vector< AmplitudeInfo* > constraints = amps[j]->constraints();
      ff << constraints.size() << endl;
      for (unsigned int k = 0; k < constraints.size(); k++){
        ff << constraints[k]->reactionName() << endl;
        ff << constraints[k]->sumName() << endl;
        ff << constraints[k]->ampName() << endl;
      }
    }
  }
  ff.close();
  for (unsigned int i = 0; i < m_reactions.size(); i++){
    m_reactions[i]->display();
    cout << "  SUM NAMES: " << endl;
    vector<string> sums = sumNameList(m_reactions[i]->reactionName());
    for (unsigned int j = 0; j < sums.size(); j++){
      cout << "\t" << j+1 << ".  " << sums[j] << endl;
    }
    vector<AmplitudeInfo*> amps = amplitudeList(m_reactions[i]->reactionName());
    for (unsigned int j = 0; j < amps.size(); j++){
      amps[j]->display();
    }
  }
  for (unsigned int k = 0; k < m_parameters.size(); k++){
    m_parameters[k]->display();
  }
}
*/

/*
void
ConfigurationInfo::read(string fileName){
  ifstream ff("test.txt");
  string type;
  ff >> type;
  if (type == "PARAMETER"){
    ff >> m_parName;
    ff >> m_value;
    ff >> m_fixed;
    ff >> m_bounded;
    ff >> m_lowerBound;
    ff >> m_upperBound;
    ff >> m_gaussianBounded;
    ff >> m_centralValue;
    ff >> m_gaussianError;
  }
  ff.close();
}
*/



void 
AmplitudeInfo::addConstraint (AmplitudeInfo* constraint){
    // don't constrain an amplitude to itself
  if ((this->reactionName() == constraint->reactionName()) &&
      (this->sumName()      == constraint->sumName()) &&
      (this->ampName()      == constraint->ampName())) return;
    // add "constraint" as a constraint
  if (!hasConstraint(constraint)) m_constraints.push_back(constraint);
    // also add all of "constraint"'s constraints as constraints
  vector<AmplitudeInfo*> constraints = constraint->constraints();
  for (unsigned int i = 0; i < constraints.size(); i++){
    if (!hasConstraint(constraints[i])) addConstraint(constraints[i]);
  }
    // also reciprocate
  if (!(constraint->hasConstraint(this))) constraint->addConstraint(this);
}


bool 
AmplitudeInfo::hasConstraint(AmplitudeInfo* constraint) const{
  bool foundConstraint = false;
  for (unsigned int i = 0; i < m_constraints.size(); i++){
    if ((m_constraints[i]->reactionName() == constraint->reactionName()) &&
        (m_constraints[i]->sumName()      == constraint->sumName()) &&
        (m_constraints[i]->ampName()      == constraint->ampName())) foundConstraint = true;
  }
  return foundConstraint;
}


void
AmplitudeInfo::removeConstraint(AmplitudeInfo* constraint){
    // remove "constraint" as a constraint
  if (hasConstraint(constraint)){
    for (unsigned int i = 0; i < m_constraints.size(); i++){
      if ((m_constraints[i]->reactionName() == constraint->reactionName()) &&
          (m_constraints[i]->sumName()      == constraint->sumName()) &&
          (m_constraints[i]->ampName()      == constraint->ampName())){
        m_constraints.erase(m_constraints.begin()+i,m_constraints.begin()+i+1);
        i = m_constraints.size();
      }
    }
  }
    // remove "constraint"'s constraints as constraints 
  vector<AmplitudeInfo*> constraints = constraint->constraints();
  for (unsigned int i = 0; i < constraints.size(); i++){
    if (hasConstraint(constraints[i])) removeConstraint(constraints[i]);
  }
    // reciprocate
  if (constraint->hasConstraint(this)) constraint->removeConstraint(this);
}


void 
AmplitudeInfo::addParameter (ParameterInfo* parameter){
  bool foundParameter = false;
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    if (m_parameters[i]->parName() == parameter->parName()) foundParameter = true;
  }
  if (!foundParameter) m_parameters.push_back(parameter);
}

void
AmplitudeInfo::removeParameter (ParameterInfo* parameter){
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    if (m_parameters[i]->parName() == parameter->parName()){
      m_parameters.erase(m_parameters.begin()+i,m_parameters.begin()+i+1);
      i = m_parameters.size();
    }
  }
}


void
ReactionInfo::clear(){
  m_dataFiles.clear();
  m_genMCFiles.clear();
  m_accMCFiles.clear();
  m_normIntFile = "";
}

void
AmplitudeInfo::clear(){
  m_factors.clear();
  m_permutations.clear();
  m_constraints.clear();
  m_value = complex<double>(0.0,0.0);
  m_real = false;
  m_fixed = false;
  m_parameters.clear();
}

void
ParameterInfo::clear(){
  m_fixed = false;
  m_bounded = false;
  m_lowerBound = 0.0;
  m_upperBound = 0.0;
  m_gaussianBounded = false;
  m_centralValue = 0.0;
  m_gaussianError = 0.0;
}





