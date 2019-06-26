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


ConfigurationInfo::~ConfigurationInfo(){
  removeParameter();
  removePDF();
  removeAmplitude();
  removeCoherentSum();
  removeReaction();
}


vector< string >
ConfigurationInfo::userKeywords() const{

  vector<string> keywords;
  for (map<string, vector< vector<string> > >::const_iterator 
       mapItr = m_userKeywordMap.begin();
       mapItr != m_userKeywordMap.end(); mapItr++){
    keywords.push_back(mapItr->first);
  }
  return keywords;
}


vector< vector<string> >
ConfigurationInfo::userKeywordArguments(const string& userKeyword) const{

  if (m_userKeywordMap.find(userKeyword) != m_userKeywordMap.end())
    return m_userKeywordMap.find(userKeyword)->second;
  return vector< vector<string> >();
}


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

vector<PDFInfo*>
ConfigurationInfo::pdfList  (const string& reactionName,
                             const string& pdfName) const {
  vector<PDFInfo*> pdfs;
  for (unsigned int i = 0; i < m_pdfs.size(); i++){
    if (((reactionName == "") || (m_pdfs[i]->reactionName() == reactionName)) &&
        ((pdfName      == "") || (m_pdfs[i]->pdfName()      == pdfName))){
      pdfs.push_back(m_pdfs[i]);
    }
  }
  return pdfs;
}

vector<TermInfo*>
ConfigurationInfo::termList  (const string& reactionName,
                              const string& sumName,
                              const string& termName) const {
  vector<TermInfo*> terms;
  vector<AmplitudeInfo*> amplitudes = amplitudeList(reactionName,sumName,termName);
  vector<PDFInfo*> pdfs = pdfList(reactionName,termName);
  for (unsigned int i = 0; i < amplitudes.size(); i++){ terms.push_back(amplitudes[i]); }
  for (unsigned int i = 0; i < pdfs.size();       i++){ terms.push_back(pdfs[i]); }
  return terms;
}


vector<ParameterInfo*>
ConfigurationInfo::parameterList  (const string& reactionName,
                                   const string& sumName,
                                   const string& termName,
                                   const string& parName) const {
  if ((reactionName == "") &&
      (sumName      == "") &&
      (termName     == "") &&
      (parName      == "")) return m_parameters;
  vector<ParameterInfo*> parameters;
  for (unsigned int i = 0; i < m_amplitudes.size(); i++){
    if (((reactionName == "") || (m_amplitudes[i]->reactionName() == reactionName)) &&
        ((sumName      == "") || (m_amplitudes[i]->sumName()      == sumName)) &&
        ((termName     == "") || (m_amplitudes[i]->ampName()      == termName))){
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
  for (unsigned int i = 0; i < m_pdfs.size(); i++){
    if (((reactionName == "") || (m_pdfs[i]->reactionName() == reactionName)) &&
        ((termName     == "") || (m_pdfs[i]->pdfName()      == termName))){
      vector<ParameterInfo*> pdfParameters = m_pdfs[i]->parameters();
      for (unsigned int j = 0; j < pdfParameters.size(); j++){
        bool foundParameter = false;
        for (unsigned int k = 0; k < parameters.size(); k++){
          if (parameters[k]->parName() == pdfParameters[j]->parName()) foundParameter = true;
        }
        if (!foundParameter){
          if ((parName == "") || (pdfParameters[j]->parName() == parName)){
            parameters.push_back(pdfParameters[j]);
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


AmplitudeInfo*
ConfigurationInfo::amplitude  (const string& fullName) const {
  for (unsigned int i = 0; i < m_amplitudes.size(); i++){
    if (m_amplitudes[i]->fullName() == fullName){
      return m_amplitudes[i];
    }
  }
  AmplitudeInfo* amplitude = NULL;
  return amplitude;
}


PDFInfo*
ConfigurationInfo::pdf  (const string& reactionName,
                         const string& pdfName) const {
  for (unsigned int i = 0; i < m_pdfs.size(); i++){
    if ((m_pdfs[i]->reactionName() == reactionName) && 
        (m_pdfs[i]->pdfName()      == pdfName)){
      return m_pdfs[i];
    }
  }
  PDFInfo* lpdf = NULL;
  return lpdf;
}

PDFInfo*
ConfigurationInfo::pdf  (const string& fullName) const {
  for (unsigned int i = 0; i < m_pdfs.size(); i++){
    if (m_pdfs[i]->fullName() == fullName){
      return m_pdfs[i];
    }
  }
  PDFInfo* lpdf = NULL;
  return lpdf;
}

TermInfo*
ConfigurationInfo::term (const string& fullName) const {
  AmplitudeInfo* amp = amplitude(fullName);  if (amp) return amp;
  PDFInfo* lpdf = pdf(fullName);  if (lpdf) return lpdf;
  TermInfo* term = NULL;
  return term;
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



void
ConfigurationInfo::addUserKeyword(const string& userKeyword, const vector<string>& arguments){

  vector< vector<string> > argList = m_userKeywordMap[userKeyword];
  argList.push_back(arguments);
  m_userKeywordMap[userKeyword] = argList;
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
  if (pdfList(reactionName).size() != 0){
    cout << "ConfigurationInfo:  problem with reaction named " << reactionName << ":" << endl;
    cout << "   mixing amplitudes and pdfs is not currently supported" << endl;
    exit(0);
  }
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


PDFInfo*
ConfigurationInfo::createPDF  (const string& reactionName, 
                               const string& pdfName){
  if (amplitudeList(reactionName).size() != 0){
    cout << "ConfigurationInfo:  problem with reaction named " << reactionName << ":" << endl;
    cout << "   mixing amplitudes and pdfs is not currently supported" << endl;
    exit(0);
  }
  PDFInfo* addpdf = pdf(reactionName,pdfName);
  if (addpdf == NULL){
    ReactionInfo* rctn = reaction(reactionName);
    if (rctn == NULL){
      cout << "ConfigurationInfo:  Trying to create pdf for unknown reaction" << endl;
      cout << "\tcreatePDF(" << reactionName << "," << pdfName << ")" << endl;
      exit(0);
    }
    addpdf = new PDFInfo(reactionName,pdfName);
    m_pdfs.push_back(addpdf);
  }
  else{
    addpdf->clear();
  }
  return addpdf;
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
ConfigurationInfo::removeUserKeyword(const string& userKeyword){

  if (userKeyword == "") m_userKeywordMap.clear();
  m_userKeywordMap.erase(userKeyword);
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
ConfigurationInfo::removePDF   (const string& reactionName,
                               const string& pdfName){
  vector<PDFInfo*> removeList = pdfList(reactionName, pdfName);
  unsigned int removeListSize = removeList.size();
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_pdfs.size(); j++){
      m_pdfs[j]->removeConstraint(removeList[i]);
    }
  }
  for (unsigned int i = 0; i < removeListSize; i++){
    for (unsigned int j = 0; j < m_pdfs.size(); j++){
      if ((removeList[i]->reactionName() == m_pdfs[j]->reactionName()) &&
          (removeList[i]->pdfName()      == m_pdfs[j]->pdfName())){
        m_pdfs[j]->clear();
        delete m_pdfs[j];
        m_pdfs.erase(m_pdfs.begin()+j,m_pdfs.begin()+j+1);
        j = m_pdfs.size();
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
    for (unsigned int j = 0; j < m_pdfs.size(); j++){
      m_pdfs[j]->removeParameter(removeList[i]);
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
    vector< TermInfo* > constraints = (**aItr).constraints();
    for( vector< TermInfo* >::const_iterator bItr = constraints.begin();
        bItr != constraints.end();
        ++bItr ){
      cMap[(**aItr).fullName()].push_back( (**bItr).fullName() );
    }
  }

  for( vector< PDFInfo* >::const_iterator aItr = m_pdfs.begin();
      aItr != m_pdfs.end();
      ++aItr ){
    cMap[(**aItr).fullName()] = vector< string >( 0 );
    vector< TermInfo* > constraints = (**aItr).constraints();
    for( vector< TermInfo* >::const_iterator bItr = constraints.begin();
        bItr != constraints.end();
        ++bItr ){
      cMap[(**aItr).fullName()].push_back( (**bItr).fullName() );
    }
  }
  
  return cMap;
}


void
ConfigurationInfo::write( const string& fileName ) const {
  
  ofstream ff(fileName.c_str());
  write( ff );
  ff.close();
}

ostream&
ConfigurationInfo::write( ostream& ff ) const {
  
  vector<ReactionInfo*>    Rs   = reactionList();
  vector<CoherentSumInfo*> Ss   = coherentSumList();
  vector<AmplitudeInfo*>   As   = amplitudeList();
  vector<PDFInfo*>         PDFs = pdfList();
  vector<ParameterInfo*>   Ps   = parameterList();
  
  ff << "### FIT CONFIGURATION ###" << endl;
  
  
  // fit
  
  ff << "fit " << fitName() << endl;
  
  
  // user keywords
  
  vector<string> keywords = userKeywords();
  for (unsigned int i = 0; i < keywords.size(); i++){
    vector< vector<string> > argslist = userKeywordArguments(keywords[i]);
    int min = 0; int max = 0;
    for (unsigned int j = 0; j < argslist.size(); j++){
      if (argslist[j].size() < min) min = argslist[j].size();
      if (argslist[j].size() > max) max = argslist[j].size();
    }
    ff << "keyword " << keywords[i] << " " << min << " " << max << endl;
    for (unsigned int j = 0; j < argslist.size(); j++){
      vector<string> args = argslist[j];
      ff << keywords[i] << " ";
      for (unsigned int k = 0; k < args.size(); k++){
        ff << args[k] << " ";
      }
      ff << endl;
    }
  }
  
  
  // reaction
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    ff << "reaction " << R->reactionName();
    vector<string> particles = R->particleList();
    for (unsigned int j = 0; j < particles.size(); j++){ ff << " " << particles[j];}
    ff << endl;}
  
  
  // sum
  
  for (unsigned int i = 0; i < Ss.size(); i++){
    CoherentSumInfo* S = Ss[i];
    ff << "sum " << S->fullName() << endl;}
  
  
  // amplitude
  
  for (unsigned int i = 0; i < As.size(); i++){
    AmplitudeInfo* A = As[i];
    vector< vector<string> > Fs = A->factors();
    for (unsigned int j = 0; j < Fs.size(); j++){
      ff << "amplitude " << A->fullName();
      vector<string> args = Fs[j];
      for (unsigned int k = 0; k < args.size(); k++){
        ff << " " << args[k];}
      ff << endl;}}
  

  // pdf
  
  for (unsigned int i = 0; i < PDFs.size(); i++){
    PDFInfo* PDF = PDFs[i];
    vector< vector<string> > Fs = PDF->factors();
    for (unsigned int j = 0; j < Fs.size(); j++){
      ff << "pdf " << PDF->fullName();
      vector<string> args = Fs[j];
      for (unsigned int k = 0; k < args.size(); k++){
        ff << " " << args[k];}
      ff << endl;}}

  
  // parameter
  
  for (unsigned int i = 0; i < Ps.size(); i++){
    ParameterInfo* P = Ps[i];
    ff << "parameter " << P->parName() << " " << P->value();
    if (P->fixed()){ ff << " fixed"; }
    else if (P->bounded()){ ff << " bounded " << P->lowerBound() 
      << " " << P->upperBound(); }
    else if (P->gaussianBounded()) { ff << " gaussian " << P->centralValue()
      << " " << P->gaussianError(); }
    ff << endl;}
  
  
  // scale
  
  for (unsigned int i = 0; i < As.size(); i++){
    AmplitudeInfo* A = As[i];
    ff << "scale " << A->fullName() << " " << A->scale() << endl;}
  for (unsigned int i = 0; i < PDFs.size(); i++){
    PDFInfo* PDF = PDFs[i];
    ff << "pdfscale " << PDF->fullName() << " " << PDF->scale() << endl;}
  
  
  // constrain
  
  for (unsigned int i = 0; i < As.size(); i++){
    AmplitudeInfo* A = As[i];
    vector<TermInfo*> constr = A->constraints();
    for (unsigned int j = 0; j < constr.size(); j++){
      ff << "constrain " << A->fullName() << " " << constr[j]->fullName() << endl;}}
  for (unsigned int i = 0; i < PDFs.size(); i++){
    PDFInfo* PDF = PDFs[i];
    vector<TermInfo*> constr = PDF->constraints();
    for (unsigned int j = 0; j < constr.size(); j++){
      ff << "pdfconstrain " << PDF->fullName() << " " << constr[j]->fullName() << endl;}}
  
  
  // initialize
  
  for (unsigned int i = 0; i < As.size(); i++){
    AmplitudeInfo* A = As[i];
    ff << "initialize " << A->fullName() << " cartesian " 
    << A->value().real() << " "
    << A->value().imag();
    if      (!A->real() &&  A->fixed()) {ff << " fixed";}
    else if ( A->real() && !A->fixed()) {ff << " real";}
    else if ( A->real() &&  A->fixed()) {ff << " fixed real";}
    ff << endl;}
  for (unsigned int i = 0; i < PDFs.size(); i++){
    PDFInfo* PDF = PDFs[i];
    ff << "pdfinitialize " << PDF->fullName() << " "
    << PDF->value();
    if (PDF->fixed()) {ff << " fixed";}
    ff << endl;}
  
  
  // permute
  
  for (unsigned int i = 0; i < As.size(); i++){
    AmplitudeInfo* A = As[i];
    vector< vector<int> > perms = A->permutations();
    for (unsigned int j = 0; j < perms.size(); j++){
      ff << "permute " << A->fullName(); 
      vector<int> perm = perms[j];
      for (unsigned int k = 0; k < perm.size(); k++){
        ff << " " << perm[k];}
      ff << endl;}}
  
  
  // data
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    string cls = R->data().first;
    vector<string> args = R->data().second;
    if (cls != ""){
      ff << "data " << R->reactionName() << " " << cls;
      for (unsigned int j = 0; j < args.size(); j++){ ff << " " << args[j];}
      ff << endl;}}

  // bkgnd
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    string cls = R->bkgnd().first;
    vector<string> args = R->bkgnd().second;
    if (cls != ""){
      ff << "bkgnd " << R->reactionName() << " " << cls;
      for (unsigned int j = 0; j < args.size(); j++){ ff << " " << args[j];}
      ff << endl;}}
  
  // genmc
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    string cls = R->genMC().first;
    vector<string> args = R->genMC().second;
    if (cls != ""){
      ff << "genmc " << R->reactionName() << " " << cls;
      for (unsigned int j = 0; j < args.size(); j++){ ff << " " << args[j];}
      ff << endl;}}
  
  
  // accmc
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    string cls = R->accMC().first;
    vector<string> args = R->accMC().second;
    if (cls != ""){
      ff << "accmc " << R->reactionName() << " " << cls;
      for (unsigned int j = 0; j < args.size(); j++){ ff << " " << args[j];}
      ff << endl;}}
  
  
  // normintfile
  
  for (unsigned int i = 0; i < Rs.size(); i++){
    ReactionInfo* R = Rs[i];
    string ni = R->normIntFile();
    if ((ni != "") && !(R->normIntFileInput())){
      ff << "normintfile " << R->reactionName() << " " << ni << endl;}
    if ((ni != "") && (R->normIntFileInput())){
      ff << "normintfile " << R->reactionName() << " " << ni << " input" << endl;}}
  
  return ff;
}

void
ConfigurationInfo::display(string fileName, bool append) const {
  
  // this is a display function - it should not modify ConfiguraitonInfo
  // but all of the accessor functions are non-const to allow
  // flexibility in manipulating the configuration info so we will
  // cast away the constness and call the private non-const diplsay
  // function
  
  ConfigurationInfo* ci = const_cast< ConfigurationInfo* >( this );
  
  ci->display( fileName, append );
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
    vector<PDFInfo*> pdfs = pdfList(m_reactions[i]->reactionName());
    for (unsigned int k = 0; k < pdfs.size(); k++){
      pdfs[k]->display(fileName,true);
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
  cout << "        DATA READER:  " << m_data.first << endl;
  for (unsigned int i = 0; i < m_data.second.size(); i++){
    cout << "\t\t\t\t" << m_data.second[i] << endl;
  }
  if( m_bkgnd.first != "" ){
    cout << "        BACKGROUND READER:  " << m_bkgnd.first << endl;
    for (unsigned int i = 0; i < m_bkgnd.second.size(); i++){
      cout << "\t\t\t\t" << m_bkgnd.second[i] << endl;
    }
  }
  cout << "      ACC MC READER:  " << m_accMC.first << endl;
  for (unsigned int i = 0; i < m_accMC.second.size(); i++){
    cout << "\t\t\t\t" << m_accMC.second[i] << endl;
  }
  cout << "      GEN MC READER:  " << m_genMC.first << endl;
  for (unsigned int i = 0; i < m_genMC.second.size(); i++){
    cout << "\t\t\t\t" << m_genMC.second[i] << endl;
  }
  cout << "  NORMALIZATION INTEGRAL FILE: " << endl;
  if (m_normIntFile != "")  cout << "\t\t    " << m_normIntFile << endl;
  if (m_normIntFileInput)   cout << "\t\t       (use as input)" << endl;
  cout << "      GPU DEVICE NUMBER:  " << m_gpuDeviceNumber << endl;

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
TermInfo::addConstraint (TermInfo* constraint){
    // don't constrain an amplitude to itself
  if ((this->fullName() == constraint->fullName())) return;
    // add "constraint" as a constraint
  if (!hasConstraint(constraint)) m_constraints.push_back(constraint);
    // also add all of "constraint"'s constraints as constraints
  vector<TermInfo*> constraints = constraint->constraints();
  for (unsigned int i = 0; i < constraints.size(); i++){
    if (!hasConstraint(constraints[i])) addConstraint(constraints[i]);
  }
    // also reciprocate
  if (!(constraint->hasConstraint(this))) constraint->addConstraint(this);
}


bool 
TermInfo::hasConstraint(TermInfo* constraint) const{
  bool foundConstraint = false;
  for (unsigned int i = 0; i < m_constraints.size(); i++){
    if ((m_constraints[i]->fullName() == constraint->fullName())) foundConstraint = true;
  }
  return foundConstraint;
}


void
TermInfo::removeConstraint(TermInfo* constraint){
    // remove "constraint" as a constraint
  if (hasConstraint(constraint)){
    for (unsigned int i = 0; i < m_constraints.size(); i++){
      if ((m_constraints[i]->fullName() == constraint->fullName())){
        m_constraints.erase(m_constraints.begin()+i,m_constraints.begin()+i+1);
        i = m_constraints.size();
      }
    }
  }
    // remove "constraint"'s constraints as constraints 
  vector<TermInfo*> constraints = constraint->constraints();
  for (unsigned int i = 0; i < constraints.size(); i++){
    if (hasConstraint(constraints[i])) removeConstraint(constraints[i]);
  }
    // reciprocate
  if (constraint->hasConstraint(this)) constraint->removeConstraint(this);
}


void 
TermInfo::addParameter (ParameterInfo* parameter){
  bool foundParameter = false;
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    if (m_parameters[i]->parName() == parameter->parName()) foundParameter = true;
  }
  if (!foundParameter) m_parameters.push_back(parameter);
}

void
TermInfo::removeParameter (ParameterInfo* parameter){
  for (unsigned int i = 0; i < m_parameters.size(); i++){
    if (m_parameters[i]->parName() == parameter->parName()){
      m_parameters.erase(m_parameters.begin()+i,m_parameters.begin()+i+1);
      i = m_parameters.size();
    }
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

  vector< vector<string> > n_factors = factors();
  vector< TermInfo* >      n_constraints = constraints();
  vector< ParameterInfo* > n_parameters = parameters();

  cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++  AMPLITUDE INFO  +++++++++++++" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "      REACTION NAME:  " << m_reactionName << endl;
  cout << "  COHERENT SUM NAME:  " << m_sumName << endl;
  cout << "     AMPLITUDE NAME:  " << m_ampName << endl;
  cout << "            FACTORS:  " << n_factors.size() << endl;
  for (unsigned int i = 0; i < n_factors.size(); i++){
    vector<string> factor = n_factors[i];
    cout << "\t\t" << i+1 << ".  ";
    for (unsigned int j = 0; j < factor.size(); j++){
      cout << " " << factor[j];
    }
    cout << endl;
  }
  cout << " EXTRA PERMUTATIONS:  " << m_permutations.size() << endl;
  for (unsigned int i = 0; i < m_permutations.size(); i++){
    vector<int> permutation = m_permutations[i];
    cout << "\t\t" << i+1 << ".  ";
    for (unsigned int j = 0; j < permutation.size(); j++){
      cout << " " << permutation[j];
    }
    cout << endl;
  }
  cout << "        CONSTRAINTS:  " << n_constraints.size() << endl;
  for (unsigned int i = 0; i < n_constraints.size(); i++){
    cout << "\t\t" << i+1 << ".  " << n_constraints[i]->fullName() << endl;
  }
  cout << "         PARAMETERS:  " << n_parameters.size() << endl;
  for (unsigned int i = 0; i < n_parameters.size(); i++){
    cout << "\t\t" << i+1 << ".  " << n_parameters[i]->parName() << endl;
  }
  cout << "      INITIAL VALUE:  " << m_value << endl;
  cout << "               REAL?  " << m_real << endl;
  cout << "              FIXED?  " << m_fixed << endl;
  cout << "              SCALE:  " << m_scale << endl;

  if (fileName != ""){
    outfile.close();
    cout.rdbuf(cout_sbuf);
  }

}



void
PDFInfo::display(string fileName, bool append){

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

  vector< vector<string> > n_factors = factors();
  vector< TermInfo* >      n_constraints = constraints();
  vector< ParameterInfo* > n_parameters = parameters();

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "+++++++++++++  PDF INFO  +++++++++++++" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "      REACTION NAME:  " << m_reactionName << endl;
  cout << "           PDF NAME:  " << m_pdfName << endl;
  cout << "            FACTORS:  " << n_factors.size() << endl;
  for (unsigned int i = 0; i < n_factors.size(); i++){
    vector<string> factor = n_factors[i];
    cout << "\t\t" << i+1 << ".  ";
    for (unsigned int j = 0; j < factor.size(); j++){
      cout << " " << factor[j];
    }
    cout << endl;
  }
  cout << "        CONSTRAINTS:  " << n_constraints.size() << endl;
  for (unsigned int i = 0; i < n_constraints.size(); i++){
    cout << "\t\t" << i+1 << ".  " << n_constraints[i]->fullName() << endl;
  }
  cout << "         PARAMETERS:  " << n_parameters.size() << endl;
  for (unsigned int i = 0; i < n_parameters.size(); i++){
    cout << "\t\t" << i+1 << ".  " << n_parameters[i]->parName() << endl;
  }
  cout << "      INITIAL VALUE:  " << m_value << endl;
  cout << "              FIXED?  " << m_fixed << endl;
  cout << "              SCALE:  " << m_scale << endl;

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







void
ReactionInfo::clear(){
  vector<string> empty;
  m_data  = pair<string, vector<string> >("",empty);
  m_bkgnd = pair<string, vector<string> >("",empty);
  m_genMC = pair<string, vector<string> >("",empty);
  m_accMC = pair<string, vector<string> >("",empty);
  m_normIntFile = "";
  m_normIntFileInput = false;
  setGPUDeviceNumber();
}

void
TermInfo::termClear(){
  m_factors.clear();
  m_constraints.clear();
  m_parameters.clear();
}

void
AmplitudeInfo::clear(){
  termClear();
  m_value = complex<double>(0.0,0.0);
  m_real = false;
  m_fixed = false;
  m_scale = "1.0";
  m_permutations.clear();
}

void
PDFInfo::clear(){
  termClear();
  m_value = 0.0;
  m_fixed = false;
  m_scale = "1.0";
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





