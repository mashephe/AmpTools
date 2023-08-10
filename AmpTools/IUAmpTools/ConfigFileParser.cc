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
#include <fstream>
#include <cctype>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <stdlib.h>

#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/report.h"

const char* ConfigFileParser::kModule = "ConfigFileParser";
const char* ConfigFileLine::kModule = "ConfigFileLine";

bool ConfigFileParser::m_verboseParsing = false;


ConfigFileParser::ConfigFileParser()
  : m_configurationInfo( NULL ), m_fitName(""), m_configFile(""){
}



ConfigFileParser::ConfigFileParser(const string& configFile)
  : m_configurationInfo( NULL ), m_fitName(""){

  readConfigFile(configFile);

}



ConfigFileParser::ConfigFileParser(istream& input)
  : m_configurationInfo( NULL ), m_fitName(""){

  readConfigFile(input);

}


void
ConfigFileParser::readConfigFile(const string& configFile){
  if (m_configurationInfo) delete m_configurationInfo;
  m_configFile = configFile;
  m_fitName = "";

  m_configFileLines = expandConfigFileLines(readConfigFileLines(configFile));
  setupConfigurationInfo();

}


void
ConfigFileParser::readConfigFile(istream& input){
  if (m_configurationInfo) delete m_configurationInfo;
  m_configFile = "stream";
  m_fitName = "";

  m_configFileLines = expandConfigFileLines(readConfigFileLines(input));
  setupConfigurationInfo();

}


vector<ConfigFileLine>
ConfigFileParser::readConfigFileLines(const string& configfile) const{

  report( DEBUG, kModule ) << "Reading from the file...  " << configfile << endl;

  vector<ConfigFileLine> configFileLines;
  int lineNumber = 0;

  ifstream in(configfile.c_str());
  if (!in.is_open()){
    report( ERROR, kModule ) << "Could not open file: " << configfile << endl;
    exit(1);
  }
  while (!in.eof()){
    string line;
    getline(in,line);
    configFileLines.push_back(ConfigFileLine(configfile,++lineNumber, line));
      if (m_verboseParsing)
      configFileLines[configFileLines.size()-1].printLine();
  }
  in.close();

  report( DEBUG, kModule ) << "Finished reading from the file...  " << configfile << endl;

  return configFileLines;

}


vector<ConfigFileLine>
ConfigFileParser::readConfigFileLines(istream& input) const{

  report( DEBUG, kModule ) << "Reading from a stream...  " << endl;

  vector<ConfigFileLine> configFileLines;
  int lineNumber = 0;

  while (!input.eof()){
    string line;
    getline(input,line);
    configFileLines.push_back(ConfigFileLine("stream",++lineNumber, line));
      if (m_verboseParsing)
      configFileLines[configFileLines.size()-1].printLine();
  }

  report( DEBUG, kModule ) <<  "Finished reading from a stream  " << endl;

  return configFileLines;

}



vector<ConfigFileLine>
ConfigFileParser::expandConfigFileLines(vector<ConfigFileLine> configFileLines) const{

  report( DEBUG, kModule ) << "Starting to expand config file lines..." << endl;


    // flush all "include" lines

  for (vector<ConfigFileLine>::iterator lineItr = configFileLines.begin();
       lineItr != configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "include"){

      if (lineItr->arguments().size() != 1){
        report( ERROR, kModule ) << "Wrong number of arguments for the include statement:  " << endl;
        lineItr->printLine();
        exit(1);
      }

        // read in the input config file lines

      string inputFile = lineItr->arguments()[0];
      vector<ConfigFileLine> inputFileLines = readConfigFileLines(inputFile);

        // replace the input line

      vector<ConfigFileLine>::iterator inputBegin = inputFileLines.begin();
      vector<ConfigFileLine>::iterator inputEnd   = inputFileLines.end();
      *lineItr = ConfigFileLine("",0,"##########  INCLUDE FILE " + inputFile + "  ########");
      configFileLines.insert(++lineItr, inputBegin, inputEnd);
      lineItr = configFileLines.begin();       

    }

  }


    // flush all definitions given by "define" lines

  for (vector<ConfigFileLine>::iterator lineItr = configFileLines.begin();
       lineItr != configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "define"){

      if (lineItr->arguments().size() == 0){
        report( ERROR, kModule ) << "Wrong number of arguments for the define statement:" << endl;
        lineItr->printLine();
        exit(1);
      }

      vector<string> arguments = lineItr->arguments();
      string         key       = arguments[0];
      vector<string> value     (arguments.begin()+1,arguments.end()); 

      for (unsigned int j = 0; j < configFileLines.size(); j++){
        configFileLines[j].flushDefinition(key,value);
      }

    }
  }


    // save the "loop" lines

  map<string, vector<string> > mapLoops;
  vector<string> vectorLoops;
  for (vector<ConfigFileLine>::iterator lineItr = configFileLines.begin();
       lineItr != configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "loop"){

      if (lineItr->arguments().size() < 3){
        report( ERROR, kModule ) << "Wrong number of arguments for the loop statement:" << endl;
        lineItr->printLine();
        exit(1);
      }

      vector<string> arguments = lineItr->arguments();
      string         key       = arguments[0];
      vector<string> value     (arguments.begin()+1,arguments.end()); 

      vectorLoops.push_back(key);
      mapLoops[key] = value;

    }
  }


    // expand the "loop" definitions

  for (vector<ConfigFileLine>::iterator lineItr = configFileLines.begin();
       lineItr != configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "loop") continue;
    if (lineItr->comment()) continue;
    vector<string> arguments = lineItr->arguments();
    if (arguments.size() < 2) continue;

    vector<string> foundLoops;
    for (unsigned int i = 0; i < arguments.size(); i++){
      for (unsigned int j = 0; j < vectorLoops.size(); j++){
        if (arguments[i] == vectorLoops[j]) foundLoops.push_back(vectorLoops[j]);
      }
    }

    if (foundLoops.size() == 0) continue;
    if (foundLoops.size() >= 2){
      for (unsigned int i = 0; i < foundLoops.size()-1; i++){
        for (unsigned int j = i+1; j < foundLoops.size(); j++){
          if (mapLoops[foundLoops[i]].size() != mapLoops[foundLoops[j]].size()){
            report( ERROR, kModule ) << "Mismatch in loop size:" << endl;
            lineItr->printLine();
            exit(1);
          }
        }
      }
    }

    vector<ConfigFileLine> newLines;
    for (unsigned int i = 0; i < mapLoops[foundLoops[0]].size(); i++){
      newLines.push_back(*lineItr);
    }

    for (unsigned int i = 0; i < newLines.size(); i++){
      for (unsigned int j = 0; j < foundLoops.size(); j++){
        vector<string> value;  value.push_back(mapLoops[foundLoops[j]][i]);
        newLines[i].flushDefinition(foundLoops[j],value);
      }
    }

    vector<ConfigFileLine>::iterator inputBegin = newLines.begin();
    vector<ConfigFileLine>::iterator inputEnd   = newLines.end();
    *lineItr = ConfigFileLine("",0,"##########  INSERTING A LOOP  ########");
    configFileLines.insert(++lineItr, inputBegin, inputEnd);
    lineItr = configFileLines.begin();       

  }

  report( DEBUG, kModule ) << "Finished expanding config file lines..." << endl;

    if (m_verboseParsing){
      for (unsigned int i = 0; i < configFileLines.size(); i++){
        configFileLines[i].printLine();
      }
    }

  return configFileLines;

}





void
ConfigFileParser::setupConfigurationInfo(){


  // Do some quick syntax checks

  report( DEBUG, kModule ) << "Begin Syntax Checking " << endl;

  checkSyntax();

  report( DEBUG, kModule ) << "Finished Syntax Checking " << endl;


  // ZEROTH PASS ("fit")
  
  report( DEBUG, kModule ) << "Starting ZEROTH PASS (finding the fit name)" << endl;
  
  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if ((*lineItr).keyword() == "fit") doFit(*lineItr);

    if ((*lineItr).keyword() == "keyword") m_userKeywords.insert((*lineItr).arguments()[0]);

  }

  report( DEBUG, kModule ) << "Finished ZEROTH PASS" << endl;


    // Create a new ConfigurationInfo object


    if (m_fitName == "")
      report( WARNING, kModule ) << "use the keyword \"fit\" to define a fit name" << endl;

  m_configurationInfo = new ConfigurationInfo(m_fitName);


      // FIRST PASS ("reaction")
  
  report( DEBUG, kModule ) << "Starting FIRST PASS (creating reactions and parameters)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if ((*lineItr).keyword() == "reaction") doReaction(*lineItr);
    if ((*lineItr).keyword() == "parameter") doParameter(*lineItr);

    if (m_userKeywords.count((*lineItr).keyword())) doKeyword(*lineItr);

  }

  report( DEBUG, kModule ) << "Finished FIRST PASS" << endl;



      // SECOND PASS ("sum" and reaction info)
  
  report( DEBUG, kModule ) << "Starting SECOND PASS (creating sums and filling reactions)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if (((*lineItr).keyword() == "datafile")  ||
        ((*lineItr).keyword() == "genmcfile") ||
        ((*lineItr).keyword() == "accmcfile") ||
        ((*lineItr).keyword() == "data") ||
        ((*lineItr).keyword() == "bkgnd") ||
        ((*lineItr).keyword() == "genmc") ||
        ((*lineItr).keyword() == "accmc"))  doData(*lineItr);

    if ((*lineItr).keyword() == "normintfile") doNormInt(*lineItr);

    if ((*lineItr).keyword() == "sum") doSum(*lineItr);

    if ((*lineItr).keyword() == "gpudevice") doGPUDevice(*lineItr);

  }

  report( DEBUG, kModule ) << "Finished SECOND PASS" << endl;



      // THIRD PASS ("amplitude")
  
  report( DEBUG, kModule ) << "Starting THIRD PASS (creating amplitudes)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if ((*lineItr).keyword() == "amplitude") doAmplitude(*lineItr);
    if ((*lineItr).keyword() == "pdf")       doPDF(*lineItr);

  }

  report( DEBUG, kModule ) << "Finished THIRD PASS" << endl;



      // FOURTH PASS (operations on amplitudes -- constrain, permute, scale)
  
  report( DEBUG, kModule ) << "Starting FOURTH PASS (filling amplitudes 1/2)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if ((*lineItr).keyword() == "constrain") doConstrain(*lineItr);
    if ((*lineItr).keyword() == "permute") doPermute(*lineItr);
    if ((*lineItr).keyword() == "scale") doScale(*lineItr);
    if ((*lineItr).keyword() == "pdfconstrain") doPDFConstrain(*lineItr);
    if ((*lineItr).keyword() == "pdfscale") doPDFScale(*lineItr);

  }

  report( DEBUG, kModule ) << "Finished FOURTH PASS" << endl;


      // FIFTH PASS (operations on amplitudes -- initialize)
  
  report( DEBUG, kModule ) << "Starting FIFTH PASS (filling amplitudes 2/2)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if ((*lineItr).keyword() == "initialize") doInitialize(*lineItr);
    if ((*lineItr).keyword() == "pdfinitialize") doPDFInitialize(*lineItr);

  }

  report( DEBUG, kModule ) << "Finished FIFTH PASS" << endl;


      // Do some quick syntax checks

  report( DEBUG, kModule ) << "Begin Syntax Checking " << endl;

  checkSyntax();

  report( DEBUG, kModule ) << "Finished Syntax Checking " << endl;

}






void
ConfigFileParser::checkSyntax() const{

  map<string, pair<int,int> > keywordParameters;
  keywordParameters["define"]        = pair<int,int>(1,500);
  keywordParameters["loop"]          = pair<int,int>(3,100);
  keywordParameters["keyword"]       = pair<int,int>(3,3);
  keywordParameters["fit"]           = pair<int,int>(1,1);
  keywordParameters["reaction"]      = pair<int,int>(3,100);
  keywordParameters["data"]          = pair<int,int>(2,100);
  keywordParameters["bkgnd"]         = pair<int,int>(2,100);
  keywordParameters["genmc"]         = pair<int,int>(2,100);
  keywordParameters["accmc"]         = pair<int,int>(2,100);
  keywordParameters["normintfile"]   = pair<int,int>(2,3);
  keywordParameters["sum"]           = pair<int,int>(2,100);
  keywordParameters["amplitude"]     = pair<int,int>(4,500);
  keywordParameters["initialize"]    = pair<int,int>(6,8);
  keywordParameters["constrain"]     = pair<int,int>(6,100);
  keywordParameters["permute"]       = pair<int,int>(5,100);
  keywordParameters["scale"]         = pair<int,int>(4,4);
  keywordParameters["pdf"]           = pair<int,int>(3,100);
  keywordParameters["pdfinitialize"] = pair<int,int>(3,4);
  keywordParameters["pdfconstrain"]  = pair<int,int>(4,100);
  keywordParameters["pdfscale"]      = pair<int,int>(3,3);
  keywordParameters["parameter"]     = pair<int,int>(2,5);
  keywordParameters["gpudevice"]     = pair<int,int>(2,2);
    // these are deprecated, but print out an error message later
  keywordParameters["datafile"]      = pair<int,int>(2,100);
  keywordParameters["genmcfile"]     = pair<int,int>(2,100);
  keywordParameters["accmcfile"]     = pair<int,int>(2,100);

  for (vector<ConfigFileLine>::const_iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (!lineItr->comment()){
      map<string, pair<int,int> >::const_iterator mapItr = keywordParameters.find(lineItr->keyword());
      if (mapItr == keywordParameters.end()){
        report( ERROR, kModule ) << "Undefined keyword:  " << lineItr->keyword() << endl;
        lineItr->printLine();
        exit(1);
      }
      else if (((int)lineItr->arguments().size() > mapItr->second.second) ||
               ((int)lineItr->arguments().size() < mapItr->second.first)){
        report( ERROR, kModule ) << "Keyword " << lineItr->keyword() <<
        " has the wrong number of arguments: " << endl;
        lineItr->printLine();
        exit(1);
      }
      else if (lineItr->keyword() == "keyword"){
        keywordParameters[lineItr->arguments()[0]] = 
        pair<int,int>(atoi(lineItr->arguments()[1].c_str()),atoi(lineItr->arguments()[2].c_str()));
      }
    }
  }

}





void
ConfigFileParser::doFit(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  m_fitName = arguments[0];
}



void
ConfigFileParser::doReaction(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  vector<string> particles (arguments.begin()+1, arguments.end());
  m_configurationInfo->createReaction(reaction,particles);
}


void
ConfigFileParser::doKeyword(const ConfigFileLine& line){
  m_configurationInfo->addUserKeyword(line.keyword(),line.arguments());
}



void
ConfigFileParser::doData(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction  = arguments[0];
  string classname = arguments[1];
  vector<string> dataargs (arguments.begin()+2, arguments.end());
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    report( ERROR, kModule ) << "Can't associate data with a reaction:  " << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "datafile"){
    report( ERROR, kModule ) << "datafile is deprecated, use data" << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "genmcfile"){
    report( ERROR, kModule ) << "genmcfile is deprecated, use genmc" << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "accmcfile"){
    report( ERROR, kModule ) << "accmcfile is deprecated, use accmc" << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "data")    rct->setData (classname, dataargs);
  if (line.keyword() == "bkgnd")   rct->setBkgnd(classname, dataargs);
  if (line.keyword() == "genmc")   rct->setGenMC(classname, dataargs);
  if (line.keyword() == "accmc")   rct->setAccMC(classname, dataargs);
}


void
ConfigFileParser::doGPUDevice(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction  = arguments[0];
  int deviceNumber = atoi((arguments[1]).c_str());
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    report( ERROR, kModule ) << "Can't associate gpudevice with a reaction:  " << endl;
    line.printLine();
    exit(1);
  }
  rct->setGPUDeviceNumber(deviceNumber);
}


void
ConfigFileParser::doNormInt(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string file = arguments[1];
  bool input = false;
  if (arguments.size() > 2 && arguments[2] == "input") input = true;
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    report( ERROR, kModule ) << "Can't associate normintfile with a reaction:  " << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "normintfile"){
    rct->setNormIntFile(file, input);
  }
}



void
ConfigFileParser::doSum(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    report( ERROR, kModule ) << "Can't associate sum with a reaction:  " << endl;
    line.printLine();
    exit(1);
  }
  for (unsigned int i = 1; i < arguments.size(); i++){
    string sum = arguments[i];
    m_configurationInfo->createCoherentSum(reaction,sum);
  }
}


void
ConfigFileParser::doParameter(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string parname = arguments[0];
  double value = atof((arguments[1]).c_str());
  ParameterInfo* parinfo = m_configurationInfo->createParameter(parname,value);
  string type("floating");
  double a = 0.0;
  double b = 0.0;
  if (arguments.size() > 2) type  = arguments[2];
  if (arguments.size() > 3) a     = atof((arguments[3]).c_str());
  if (arguments.size() > 4) b     = atof((arguments[4]).c_str());
  if (type == "floating"){
  }
  else if (type == "fixed"){
    parinfo->setFixed(true);
  }
  else if (type == "bounded"){
    parinfo->setBounded(true);
    parinfo->setLowerBound(a);
    parinfo->setUpperBound(b);
  }
  else if (type == "gaussian"){
    parinfo->setGaussianBounded(true);
    parinfo->setCentralValue(a);
    parinfo->setGaussianError(b);
  }
  else{
    report( ERROR, kModule ) << "parameter type must be floating, fixed, bounded, or gaussian  " << endl;
    line.printLine();
    exit(1);
  }    
}


void
ConfigFileParser::doAmplitude(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction  = arguments[0];
  string sumname   = arguments[1];
  string ampname   = arguments[2];
  vector<string> ampargs (arguments.begin()+3, arguments.end());
  AmplitudeInfo* ampinfo = m_configurationInfo->amplitude(reaction,sumname,ampname);
  if (!ampinfo) ampinfo = m_configurationInfo->createAmplitude(reaction,sumname,ampname);
  ampinfo->addFactor(ampargs);
  for (unsigned int i = 1; i < ampargs.size(); i++){
    unsigned int j = ampargs[i].size()-1;
    if ((ampargs[i][0] == '[') && (ampargs[i][j] == ']')){
      string parname("");
      for (unsigned int k = 1; k < j; k++){
        parname += ampargs[i][k];
      }
      ParameterInfo* parinfo = m_configurationInfo->parameter(parname);
      if (!parinfo){
        report( ERROR, kModule ) << "can't find parameter " << parname << endl;
        line.printLine();
        exit(1);
      }
      ampinfo->addParameter(parinfo);
    }
  }
}


void
ConfigFileParser::doPDF(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction  = arguments[0];
  string pdfname   = arguments[1];
  vector<string> pdfargs (arguments.begin()+2, arguments.end());
  PDFInfo* pdfinfo = m_configurationInfo->pdf(reaction,pdfname);
  if (!pdfinfo) pdfinfo = m_configurationInfo->createPDF(reaction,pdfname);
  pdfinfo->addFactor(pdfargs);
  for (unsigned int i = 1; i < pdfargs.size(); i++){
    unsigned int j = pdfargs[i].size()-1;
    if ((pdfargs[i][0] == '[') && (pdfargs[i][j] == ']')){
      string parname("");
      for (unsigned int k = 1; k < j; k++){
        parname += pdfargs[i][k];
      }
      ParameterInfo* parinfo = m_configurationInfo->parameter(parname);
      if (!parinfo){
        report( ERROR, kModule ) << "can't find parameter " << parname << endl;
        line.printLine();
        exit(1);
      }
      pdfinfo->addParameter(parinfo);
    }
  }
}


void
ConfigFileParser::doConstrain(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  if (arguments.size()%3 != 0){
    report( ERROR, kModule ) << "wrong number of arguments for constrain keyword " << endl;
    line.printLine();
    exit(1);
  }
  string reaction1  = arguments[0];
  string sumname1   = arguments[1];
  string ampname1   = arguments[2];
  for (unsigned int i = 1; i < arguments.size()/3; i++){
    string reaction2 = arguments[i*3];
    string sumname2  = arguments[i*3+1];
    string ampname2  = arguments[i*3+2];
    AmplitudeInfo* amplitude1 = m_configurationInfo->amplitude(reaction1,sumname1,ampname1);
    AmplitudeInfo* amplitude2 = m_configurationInfo->amplitude(reaction2,sumname2,ampname2);
    if ((!amplitude1) || (!amplitude2)){
      report( ERROR, kModule ) << "trying to constrain nonexistent amplitude " << endl;
      line.printLine();
      exit(1);
    }
    amplitude1->addConstraint(amplitude2);
  }
}

void
ConfigFileParser::doPDFConstrain(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  if (arguments.size()%2 != 0){
    report( ERROR, kModule ) << "wrong number of arguments for pdfconstrain keyword " << endl;
    line.printLine();
    exit(1);
  }
  string reaction1  = arguments[0];
  string pdfname1   = arguments[1];
  for (unsigned int i = 1; i < arguments.size()/2; i++){
    string reaction2 = arguments[i*2];
    string pdfname2  = arguments[i*2+1];
    PDFInfo* pdf1 = m_configurationInfo->pdf(reaction1,pdfname1);
    PDFInfo* pdf2 = m_configurationInfo->pdf(reaction2,pdfname2);
    if ((!pdf1) || (!pdf2)){
      report( ERROR, kModule ) << "trying to constrain nonexistent pdf " << endl;
      line.printLine();
      exit(1);
    }
    pdf1->addConstraint(pdf2);
  }
}


void
ConfigFileParser::doPermute(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string sumname  = arguments[1];
  string ampname  = arguments[2];
  vector<string> permutation(arguments.begin()+3, arguments.end());
  AmplitudeInfo* amplitude = m_configurationInfo->amplitude(reaction,sumname,ampname);
  if (!amplitude){
    report( ERROR, kModule ) << "trying to permute nonexistent amplitude " << endl;
    line.printLine();
    exit(1);
  }
  vector<string> particleList = m_configurationInfo->reaction(reaction)->particleList();
  if (permutation.size() != particleList.size()){
    report( ERROR, kModule ) << "wrong number of arguments for permute keyword " << endl;
    line.printLine();
    exit(1);
  }
  vector<int> intpermutation;
  for (unsigned int i = 0; i < permutation.size(); i++){
    for (unsigned int j = 0; j < permutation[i].size(); j++){
      if (!isdigit(permutation[i][j])){
        report( ERROR, kModule ) << "particle index is not an unsigned integer " << endl;
        line.printLine();
        exit(1);
      }
    }
    int ipart = atoi(permutation[i].c_str());
    if (ipart < 0 || ipart >= (int)particleList.size()){
      report( ERROR, kModule ) << "particle index is out of bounds " << endl;
      line.printLine();
      exit(1);
    }
    intpermutation.push_back(ipart);
  }
  for (unsigned int i = 0; i < intpermutation.size(); i++){
    for (unsigned int j = i+1; j < intpermutation.size(); j++){
      if (intpermutation[i] == intpermutation[j]){
        report( ERROR, kModule ) << "particle index repeated " << endl;
        line.printLine();
        exit(1);
      }
    }
  }
  amplitude->addPermutation(intpermutation);
}


void
ConfigFileParser::doInitialize(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string sumname  = arguments[1];
  string ampname  = arguments[2];
  string type     = arguments[3];
  double value1   = atof(arguments[4].c_str());  
  double value2   = atof(arguments[5].c_str());  
  string fixtype1("floating");
  string fixtype2("");
  if (arguments.size() >= 7) fixtype1 = arguments[6];
  if (arguments.size() == 8) fixtype2 = arguments[7];
  AmplitudeInfo* amplitude = m_configurationInfo->amplitude(reaction,sumname,ampname);
  if (!amplitude){
    report( ERROR, kModule ) << "trying to initialize nonexistent amplitude " << endl;
    line.printLine();
    exit(1);
  }
  if (type == "cartesian"){
    amplitude->setValue(complex<double>(value1,value2));
  }
  else if (type == "polar"){
    amplitude->setValue(polar(value1,value2));
  }
  else if (type == "events"){
    report( ERROR, kModule ) << "initializing with events is not yet re-implemented " << endl;
    line.printLine();
    exit(1);
    /*
    string normIntFile = m_configurationInfo->reaction(reaction)->normIntFile();
    if (normIntFile == ""){
     report( ERROR, kModule ) << "ConfigFileParser ERROR:  initializing with events, but no normalization integral " << endl;
     report( ERROR, kModule ) << "                         file has been specified for this reaction " << endl;
      line.printLine();
      exit(1);
    }
    ifstream test(normIntFile.c_str());
    if (!test){
     report( ERROR, kModule ) << "ConfigFileParser ERROR:  initializing with events, but can't find the normalization integral " << endl;
     report( ERROR, kModule ) << "                         file specified for this reaction (" << normIntFile << ")" << endl;
      line.printLine();
      exit(1);
    }
    NormIntInterface normint(normIntFile);
    if (!normint.hasAmpInt(ampname,ampname)){
     report( ERROR, kModule ) << "ConfigFileParser ERROR:  can't find the right amplitude in the normalization integral file " << endl;
      line.printLine();
      exit(1);
    }
    value1 = sqrt(value1/(abs(normint.ampInt(ampname,ampname))));
    amplitude->setValue(polar(value1,value2));
    */
  }
  else{
    report( ERROR, kModule ) << "initialize must use cartesian, polar, or events  " << endl;
    line.printLine();
    exit(1);
  }
  if (fixtype1 == "floating") {
  }
  else if (fixtype1 == "real"){
    amplitude->setReal(true);
  }
  else if (fixtype1 == "fixed"){
    amplitude->setFixed(true);
  }
  else{
    report( ERROR, kModule ) << "initialize must use floating, fixed, or real  " << endl;
    line.printLine();
    exit(1);
  }
  if (fixtype2 == "") {
  }
  else if (fixtype2 == "real"){
    amplitude->setReal(true);
  }
  else if (fixtype2 == "fixed"){
    amplitude->setFixed(true);
  }
  else if (fixtype2 == "floating"){
    amplitude->setReal(false);
    amplitude->setFixed(false);
  }
  else{
    report( ERROR, kModule ) << "initialize must use floating, fixed, or real  " << endl;
    line.printLine();
    exit(1);
  }
}



void
ConfigFileParser::doPDFInitialize(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string pdfname  = arguments[1];
  double value    = atof(arguments[2].c_str());  
  string fixtype("floating");
  if (arguments.size() >= 4) fixtype = arguments[3];
  PDFInfo* pdf = m_configurationInfo->pdf(reaction,pdfname);
  if (!pdf){
    report( ERROR, kModule ) << "trying to initialize nonexistent pdf " << endl;
    line.printLine();
    exit(1);
  }
  pdf->setValue(value);
       if (fixtype == "floating"){ pdf->setFixed(false); }
  else if (fixtype == "fixed")   { pdf->setFixed(true); }
  else{
    report( ERROR, kModule ) << "pdfinitialize must use floating or fixed  " << endl;
    line.printLine();
    exit(1);
  }
}


void
ConfigFileParser::doScale(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string sumname  = arguments[1];
  string ampname  = arguments[2];
  string value    = arguments[3];
  AmplitudeInfo* amplitude = m_configurationInfo->amplitude(reaction,sumname,ampname);
  if (!amplitude){
    report( ERROR, kModule ) << "trying to scale nonexistent amplitude " << endl;
    line.printLine();
    exit(1);
  }
  if ((value.size() > 0) && (value[0] == '[') and (value[value.size()-1] == ']')){
    string parname("");
    for (unsigned int k = 1; k < value.size()-1; k++){
      parname += value[k];
    }
    ParameterInfo* parinfo = m_configurationInfo->parameter(parname);
    if (!parinfo){
      report( ERROR, kModule ) << "can't find parameter " << parname << endl;
      line.printLine();
      exit(1);
    }
    amplitude->addParameter(parinfo);
  }
  amplitude->setScale(value);
}


void
ConfigFileParser::doPDFScale(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  string pdfname  = arguments[1];
  string value    = arguments[2];
  PDFInfo* pdf = m_configurationInfo->pdf(reaction,pdfname);
  if (!pdf){
    report( ERROR, kModule ) << "trying to scale nonexistent pdf " << endl;
    line.printLine();
    exit(1);
  }
  if ((value.size() > 0) && (value[0] == '[') and (value[value.size()-1] == ']')){
    string parname("");
    for (unsigned int k = 1; k < value.size()-1; k++){
      parname += value[k];
    }
    ParameterInfo* parinfo = m_configurationInfo->parameter(parname);
    if (!parinfo){
      report( ERROR, kModule ) << "can't find parameter " << parname << endl;
      line.printLine();
      exit(1);
    }
    pdf->addParameter(parinfo);
  }
  pdf->setScale(value);
}


void
ConfigFileParser::displayConfigFile() const{
  for (unsigned int i = 0; i < m_configFileLines.size(); i++){
    m_configFileLines[i].printLine();
  }
}












ConfigFileLine::ConfigFileLine(const string& fileName, int lineNumber, const string& line){

  m_fileName   = fileName;
  m_lineNumber = lineNumber;
  m_line       = line;
  m_keyword    = "";
  m_comment    = false;

    // replace all "::" with " "
  string newline(line);
  while (newline.find("::") != string::npos){
    newline.replace(newline.find("::"),2," ");
  }

    // parse the line into words
  vector<string> words;
  string word("");
  for (unsigned int j = 0; j < newline.size(); j++){
    if (!isspace(newline[j])){
      word += newline[j];
      if ((j == (newline.size()-1))&&(!word.empty())){
        words.push_back(word);
        word = "";
      }
    } 
    else if (!word.empty()){
      words.push_back(word);
      word = "";
    }
  }

    // record the keyword
  if (words.size() > 0) m_keyword = words[0];

    // check if this is a comment line
  if (m_keyword.empty())        m_comment = true;
  else if (m_keyword[0] == '#') m_comment = true;
  if (m_comment) m_keyword = "";

    // record the arguments
  if (!m_comment){
    for (unsigned int i = 1; i < words.size(); i++){
      m_arguments.push_back(words[i]);
    }
  }
}


void
ConfigFileLine::flushDefinition(const string& word, const vector<string>& definition){
  if (comment()) return;
  bool hasWord = false;
  for (unsigned int i = 0; i < m_arguments.size(); i++){ 
    if (m_arguments[i] == word) hasWord = true;
  }
  if (!hasWord) return;
  if ((m_keyword == "define") && (m_arguments[0] == word)) return;
  vector<string> newArguments;
  for (unsigned int i = 0; i < m_arguments.size(); i++){
    if (m_arguments[i] != word){
      newArguments.push_back(m_arguments[i]);
    }
    else{
      for (unsigned int j = 0; j < definition.size(); j++){
        newArguments.push_back(definition[j]);
      }
    }
  }
  m_arguments = newArguments;
}

void
ConfigFileLine::printLine() const {
  
  report( INFO, kModule ) << "(" << m_lineNumber << ") " <<
  m_fileName   << " >>   " <<
  m_line       << endl;
}

void
ConfigFileLine::printArguments() const {
  
  report( INFO, kModule ) << "KEYWORD:  " << m_keyword << endl;
  report( INFO, kModule ) << "ARGUMENTS: " << endl;
  for ( unsigned int i = 0; i < m_arguments.size(); i++){

    report( INFO, kModule ) << m_arguments[i] << endl;
  }
}
