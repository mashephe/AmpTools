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


ConfigFileParser::ConfigFileParser(const string& configFile, bool verboseParsing)
  : m_configurationInfo( NULL ), m_fitName(configFile)
{


      // Initial Setup

    cout << endl;
    cout << "ConfigFileParser INFO:  Setting up new ConfigFileParser with file...  " 
         << configFile << endl;

  m_configFile = configFile;
  m_configFileLines.clear();
  m_verboseParsing = verboseParsing;



      // Read the config file and expand the include files

    cout << "ConfigFileParser INFO:  Begin Parsing Files " << endl;

  getConfigFileLines();

    cout << "ConfigFileParser INFO:  Finished Parsing Files " << endl;



      // Do some quick syntax checks

    cout << "ConfigFileParser INFO:  Begin Syntax Checking " << endl;

  checkSyntax();

    cout << "ConfigFileParser INFO:  Finished Syntax Checking " << endl;



      // Fill the ConfigurationInfo object

    cout << "ConfigFileParser INFO:  Begin Filling the ConfigurationInfo Object " << endl;

  getConfigurationInfo();

    cout << "ConfigFileParser INFO:  Finished Filling the ConfigurationInfo Object " << endl;



    cout << "ConfigFileParser INFO:  Finished setting up new ConfigFileParser" << endl;
    cout << endl;

}



ConfigFileParser::~ConfigFileParser(){

  if( m_configurationInfo != NULL )
    delete m_configurationInfo;
}




ConfigurationInfo*
ConfigFileParser::getConfigurationInfo(){


    // Check if the ConfigurationInfo object has already been filled.

  if (m_configurationInfo) return m_configurationInfo;



      // ZEROTH PASS ("fit")
  
    cout << "ConfigFileParser INFO:  Starting ZEROTH PASS (finding the fit name)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (m_verboseParsing) (*lineItr).printLine();

    if ((*lineItr).keyword() == "fit") doFit(*lineItr);

  }

    cout << "ConfigFileParser INFO:  Finished ZEROTH PASS" << endl;


    // Create a new ConfigurationInfo object

  m_configurationInfo = new ConfigurationInfo(m_fitName);



      // FIRST PASS ("reaction")
  
    cout << "ConfigFileParser INFO:  Starting FIRST PASS (creating reactions and parameters)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (m_verboseParsing) (*lineItr).printLine();

    if ((*lineItr).keyword() == "reaction") doReaction(*lineItr);
    if ((*lineItr).keyword() == "parameter") doParameter(*lineItr);

  }

    cout << "ConfigFileParser INFO:  Finished FIRST PASS" << endl;



      // SECOND PASS ("sum" and reaction info)
  
    cout << "ConfigFileParser INFO:  Starting SECOND PASS (creating sums and filling reactions)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (m_verboseParsing) (*lineItr).printLine();

    if (((*lineItr).keyword() == "datafile")  ||
        ((*lineItr).keyword() == "genmcfile") ||
        ((*lineItr).keyword() == "accmcfile") ||
        ((*lineItr).keyword() == "normintfile")) doFileName(*lineItr);

    if ((*lineItr).keyword() == "sum") doSum(*lineItr);

  }

    cout << "ConfigFileParser INFO:  Finished SECOND PASS" << endl;



      // THIRD PASS ("amplitude")
  
    cout << "ConfigFileParser INFO:  Starting THIRD PASS (creating amplitudes)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (m_verboseParsing) (*lineItr).printLine();

    if ((*lineItr).keyword() == "amplitude") doAmplitude(*lineItr);

  }

    cout << "ConfigFileParser INFO:  Finished THIRD PASS" << endl;



      // FOURTH PASS (operations on amplitudes)
  
    cout << "ConfigFileParser INFO:  Starting FOURTH PASS (filling amplitudes)" << endl;

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (m_verboseParsing) (*lineItr).printLine();

    if ((*lineItr).keyword() == "constrain") doConstrain(*lineItr);

    if ((*lineItr).keyword() == "permute") doPermute(*lineItr);

    if ((*lineItr).keyword() == "initialize") doInitialize(*lineItr);

    if ((*lineItr).keyword() == "scale") doScale(*lineItr);

  }

    cout << "ConfigFileParser INFO:  Finished FOURTH PASS" << endl;


      // Return the filled ConfigurationInfo object

  return m_configurationInfo;

}






vector<ConfigFileLine>
ConfigFileParser::getConfigFileLines(const string& configfile) const{

  cout << "ConfigFileParser INFO:  Begin parsing the file...  " << configfile << endl;

  vector<ConfigFileLine> configFileLines;
  int lineNumber = 0;

  ifstream in(configfile.c_str());
  if (!in.is_open()){
    cout << "ConfigFileParser ERROR:  Could not open file: " << configfile << endl;
    exit(1);
  }
  while (!in.eof()){
    string line;
    getline(in,line);
    configFileLines.push_back(ConfigFileLine(configfile,++lineNumber, line));
  }
  in.close();

  return configFileLines;

}


vector<ConfigFileLine>
ConfigFileParser::getConfigFileLines(){


    // check if the file is already parsed

  if (m_configFileLines.size() != 0) return m_configFileLines;


    // read in the initial config file

  m_configFileLines = getConfigFileLines(m_configFile);


    // flush all "include" lines

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "include"){

      if (lineItr->arguments().size() != 1){
        cout << "ConfigFileParser ERROR:  Wrong number of arguments for the include statement:" << endl;
        lineItr->printLine();
        exit(1);
      }

        // read in the input config file lines

      string inputFile = lineItr->arguments()[0];
      vector<ConfigFileLine> inputFileLines = getConfigFileLines(inputFile);

        // replace the input line

      vector<ConfigFileLine>::iterator inputBegin = inputFileLines.begin();
      vector<ConfigFileLine>::iterator inputEnd   = inputFileLines.end();
      *lineItr = ConfigFileLine("",0,"##########  INCLUDE FILE " + inputFile + "  ########");
      m_configFileLines.insert(++lineItr, inputBegin, inputEnd);
      lineItr = m_configFileLines.begin();       

    }

  }


    // flush all definitions given by "define" lines

  for (vector<ConfigFileLine>::iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){

    if (lineItr->keyword() == "define"){

      if (lineItr->arguments().size() == 0){
        cout << "ConfigFileParser ERROR:  Wrong number of arguments for the define statement:" << endl;
        lineItr->printLine();
        exit(1);
      }

      vector<string> arguments = lineItr->arguments();
      string         key       = arguments[0];
      vector<string> value     (arguments.begin()+1,arguments.end()); 

      for (unsigned int j = 0; j < m_configFileLines.size(); j++){
        m_configFileLines[j].flushDefinition(key,value);
      }

    }
  }

  return m_configFileLines;

}




void
ConfigFileParser::checkSyntax() const{

  map<string, pair<int,int> > keywordParameters;
  keywordParameters["define"]        = pair<int,int>(1,100);
  keywordParameters["keyword"]       = pair<int,int>(3,3);
  keywordParameters["fit"]           = pair<int,int>(1,1);
  keywordParameters["reaction"]      = pair<int,int>(3,100);
  keywordParameters["datafile"]      = pair<int,int>(2,100);
  keywordParameters["genmcfile"]     = pair<int,int>(2,100);
  keywordParameters["accmcfile"]     = pair<int,int>(2,100);
  keywordParameters["normintfile"]   = pair<int,int>(2,2);
  keywordParameters["sum"]           = pair<int,int>(2,100);
  keywordParameters["amplitude"]     = pair<int,int>(4,100);
  keywordParameters["initialize"]    = pair<int,int>(6,7);
  keywordParameters["constrain"]     = pair<int,int>(6,100);
  keywordParameters["permute"]       = pair<int,int>(5,100);
  keywordParameters["parameter"]     = pair<int,int>(2,5);
  keywordParameters["scale"]         = pair<int,int>(4,4);

  for (vector<ConfigFileLine>::const_iterator lineItr = m_configFileLines.begin();
       lineItr != m_configFileLines.end(); ++lineItr){
    if (!lineItr->comment()){
      map<string, pair<int,int> >::const_iterator mapItr = keywordParameters.find(lineItr->keyword());
      if (mapItr == keywordParameters.end()){
	cout << "ConfigFileParser ERROR:  Undefined keyword:  " << lineItr->keyword() << endl;
	lineItr->printLine();
	exit(1);
      }
      else if (((int)lineItr->arguments().size() > mapItr->second.second) ||
               ((int)lineItr->arguments().size() < mapItr->second.first)){
	cout << "ConfigFileParser ERROR:  Keyword " << lineItr->keyword() << 
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
ConfigFileParser::doFileName(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  vector<string> files (arguments.begin()+1, arguments.end());
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    cout << "ConfigFileParser ERROR:  Can't associate filename with a reaction:  " << endl;
    line.printLine();
    exit(1);
  }
  if (line.keyword() == "datafile"){
    for (unsigned int i = 0; i < files.size(); i++){
      rct->addDataFile(files[i]);
    }
  }
  if (line.keyword() == "genmcfile"){
    for (unsigned int i = 0; i < files.size(); i++){
      rct->addGenMCFile(files[i]);
    }
  }
  if (line.keyword() == "accmcfile"){
    for (unsigned int i = 0; i < files.size(); i++){
      rct->addAccMCFile(files[i]);
    }
  }
  if (line.keyword() == "normintfile"){
    for (unsigned int i = 0; i < files.size(); i++){
      rct->setNormIntFile(files[i]);
    }
  }
}



void
ConfigFileParser::doSum(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  string reaction = arguments[0];
  ReactionInfo* rct = m_configurationInfo->reaction(reaction);
  if (!rct){
    cout << "ConfigFileParser ERROR:  Can't associate sum with a reaction:  " << endl;
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
    cout << "ConfigFileParser ERROR:  parameter type must be floating, fixed, bounded, or gaussian  " << endl;
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
        cout << "ConfigFileParser ERROR:  can't find parameter " << parname << endl;
        line.printLine();
        exit(1);
      }
      ampinfo->addParameter(parinfo);
    }
  }
}


void
ConfigFileParser::doConstrain(const ConfigFileLine& line){
  vector<string> arguments = line.arguments();
  if (arguments.size()%3 != 0){
    cout << "ConfigFileParser ERROR:  wrong number of arguments for constrain keyword " << endl;
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
      cout << "ConfigFileParser ERROR:  trying to constrain nonexistent amplitude " << endl;
      line.printLine();
      exit(1);
    }
    amplitude1->addConstraint(amplitude2);
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
    cout << "ConfigFileParser ERROR:  trying to permute nonexistent amplitude " << endl;
    line.printLine();
    exit(1);
  }
  vector<string> particleList = m_configurationInfo->reaction(reaction)->particleList();
  if (permutation.size() != particleList.size()){
    cout << "ConfigFileParser ERROR:  wrong number of arguments for permute keyword " << endl;
    line.printLine();
    exit(1);
  }
  vector<int> intpermutation;
  for (unsigned int i = 0; i < permutation.size(); i++){
    for (unsigned int j = 0; j < permutation[i].size(); j++){
      if (!isdigit(permutation[i][j])){
        cout << "ConfigFileParser ERROR:  particle index is not an unsigned integer " << endl;
        line.printLine();
        exit(1);
      }
    }
    int ipart = atoi(permutation[i].c_str());
    if (ipart < 0 || ipart >= (int)particleList.size()){
      cout << "ConfigFileParser ERROR:  particle index is out of bounds " << endl;
      line.printLine();
      exit(1);
    }
    intpermutation.push_back(ipart);
  }
  for (unsigned int i = 0; i < intpermutation.size(); i++){
    for (unsigned int j = i+1; j < intpermutation.size(); j++){
      if (intpermutation[i] == intpermutation[j]){
        cout << "ConfigFileParser ERROR:  particle index repeated " << endl;
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
  string fixtype("floating");
  if (arguments.size() == 7) fixtype = arguments[6];
  AmplitudeInfo* amplitude = m_configurationInfo->amplitude(reaction,sumname,ampname);
  if (!amplitude){
    cout << "ConfigFileParser ERROR:  trying to initialize nonexistent amplitude " << endl;
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
    cout << "ConfigFileParser ERROR:  initializing with events is not yet re-implemented " << endl;
    line.printLine();
    exit(1);
    /*
    string normIntFile = m_configurationInfo->reaction(reaction)->normIntFile();
    if (normIntFile == ""){
      cout << "ConfigFileParser ERROR:  initializing with events, but no normalization integral " << endl;
      cout << "                         file has been specified for this reaction " << endl;
      line.printLine();
      exit(1);
    }
    ifstream test(normIntFile.c_str());
    if (!test){
      cout << "ConfigFileParser ERROR:  initializing with events, but can't find the normalization integral " << endl;
      cout << "                         file specified for this reaction (" << normIntFile << ")" << endl;
      line.printLine();
      exit(1);
    }
    NormIntInterface normint(normIntFile);
    if (!normint.hasAmpInt(ampname,ampname)){
      cout << "ConfigFileParser ERROR:  can't find the right amplitude in the normalization integral file " << endl;
      line.printLine();
      exit(1);
    }
    value1 = sqrt(value1/(abs(normint.ampInt(ampname,ampname))));
    amplitude->setValue(polar(value1,value2));
    */
  }
  else{
    cout << "ConfigFileParser ERROR:  initialize must use cartesian, polar, or events  " << endl;
    line.printLine();
    exit(1);
  }
  if (fixtype == "floating") {
  }
  else if (fixtype == "real"){
    amplitude->setReal(true);
  }
  else if (fixtype == "fixed"){
    amplitude->setFixed(true);
  }
  else{
    cout << "ConfigFileParser ERROR:  initialize must use floating, fixed, or real  " << endl;
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
    cout << "ConfigFileParser ERROR:  trying to scale nonexistent amplitude " << endl;
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
      cout << "ConfigFileParser ERROR:  can't find parameter " << parname << endl;
      line.printLine();
      exit(1);
    }
    amplitude->addParameter(parinfo);
  }
  amplitude->setScale(value);
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


