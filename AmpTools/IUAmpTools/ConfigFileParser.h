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


/**
 *  The ConfigFileParser class fills a ConfigurationInfo object from either
 *  a file or a stream, depending on the constructor chosen.  It also 
 *  provides access to the raw config file information through the 
 *  getConfigFileLines method.
 *
 *  Config files follow these rules:
 *
 * #####################################
 * ####    THIS IS A CONFIG FILE    ####
 * #####################################
 * ##
 * ##  Blank lines or lines beginning with a "#" are ignored.
 * ##
 * ##  Double colons (::) are treated like a space.
 * ##  This is sometimes useful for grouping (for example,
 * ##  grouping strings like "reaction::sum::amplitudeName")
 * ##
 * ##  All non-comment lines must begin with one of the following keywords.
 * ##
 * ##  (note:  <word> means necessary 
 * ##       (word) means optional)
 * ##
 * ##  include     <file>
 * ##  define      <word> (defn1) (defn2) (defn3) ...
 * ##  fit         <fitname>
 * ##  keyword     <keyword> <min arguments> <max arguments>
 * ##  reaction     <reaction> <particle1> <particle2> (particle3) ...
 * ##  data         <reaction> <class> (arg1) (arg2) (arg3) ...
 * ##  genmc        <reaction> <class> (arg1) (arg2) (arg3) ...
 * ##  accmc        <reaction> <class> (arg1) (arg2) (arg3) ...
 * ##  normintfile  <reaction> <file> ("input")
 * ##  sum          <reaction> <sum> (sum2) (sum3) ...
 * ##  amplitude    <reaction> <sum> <amp> <class> (arg1) (arg2) ([par]) ... 
 * ##  initialize    <reaction> <sum> <amp> <"events"/"polar"/"cartesian">
 * ##       <value1> <value2> ("fixed"/"real")
 * ##  scale        <reaction> <sum> <amp> <value or [parameter]>
 * ##  constrain    <reaction1> <sum1> <amp1> <reaction2> <sum2> <amp2> ...
 * ##  permute      <reaction> <sum> <amp> <index1> <index2> ...
 * ##  parameter    <par> <value> ("fixed"/"bounded"/"gaussian") 
 * ##       (lower/central) (upper/error)
 * ##    DEPRECATED:
 * ##  datafile      <reaction> <file> (file2) (file3) ...
 * ##  genmcfile     <reaction> <file> (file2) (file3) ...
 * ##  accmcfile     <reaction> <file> (file2) (file3) ...
 * ##
 * #####################################
 * 
 */

#if !defined(CONFIGFILEPARSER)
#define CONFIGFILEPARSER

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <complex>
#include "IUAmpTools/ConfigurationInfo.h"

using namespace std;

class ConfigFileLine;


class ConfigFileParser
{

  public:

      /**
       *  Default constructor.
       */

    ConfigFileParser();


      /**
       *  A constructor that takes the name of a file as input.
       */

    ConfigFileParser(const string& configFile);


      /**
       *  A constructor that takes a stream as input.
       */

    ConfigFileParser(istream& input);


      /**
       *  A method to read configuration information from a file.  Use this 
       *  in conjunction with the default constructor as an alternative to 
       *  the ConfigFileParser(string configFile) constructor.
       */

    void readConfigFile(const string& configFile);


      /**
       *  A method to read configuration information from a stream.  Use this 
       *  in conjunction with the default constructor as an alternative to 
       *  the ConfigFileParser(istream input) constructor.
       */

    void readConfigFile(istream& input);


      /**
       *  The destructor.
       */

    ~ConfigFileParser() {}


      /**
       *  Returns a pointer to a filled ConfigurationInfo object.
       */

    ConfigurationInfo* getConfigurationInfo() {return m_configurationInfo;}


      /**
       *  Returns a vector of ConfigFileLine, which includes all parsed information
       *   (with expanded "include" and "define" statements).
       */

    const vector<ConfigFileLine>& getConfigFileLines() const {return m_configFileLines;}


      /**
       *  Displays the final parsed vector of ConfigFileLine.
       */

    void displayConfigFile() const;


      /**
       *  Control the level of output printed while parsing (useful for debugging).
       */

    static void setVerboseParsing(bool verboseParsing = false) 
                                    {m_verboseParsing = verboseParsing;}



  private:


      // read from a given file

    vector<ConfigFileLine> readConfigFileLines(const string& configFile) const;


      // read from a stream

    vector<ConfigFileLine> readConfigFileLines(istream& input) const;


      // expand the "include" and "define" statements

    vector<ConfigFileLine> expandConfigFileLines(vector<ConfigFileLine> configFileLines) const;


      // set up the ConfigurationInfo object

    void setupConfigurationInfo();


      // Do checks on the syntax, check keywords, etc.

    void checkSyntax() const;


      // Do the setup for each keyword

    void doFit           (const ConfigFileLine& line);
    void doKeyword       (const ConfigFileLine& line);
    void doReaction      (const ConfigFileLine& line);
    void doParameter     (const ConfigFileLine& line);
    void doData          (const ConfigFileLine& line);
    void doNormInt       (const ConfigFileLine& line);
    void doSum           (const ConfigFileLine& line);
    void doAmplitude     (const ConfigFileLine& line);
    void doInitialize    (const ConfigFileLine& line);
    void doPermute       (const ConfigFileLine& line);
    void doConstrain     (const ConfigFileLine& line);
    void doScale         (const ConfigFileLine& line);


      // Member data

    string                  m_fitName;
    string                  m_configFile;
    set<string>             m_userKeywords;
    static bool             m_verboseParsing;
    vector<ConfigFileLine>  m_configFileLines;
    ConfigurationInfo*      m_configurationInfo;


};


inline istream& operator>>( istream& input, ConfigFileParser& parser ){
  parser.readConfigFile( input );  return input;
}


/**
 *  The ConfigFileLine class holds a line of a parsed config file.
 *
 *    A line of the config file has the form:
 *
 *        keyword  argument1  argument2  argument3 .....
 *
 *    comment() = true for lines starting with # or empty lines
 */


class ConfigFileLine
{

  public:

    ConfigFileLine(const string& fileName, int lineNumber, const string& line);

    string          line()           const {return m_line;}
    string          fileName()       const {return m_fileName;}
    int             lineNumber()     const {return m_lineNumber;}
    string          keyword()        const {return m_keyword;}
    vector<string>  arguments()      const {return m_arguments;}
    bool            comment()        const {return m_comment;}

    void            printLine()      const {cout << "(" << m_lineNumber << ") " << 
                                                    m_fileName   << " >>   " << 
                                                    m_line       << endl;}

    void            printArguments() const {cout << "KEYWORD:  " << m_keyword << endl;
                                            cout << "ARGUMENTS: " << endl;
                                            for ( unsigned int i = 0; i < m_arguments.size(); i++){
                                              cout << m_arguments[i] << endl;}}

    void            flushDefinition(const string& word, const vector<string>& definition);

  private:

    string          m_line;
    string          m_fileName;
    int             m_lineNumber;
    string          m_keyword;
    vector<string>  m_arguments;
    bool            m_comment;

};


#endif
