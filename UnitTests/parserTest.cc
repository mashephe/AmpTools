#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <stdio.h>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/AmpToolsInterface.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzAmp/Constraint.h"

using namespace std;

class unitTest {
    public:
    bool passed = true;
    vector<string> failedTests;
    vector<string> passedTests;
    void add(bool expr, string name) {
        passed = passed && expr;
        if (expr) {
            passedTests.push_back(name);
        } else {
            failedTests.push_back(name);
        }
    }
    bool summary() {
        assert(passed || failedTests.size() > 0);
        if (passed) {
            cout << "All unit tests passed." << endl;
        } else {
            cout << "The following unit tests failed:" << endl;
            for (const string& failedTest : failedTests) {
                cout << "* " << failedTest << endl;
            }
            if (passedTests.size() > 0) {
            cout << "The following unit tests were successful:" << endl;
            for (const string& passedTest : passedTests) {
                cout << "* " << passedTest << endl;
            }
            }
        }
        return passed;
    }
};

vector<string> read_model(string model_file_name) {
    ifstream fin;
    fin.open(model_file_name);
    vector<string> lines;
    string line;
    while (getline(fin, line)) {
        lines.push_back(line);
    }
    return lines;
}

bool testConfigFileParser(ConfigFileParser parser, string target_name) {
    cout << "________________________________________" << endl;
    cout << "Testing " << target_name << ":" << endl;
    cout << "________________________________________" << endl;
    
    unitTest unit_test;

    vector<ConfigFileLine> cfgLines = parser.getConfigFileLines();
 
    //checking that the file was accurately parsed
    
    vector<string> cfgStrings;
    string cfgLine;
    for (const ConfigFileLine cfgLine : cfgLines) {
        cfgStrings.push_back(cfgLine.line());
    }


    vector<string> modelStrings = read_model("models/parsedConfig.txt");
    unit_test.add(modelStrings.size() == cfgStrings.size(), "Number of lines in parsed file matches model");
    unit_test.add(modelStrings == cfgStrings, "Parsed file exactly matches model");
    bool result = unit_test.summary();
    cout << endl;

    return result;
}

int main() {
  string cfgname = "parserTest.cfg";
  vector<bool> results;
  
  // Testing ConfigFileParser(const string& configFile)
  ConfigFileParser parser_from_string(cfgname);
  results.push_back(testConfigFileParser(parser_from_string, "ConfigFileParser(const string& configFile)"));

  // Testing ConfigFileParser(istream& input)
  ifstream fin;
  fin.open(cfgname);
  ConfigFileParser parser_from_istream(fin);
  results.push_back(testConfigFileParser(parser_from_istream, "ConfigFileParser(istream& input)"));
  fin.close();

  // Testing ConfigFileParser() and readConfigFile(const string& configFile)
  ConfigFileParser default_with_read_string;
  default_with_read_string.readConfigFile(cfgname);
  results.push_back(testConfigFileParser(default_with_read_string, "ConfigFileParser() and readConfigFile(const string& configFile)"));

  // Testing ConfigFileParser() and readConfigFile(istream& input)
  fin.open(cfgname);
  ConfigFileParser default_with_read_istream;
  default_with_read_istream.readConfigFile(fin);
  results.push_back(testConfigFileParser(default_with_read_istream, "ConfigFileParser() and readConfigFile(istream& input)"));
  fin.close();

  bool result;
  for (const bool result : results) {
    if (!result) {
        throw runtime_error("Unit Tests Failed. See previous logs for more information.");
    }
  }
  return 0;
}