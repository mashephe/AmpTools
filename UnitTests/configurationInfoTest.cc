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

bool testConfigurationInfo(ConfigurationInfo* cfgInfo, string target_name) {

    cout << "________________________________________" << endl;
    cout << "Testing " << target_name << ":" << endl;
    cout << "________________________________________" << endl;
    
    unitTest unit_test;
    ifstream fin;
    fin.open("models/configurationInfo.txt");

    string model_fitname;
    fin >> model_fitname;
    unit_test.add(model_fitname == cfgInfo->fitName(), "Fit name matches model fit name");
    
    string model_fitOutputFileName;
    fin >> model_fitOutputFileName;
    unit_test.add(model_fitOutputFileName== cfgInfo->fitOutputFileName(), "Fit output file name matches model fit output file name");

    vector<string> kwds = cfgInfo->userKeywords();
    int model_kwds_size;
    fin >> model_kwds_size;
    unit_test.add(kwds.size() == model_kwds_size, "Keywords vector size matches model vector size");
    vector<ReactionInfo*> rinfo = cfgInfo->reactionList();
    int model_reaction_size;
    fin >> model_reaction_size;
    unit_test.add(rinfo.size() == model_reaction_size, "Reaction vector size matches model vector size");
    
    vector<AmplitudeInfo*> ainfo= cfgInfo->amplitudeList();
    int model_amplitude_size;
    fin >> model_amplitude_size;
    unit_test.add(ainfo.size() == model_amplitude_size, "Amplitude vector size matches model vector size");
    
    vector<CoherentSumInfo*> csinfo= cfgInfo->coherentSumList();
    int model_cs_size;
    fin >> model_cs_size;
    unit_test.add(csinfo.size() == model_cs_size, "Coherent sum vector size matches model vector size");
    
    vector<Neg2LnLikContribInfo*> llinfo= cfgInfo->neg2LnLikContribList();
    int model_ll_size;
    fin >> model_ll_size;
    unit_test.add(llinfo.size() == model_ll_size, "Neg2LnLikContrib vector size matches model vector size");
    
    vector<PDFInfo*> pdfinfo= cfgInfo->pdfList();
    int model_pdf_size;
    fin >> model_pdf_size;
    unit_test.add(pdfinfo.size() == model_pdf_size, "PDF vector size matches model vector size");
    
    vector<TermInfo*> tinfo= cfgInfo->termList();
    int model_term_size;
    fin >> model_term_size;
    unit_test.add(tinfo.size() == model_term_size, "Term vector size matches model vector size");
    
    vector<ParameterInfo*> pinfo = cfgInfo->parameterList();
    int model_parameter_size;
    fin >> model_parameter_size;
    unit_test.add(pinfo.size() == model_parameter_size, "Parameter vector size matches model vector size");
    
    bool result = unit_test.summary();
    return result;
}

int main() {
  string cfgname = "parserTest.cfg";
  vector<bool> results;
  ConfigFileParser parser(cfgname);
  ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();

  results.push_back(testConfigurationInfo(cfgInfo, "configurationInfo from ConfigFileParser.getConfigurationInfo()"));
  bool result;
  for (const bool result : results) {
    if (!result) {
        throw runtime_error("Unit Tests Failed. See previous logs for more information.");
    }
  }
  return 0;
}