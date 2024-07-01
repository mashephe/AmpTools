#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <mpi.h>
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpToolsMPI/AmpToolsInterfaceMPI.h"
#include "IUAmpToolsMPI/DataReaderMPI.h"
#include "IUAmpTools/FitResults.h"
#include "DalitzDataIO/DalitzDataReader.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzAmp/Constraint.h"
#include "IUAmpTools/report.h"


static const char* kModule = "fitResultsTestMPI";

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

bool testFitResults(const FitResults* fitResults) {
    string fitResultsFile = "models/fitResults.txt";
    unitTest unit_test;
    ifstream fin;
    fin.open(fitResultsFile);
    double intensity_first;
    double intensity_second;
    double pd_first;
    double pd_second;
    double ppBase_real;
    double ppBase_imag;
    double ppConstrained_real;
    double ppConstrained_imag;
    double ppSymm_real;
    double ppSymm_imag;
    double bestMinimum;
    int num_parameters;
    fin >> intensity_first;
    fin >> intensity_second;
    fin >> pd_first;
    fin >> pd_second;
    fin >> ppBase_real;
    fin >> ppBase_imag;
    fin >> ppConstrained_real;
    fin >> ppConstrained_imag;
    fin >> ppSymm_real;
    fin >> ppSymm_imag;
    fin >> bestMinimum;
    fin >> num_parameters;
    unit_test.add(abs(intensity_first-fitResults->intensity().first)<= 1, "Intensity matches model");
    unit_test.add(abs(intensity_second-fitResults->intensity().second)<= 1e1, "Intensity error matches model");
    unit_test.add(abs(pd_first-fitResults->phaseDiff("base::s1::R12","base::s1::R13").first)<= 1e-02, "Phase difference between amplitudes matches model");
    unit_test.add(abs(pd_second-fitResults->phaseDiff("base::s1::R12","base::s1::R13").second)<= 1, "Phase difference error between amplitudes matches model");
    unit_test.add(abs(ppBase_real-fitResults->productionParameter("base::s1::R12").real())<= 1e-02, "Real part of base reaction production parameter matches model");
    unit_test.add(abs(ppBase_imag-fitResults->productionParameter("base::s1::R12").imag())<= 1e-02, "Imaginary part of base reaction production parameter matches model");
    unit_test.add(abs(ppConstrained_real-fitResults->productionParameter("constrained::s2::RC12").real())<= 1e-02, "Real part of constrained reaction production parameter matches model");
    unit_test.add(abs(ppConstrained_imag-fitResults->productionParameter("constrained::s2::RC12").imag())<= 1e-02, "Imaginary part of constrained reaction production parameter matches model");
    unit_test.add(abs(ppSymm_real-fitResults->productionParameter("symmetrized_explicit::s4::RSE12").real())<= 1e-02, "Real part of symmetrized reaction production parameter matches model");
    unit_test.add(abs(ppSymm_imag-fitResults->productionParameter("symmetrized_explicit::s4::RSE12").imag())<= 1e-2, "Imaginary part of symmetrized reaction production parameter matches model");
    unit_test.add(abs(bestMinimum-fitResults->bestMinimum())<=1e-2,"Best minimum matches model");
    vector<string> parNames = fitResults->parNameList();
    int sz = parNames.size();
    unit_test.add(abs(num_parameters-sz)==0, "Number of parameter names matches model");
    vector<double> parVals = fitResults->parValueList();
    for (int i = 0; i< sz; i++) {
        double parValModel;
        fin >> parValModel;
        unit_test.add(abs(parValModel-parVals[i])<=1e-02, parNames[i] + " value matches model value");    
    }

    return unit_test.summary();
}

int main( int argc, char* argv[] ) {
    MPI_Init( &argc, &argv );
    int rank;
    int size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    string cfgname = "parserTest.cfg";
    ConfigFileParser parser(cfgname);
    ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
    AmpToolsInterfaceMPI::registerAmplitude(BreitWigner());
    AmpToolsInterfaceMPI::registerNeg2LnLikContrib(Constraint());
    AmpToolsInterfaceMPI::registerDataReader(DataReaderMPI<DalitzDataReader>());
    AmpToolsInterfaceMPI ATI(cfgInfo);
    if (rank==0) {
        vector<bool> results;
    cout << "________________________________________" << endl;
    cout << "Testing FitResults from AmpToolsInterface:" << endl;
    cout << "________________________________________" << endl;
    

    MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
    fitManager->setStrategy(1);

    fitManager->migradMinimization();
    ATI.finalizeFit();
    const FitResults* fitResults = ATI.fitResults();
    results.push_back(testFitResults(fitResults));

    cout << "________________________________________" << endl;
    cout << "Testing FitResults from file:" << endl;
    cout << "________________________________________" << endl;

    FitResults fitResults_from_file("fitTest.fit");
    const FitResults* fr_ff = &fitResults_from_file;
    results.push_back(testFitResults(fr_ff));
    for (const bool result : results) {
    if (!result) {
        throw runtime_error("Unit Tests Failed. See previous logs for more information.");
    }
  }
    }

    ATI.exitMPI();
    MPI_Finalize();

    return 0;
}