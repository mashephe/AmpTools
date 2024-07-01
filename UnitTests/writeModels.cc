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
#include "IUAmpTools/FitResults.h"
#include "DalitzAmp/BreitWigner.h"
#include "DalitzAmp/Constraint.h"

#include "IUAmpTools/report.h"
static const char* kModule = "writeModels";
using namespace std;

int main() {
    ofstream fout;
    string cfgname("parserTest.cfg");
    ConfigFileParser parser(cfgname);
    ConfigurationInfo* cfgInfo = parser.getConfigurationInfo();
    AmpToolsInterface::registerAmplitude(BreitWigner());
    AmpToolsInterface::registerNeg2LnLikContrib(Constraint());
    AmpToolsInterface::registerDataReader(DalitzDataReader());
    AmpToolsInterface ATI(cfgInfo);


    //ConfigFileParser
    vector<ConfigFileLine> cfgLines = parser.getConfigFileLines();
    fout.open("models/parsedConfig.txt");
    vector<string> cfgStrings;
    for (const ConfigFileLine cfgLine : cfgLines) {
        fout << cfgLine.line() << "\n";
    }
    fout.close();

    //ConfigurationInfo
    fout.open("models/configurationInfo.txt");
    fout << cfgInfo->fitName() << "\n";
    fout << cfgInfo->fitOutputFileName() << "\n";
    vector<string> kwds = cfgInfo->userKeywords();
    fout << kwds.size()<<"\n";
    vector<ReactionInfo*> rinfo = cfgInfo->reactionList();
    fout << rinfo.size() << "\n";
    vector<AmplitudeInfo*> ainfo= cfgInfo->amplitudeList();
    fout << ainfo.size() << "\n";
    vector<CoherentSumInfo*> csinfo= cfgInfo->coherentSumList();
    fout << csinfo.size() << "\n";
    vector<Neg2LnLikContribInfo*> llinfo= cfgInfo->neg2LnLikContribList();
    fout << llinfo.size() << "\n";
    vector<PDFInfo*> pdfinfo= cfgInfo->pdfList();
    fout << pdfinfo.size() << "\n";
    vector<TermInfo*> tinfo= cfgInfo->termList();
    fout << tinfo.size() << "\n";
    vector<ParameterInfo*> pinfo = cfgInfo->parameterList();
    fout << pinfo.size() << "\n"; 

    fout.close();

    //AmpToolsInterface
    string ATIFile = "models/AmpToolsInterface.txt";
    fout.open(ATIFile);
    double neg2LL_before = ATI.likelihood();
    fout << neg2LL_before << "\n";

    MinuitMinimizationManager* fitManager = ATI.minuitMinimizationManager();
    fitManager->setStrategy(1);

    fitManager->migradMinimization();


    double neg2LL_after = ATI.likelihood();
    fout << neg2LL_after << "\n";
    fout << ATI.likelihood("base") << "\n";
    fout << ATI.likelihood("constrained") << "\n";
    fout << ATI.likelihood("symmetrized_implicit") << "\n";
    fout << ATI.likelihood("symmetrized_explicit") << "\n";
    ATI.finalizeFit();
    fout.close();

    //fitResults

    const FitResults* fitResults = ATI.fitResults();
    fout.open("models/fitResults.txt");
    pair<double, double> intensity = fitResults->intensity();
    fout << intensity.first << "\n";
    fout << intensity.second << "\n";
    pair<double, double> pd = fitResults->phaseDiff("base::s1::R12","base::s1::R13");
    fout << pd.first << "\n";
    fout << pd.second << "\n";
    complex<double> ppBase = fitResults->productionParameter("base::s1::R12");
    fout << ppBase.real() << "\n";
    fout << ppBase.imag() << "\n";
    complex<double> ppConstrained = fitResults->productionParameter("constrained::s2::RC12");
    fout << ppConstrained.real() << "\n";
    fout << ppConstrained.imag() << "\n";
    complex<double> ppSymm = fitResults->productionParameter("symmetrized_explicit::s4::RSE12");
    fout << ppSymm.real() << "\n";
    fout << ppSymm.imag() << "\n";
    double bestMinimum = fitResults->bestMinimum();
    fout << bestMinimum << "\n";
    vector<string> parNames = fitResults->parNameList();
    fout << parNames.size() << "\n";
    vector<double> parVals = fitResults->parValueList();
    for (const double i : parVals) {
        fout << i << "\n";
    }
    fout.close();
    return 0;
}