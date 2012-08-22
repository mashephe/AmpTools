#if !(defined GAMMAKKPLOTGENERATOR)
#define GAMMAKKPLOTGENERATOR

#include <vector>
#include <string>
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "gammaKKDataIO/gammaKKDataReader.h"


using namespace std;

class gammaKKPlotGenerator : public PlotGenerator
{
    
public:
    
  gammaKKPlotGenerator( const ConfigurationInfo* cfgInfo,
                       const string& parFile );

  virtual ~gammaKKPlotGenerator();

  const vector< string >& availablePlots() const;

  enum HistIndex { 
    khm12 = 0, khm13, khm23,
    // angles of photon in CM frame
    khPhotonCosTheta, khPhotonPhi,
    // angles of 1st K in KK CM frame
    khKaonCosTheta, khKaonPhi,
    kNumHist 
  };
    
private:
        
  vector< Histogram > fillProjections( const string& fsName,
                                          PlotType type );

  void registerPhysics( AmplitudeManager* ampManager );

  vector< string > m_histTitles;

  mutable vector< Histogram > m_histVect;

  void clearHistograms();
  void bookHistograms();

  map< string, AmpVecs > m_dataMap;
  map< string, AmpVecs > m_accMCMap;
  map< string, AmpVecs > m_genMCMap;

};

#endif
