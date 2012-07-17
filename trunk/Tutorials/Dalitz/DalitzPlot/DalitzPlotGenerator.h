#if !(defined DALITZPLOTGENERATOR)
#define DALITZPLOTGENERATOR

#include <vector>
#include <string>
#include "IUAmpTools/ConfigurationInfo.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/Histogram.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "DalitzDataIO/DalitzDataReader.h"


using namespace std;

class DalitzPlotGenerator : public PlotGenerator
{
    
public:
    
  DalitzPlotGenerator( const ConfigurationInfo* cfgInfo,
                       const string& parFile );

  virtual ~DalitzPlotGenerator();

  const vector< string >& availablePlots() const;

  enum HistIndex { 
    khm12 = 0, khm13, khm23, 
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
