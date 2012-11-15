#if !(defined PLOTFACTORY)
#define PLOTFACTORY

#include <stack>
#include <list>
#include <string>

#include "TCanvas.h"
#include "TH1.h"

#include "AmpPlotter/Plot.h"
#include "AmpPlotter/PlotComponentManager.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/ConfigFileParser.h"

using namespace std;

class PlotFactory
{
	
 public:
	
  PlotFactory( PlotGenerator& plotGen,
	       bool boringPlots = false );
    
  ~PlotFactory(){};
	
  PlotComponentManager& componentManager() { return m_componentManager; }
  PlotGenerator& generator() { return m_plotGen; }

  void setPlot( unsigned int plotIndex );
  void setPad( unsigned int padNo );
  void clearCanvas( void );
  void chooseCanvLayout( unsigned int canvOpt );

  void drawPlot( void );
  void printPlot( const char* fileName );
  void printPlot( string fileName ){ printPlot( fileName.c_str() ); }
	
  void setTitle( const string& title ) { m_title = title; }
  void setShowTitle( bool showTitle ) { m_showTitle = showTitle; }

  TCanvas* getCanvas(){ return m_canvas; }

  unsigned int getMaxNoPads(){ return m_MaxNoPads; }
		
  const vector< string >& availablePlots() const { return m_availablePlots; }
	
 private:
		
  map< string, map< string, int > >
    buildReactionAmpMap( const vector< string >& indexFiles, bool havePars ) const;
    
    
  void cleanObjectCache( void );
		
  PlotGenerator& m_plotGen;
  PlotComponentManager m_componentManager;
	
  TCanvas* m_canvas;
  unsigned int m_currentPad;
  unsigned int m_MaxNoPads;
  unsigned int m_currentNoPads;
	
  string m_title;
  bool m_showTitle;
 
  vector< string > m_availablePlots;
  unsigned int m_currentPlot;
	
  //stack<TObject*> m_activeObjects;
  // edited this - not very elegant - an array of stacks - 
  // one for each pad on the canvas. Hardwired to 9 
  // i.e. the maximum number of pads.
  stack<TObject*> m_activeObjects[9];
    
  bool m_havePars;
};

#endif
