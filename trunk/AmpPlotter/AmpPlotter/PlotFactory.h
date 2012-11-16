//******************************************************************************
// This file is part of AmpPlotter, a GUI interface to AmpTools fits
// 
// Copyright Trustees of Indiana University 2012, all rights reserved
// 
// This software written by Matthew Shepherd at Indiana University, Bloomington
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
