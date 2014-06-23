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

#if !(defined PLOTCOMPONENT)
#define PLOTCOMPONENT

#include <string>
#include <list>

#include "IUAmpTools/PlotGenerator.h"

class TH1F;
class AnalysisBin;
class ComponentGroup;

using namespace std;

class PlotComponent
{
	
public:
		
	PlotComponent( const string& title, 
                   const string& reaction,
                   PlotGenerator& pltGen ) : 
  m_pltGen( pltGen ),
  m_title( title ),
  m_reaction( reaction ),
  m_dataOK( false ),
  m_fillColor(   10 ),
  m_fillStyle(    0 ),
  m_markerStyle( 20 ),
  m_markerSize( 0.5 ){}
	
	virtual ~PlotComponent();

	virtual TH1* deliverPlot( unsigned int plotIndex,
                            double scale = 1.0 ) = 0;
    
	bool isAvailable() const;
	bool isEnabled() const;
	
	int fillStyle() const { return m_fillStyle; }
	int fillColor() const { return m_fillColor; }
	
	double markerSize() const { return m_markerSize; }
	int    markerStyle() const { return m_markerStyle; }
	
	string title() const { return m_title; }
  string reaction() const { return m_reaction; }
  
  virtual string amp() const { return ""; }
  
  PlotGenerator& generator() { return m_pltGen; }
  
	void setDataOK( void ) { m_dataOK = true; }
	
	void setFillStyle( int style ) { m_fillStyle = style; }
	void setFillColor( int color ) { m_fillColor = color; }
	
	void setMarkerStyle( int style ) { m_markerStyle = style; }
	void setMarkerSize( double size ) { m_markerSize = size; }
	
	void addGroup( ComponentGroup* group );
	void deleteGroup( ComponentGroup* group );
	
private:
		
	PlotComponent( const PlotComponent& );
	const PlotComponent& operator=( const PlotComponent& );
			
    PlotGenerator& m_pltGen;
    
	string  m_title;
    string       m_reaction;
	bool         m_dataOK;
	
	list<ComponentGroup*> m_groups;
	
	int m_fillColor;
	int m_fillStyle;
	
	int    m_markerStyle;
	double m_markerSize;
};

#endif
