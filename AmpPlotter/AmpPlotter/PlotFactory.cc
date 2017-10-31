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

#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TH1.h"

#include "AmpPlotter/PlotFactory.h"
#include "IUAmpTools/PlotGenerator.h"
#include "IUAmpTools/ConfigFileParser.h"

using namespace std;

PlotFactory::PlotFactory( PlotGenerator& plotGen,
                          bool boringPlots ) :
  m_plotGen( plotGen ),
  m_componentManager( plotGen.reactions(), plotGen, boringPlots ),
  m_canvas(  ),
  m_title( "" ),
  m_showTitle( false ),
  m_availablePlots( plotGen.availablePlots() ),
  m_currentPlot( 0 )
{
  // preserve your sanity -- trust yourself rather than ROOT
  // to handle your pointers properly
  TH1::AddDirectory( kFALSE ); 

  m_canvas = new TCanvas( "canvas", "Amplitude Projection Plotter", 500, 500 );   
  m_canvas->Divide(1,1);
  m_currentPad = 1;
  m_currentNoPads = 1;
  m_MaxNoPads = 9;   // if you change m_MaxNoPads remember to change size of m_activeObjects.
}

void
PlotFactory::setPlot( unsigned int plotIndex ){

  m_currentPlot = plotIndex;
}

void
PlotFactory::setPad( unsigned int padNo ){

  m_currentPad = padNo;
  if( m_currentPad > m_currentNoPads )m_currentPad = m_currentNoPads;

}

void
PlotFactory::chooseCanvLayout( unsigned int canvOpt ){

  clearCanvas(); // want to clear out cache of hists for each pad.
  m_canvas->Clear();
  switch( canvOpt ){
  case 1:
    m_canvas->Divide(1,1);
    m_currentNoPads = 1;
    break;
  case 2:
    m_canvas->Divide(2,1);
    m_currentNoPads = 2;
    break;
  case 3:
    m_canvas->Divide(2,2);
    m_currentNoPads = 4;
    break;
  case 4:
    m_canvas->Divide(3,2);
    m_currentNoPads = 6;
    break;
  case 5:
    m_canvas->Divide(2,3);
    m_currentNoPads = 6;
    break;
  case 6:
    m_canvas->Divide(3,3);
    m_currentNoPads = 9;
    break;
  case 7:
    m_canvas->Divide(2,4);
    m_currentNoPads = 8;
    break;
  case 8:
    m_canvas->Divide(4,2);
    m_currentNoPads = 8;
    break;
  }
  m_currentPad = 1;

}

void
PlotFactory::drawPlot( void )
{	
  string plotName = m_availablePlots[m_currentPlot];
	
  // flush out the objects from the last plot
  m_canvas->GetPad(m_currentPad)->Clear();
  m_canvas->cd(m_currentPad);
  cleanObjectCache();
	
  bool haveGen = false;
  bool haveAcc = false;
  bool is2D = false;
    
  // let's first get and stack up the acc MC components
  // get a list of the individual contributions --
  THStack* accStack = new THStack( "acc_stack", m_title.c_str() );
  m_activeObjects[m_currentPad-1].push( accStack );

  // background should be on the bottom of the acc_stack only
  const list<PlotComponent*>& bkgndComponentList =
  m_componentManager.bkgndGroup()->componentList();
  for( list<PlotComponent*>::const_iterator part = bkgndComponentList.begin();
      part != bkgndComponentList.end(); ++part ){
    
    // if the part has data available and enabled ask it to deliver
    // plots
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;

    TH1* h = (*part)->deliverPlot( m_currentPlot );
				
    // skip over empty default histograms that may be returned
    if( h->GetEntries() == 0 ) continue;
    
    haveAcc = true;
    
    is2D = ( h->GetDimension() == 2 );
    h->SetTitleOffset( m_showTitle ? 1 : 100 );
    accStack->Add( h, "hist" );
  }
  
  // now stack on the accMC components
  const list<PlotComponent*>& accComponentList =
    m_componentManager.accMCGroup()->componentList();
  for( list<PlotComponent*>::const_iterator part = accComponentList.begin();
       part != accComponentList.end(); ++part ){
        
    // if the part has data available and enabled ask it to deliver
    // plots
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;
    
    TH1* h = (*part)->deliverPlot( m_currentPlot );
				
    // skip over empty default histograms that may be returned
    if( h->GetEntries() == 0 ) continue;

    haveAcc = true;

    is2D = ( h->GetDimension() == 2 );
    h->SetTitleOffset( m_showTitle ? 1 : 100 );
    accStack->Add( h, "hist" );
  }

  // now the generated components
  THStack* genStack = new THStack( "gen_stack", m_title.c_str() );
  m_activeObjects[m_currentPad-1].push( genStack );
    
  const list<PlotComponent*>& genComponentList = 
    m_componentManager.genMCGroup()->componentList();
  for( list<PlotComponent*>::const_iterator part = genComponentList.begin();
       part != genComponentList.end(); ++part ){
        
    // if the part has data available and enabled ask it to deliver
    // plots
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;
        
    TH1* h = (*part)->deliverPlot( m_currentPlot );
        
    // skip over empty default histograms that may be returned
    if( h->GetEntries() == 0 ) continue;
    
    is2D = ( h->GetDimension() == 2 );
    haveGen = true;
    h->SetTitleOffset( m_showTitle ? 1 : 100 );
    genStack->Add( h, "hist" );
  }
    	
  // make new histo pointer for the data
  TH1* data = 0;
	
  const list<PlotComponent*>& dataList = 
    m_componentManager.dataGroup()->componentList();
	
  for( list<PlotComponent*>::const_iterator part = dataList.begin();
       part != dataList.end(); ++part ){
		
    // if the part has data available and enabled ask it to deliver
    // plots for all of the anlaysis bins we care about
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;
					
    TH1* h = (*part)->deliverPlot( m_currentPlot );
			
    // skip over empty default histograms that may be returned
    if( h->GetEntries() == 0 ) continue;

    is2D = ( h->GetDimension() == 2 );

    if( ! data ){ // this is the first of many allocate memory
				
      data = h;
      m_activeObjects[m_currentPad-1].push( data );
    }
    else{
				
      data->Add( h );
    }
  }

  if( data ){
        
    data->SetStats( 0 );
    data->SetMinimum( 0 );
    data->SetTitle( m_availablePlots[m_currentPlot].c_str() );
    if( is2D ){
      data->Draw( "SCAT" );
      if( haveAcc ) accStack->Draw( "CONTZSAME" );
      if( haveGen ) genStack->Draw( "CONT3SAME" );
      data->Draw( "SCATSAME" );
    }
    else{
      data->Draw( "E" );
      if( haveGen ) genStack->Draw( "SAME" );
      if( haveAcc ) accStack->Draw( "SAME" );
      data->Draw( "ESAME" );
    }
  }
  if( !data && ( haveGen || haveAcc ) ){
			
    if( is2D ){
      
      if( haveGen ) genStack->Draw( "CONT3" );
      if( haveGen && haveAcc ) accStack->Draw( "CONTZSAME" );
      else if( haveAcc ) accStack->Draw( "CONTZ" );
    }
    else{
      
      if( haveGen ) genStack->Draw();
      if( haveGen && haveAcc ) accStack->Draw( "SAME" );
      else if( haveAcc ) accStack->Draw();
    }
  }

  m_canvas->cd();
  m_canvas->Update();
}

void
PlotFactory::printPlot( const char* fileName ){

  drawPlot();
  m_canvas->Print( fileName );
}

void
PlotFactory::cleanObjectCache( void ){
  
  while( ! m_activeObjects[m_currentPad-1].empty() ){
    m_activeObjects[m_currentPad-1].top()->Delete();
    m_activeObjects[m_currentPad-1].pop();
  }

}

void
PlotFactory::clearCanvas( void ){
  // Clear the canvas but make sure you clean out the
  // the stack of objects for each pad as well.

  for( int i=0; i<m_MaxNoPads; i++ ){
    setPad(i+1);
    cleanObjectCache();
    if( i <= m_currentNoPads ){
      m_canvas->GetPad(m_currentPad)->Clear();
      m_canvas->Update();
    }
  }
}

