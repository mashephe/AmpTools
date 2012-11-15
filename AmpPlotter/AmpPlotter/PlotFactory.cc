
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
    
  // let's first get and stack up the acc MC components
  // get a list of the individual contributions
  THStack* accStack = new THStack( "acc_stack", m_title.c_str() );
  m_activeObjects[m_currentPad-1].push( accStack );

  const list<PlotComponent*>& accComponentList = 
    m_componentManager.accMCGroup()->componentList();
  for( list<PlotComponent*>::const_iterator part = accComponentList.begin();
       part != accComponentList.end(); ++part ){
        
    // if the part has data available and enabled ask it to deliver
    // plots
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;
         
    TH1F h = (*part)->deliverPlot( m_currentPlot );
				
    // skip over empty default histograms that may be returned
    if( h.GetEntries() == 0 ) continue;
 
    haveAcc = true;
    h.SetTitleOffset( m_showTitle ? 1 : 100 );
    accStack->Add( (TH1F*)h.Clone(), "hist" );
        
    // don't put h on the active objects stack since THStack
    // decides to take ownership of this histogram once
    // you add it        
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
        
    TH1F h = (*part)->deliverPlot( m_currentPlot );
        
    // skip over empty default histograms that may be returned
    if( h.GetEntries() == 0 ) continue;
            
    haveGen = true;
    h.SetTitleOffset( m_showTitle ? 1 : 100 );
    genStack->Add( (TH1F*)h.Clone(), "hist" );
        
    // don't put h on the active objects stack since THStack
    // decides to take ownership of this histogram once
    // you add it        
  }
    	
  // make new histo pointer for the data
  TH1F* data = 0;
	
  const list<PlotComponent*>& dataList = 
    m_componentManager.dataGroup()->componentList();
	
  for( list<PlotComponent*>::const_iterator part = dataList.begin();
       part != dataList.end(); ++part ){
		
    // if the part has data available and enabled ask it to deliver
    // plots for all of the anlaysis bins we care about
    if( ! ( (*part)->isAvailable() && (*part)->isEnabled() ) ) continue;
					
    TH1F h = (*part)->deliverPlot( m_currentPlot );
			
    // skip over empty default histograms that may be returned
    if( h.GetEntries() == 0 ) continue;

    if( ! data ){ // this is the first of many allocate memory
				
      data = (TH1F*)h.Clone();
      m_activeObjects[m_currentPad-1].push( data );
    }
    else{
				
      data->Add( &h );
    }
  }

  if( data ){
        
    data->SetStats( 0 );
    data->SetMinimum( 0 );
    data->SetTitle( m_availablePlots[m_currentPlot].c_str() );
    data->Draw( "E" );
    if( haveGen ) genStack->Draw( "SAME" );
    if( haveAcc ) accStack->Draw( "SAME" );
    data->Draw( "ESAME" );
  }
  if( !data && ( haveGen || haveAcc ) ){
			
    if( haveGen ) genStack->Draw();
    if( haveGen && haveAcc ) accStack->Draw( "SAME" );
    else if( haveAcc ) accStack->Draw(); 
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

