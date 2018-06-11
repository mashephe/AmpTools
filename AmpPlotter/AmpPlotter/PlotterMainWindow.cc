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

#include "TROOT.h"
#include "TGWindow.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGLabel.h"
#include "TGTab.h"
#include "TGClient.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/ComponentGroup.h"
#include "AmpPlotter/PlotFactory.h"
#include "AmpPlotter/PlotComponentManager.h"

using namespace std;

PlotterMainWindow::PlotterMainWindow( const TGWindow* win, 
                                     PlotFactory& factory ) :
TGMainFrame( win, kWidth, kHeight ),
m_factory( factory ),
m_componentManager( factory.componentManager() ),
m_generator( factory.generator() )
{
  
  
  m_mainFrame = new TGVerticalFrame( this, kWidth, kHeight, kChildFrame );
	
  // divide main frame  into menubar and "the rest"
  m_menuFrame = new TGHorizontalFrame( m_mainFrame, kWidth, 2*kHeight/15, kChildFrame );
  m_restFrame = new TGHorizontalFrame( m_mainFrame, kWidth, 13*kHeight/15, kChildFrame );
	
  // setup the right and left frames
  m_rightFrame = new TGVerticalFrame( m_restFrame, kWidth/2, 13*kHeight/15, kChildFrame );
  m_leftFrame = new TGVerticalFrame( m_restFrame, kWidth/2, 13*kHeight/15, kChildFrame );
  
  TGLayoutHints mainLayoutHints( kLHintsCenterY | kLHintsExpandX, 5, 5, 5, 5 );
  
  // the top group of checkboxes for the decay reaction
  m_reactionFrame = new TGGroupFrame( m_rightFrame, "Reaction", kVerticalFrame );
  TGLayoutHints reactionLayoutHints( kLHintsBottom | kLHintsExpandX, 0, 0, 5, 0 );
  
  const std::vector<ComponentGroup*>& reactions = m_componentManager.reactionGroups();
  for( std::vector<ComponentGroup*>::const_iterator reaction = 
      reactions.begin(); reaction != reactions.end(); ++reaction ){
		
    m_reactionButtons.
    push_back( new TGCheckButton( m_reactionFrame, (*reaction)->title().c_str(),
                                 ( reaction - reactions.begin() ) | kReaction ) );
		
    m_reactionButtons.back()->Associate( this );
    m_reactionFrame->AddFrame( m_reactionButtons.back(), &reactionLayoutHints );
  }
	
  m_typeFrame = new TGGroupFrame( m_rightFrame, "Plots to Include", kVerticalFrame );
  m_typeButtons.resize( 3 );
  
  m_typeButtons[kGenMC] = new TGCheckButton( m_typeFrame, "Generated Monte Carlo", kGenMC | kType );
  m_typeButtons[kGenMC]->Associate( this );
  m_typeFrame->AddFrame( m_typeButtons[kGenMC], &reactionLayoutHints );
  
  m_typeButtons[kAccMC] = new TGCheckButton( m_typeFrame, "Accepted Monte Carlo", kAccMC | kType );
  m_typeButtons[kAccMC]->Associate( this );
  m_typeFrame->AddFrame( m_typeButtons[kAccMC], &reactionLayoutHints );

  m_typeButtons[kBkgnd] = new TGCheckButton( m_typeFrame, "Background (if used)", kBkgnd | kType );
  m_typeButtons[kBkgnd]->Associate( this );
  m_typeFrame->AddFrame( m_typeButtons[kBkgnd], &reactionLayoutHints );
  
  m_typeButtons[kData] = new TGCheckButton( m_typeFrame, "Data", kData | kType );
  m_typeButtons[kData]->Associate( this );
  m_typeFrame->AddFrame( m_typeButtons[kData], &reactionLayoutHints );
  
  m_ampSumTab = new TGTab( m_leftFrame );
  m_ampSumTab->AddTab( "Amplitudes" );
  m_ampSumTab->AddTab( "Coherent Sums" );
  
  // layout the frame that holds the amplitude selectors
  
  m_ampFrame = m_ampSumTab->GetTabContainer( "Amplitudes" );
  
  TGLayoutHints buttonLayoutHints( kLHintsExpandX, 2, 2, 2, 2 );
  
  m_selClrAmpFrame = new TGHorizontalFrame( m_ampFrame, kWidth/2, 1, kFitHeight );
  m_selectAllAmpButton = new TGTextButton( m_selClrAmpFrame, "&Select All", kSelectAllAmp );
  m_selectAllAmpButton->Associate( this );
  m_selClrAmpFrame->AddFrame( m_selectAllAmpButton, &buttonLayoutHints );
  
  m_clearAmpButton = new TGTextButton( m_selClrAmpFrame, "&Clear", kClearAmp );
  m_clearAmpButton->Associate( this );
  m_selClrAmpFrame->AddFrame( m_clearAmpButton, &buttonLayoutHints );
  m_ampFrame->AddFrame( m_selClrAmpFrame, &reactionLayoutHints );
  
  m_ampSelectBox = new TGListBox( m_ampFrame, kAmpSelect );
  m_ampSelectBox->SetMultipleSelections( true );
  const std::vector<string>& amps = m_generator.uniqueAmplitudes();
  for( std::vector<string>::const_iterator amp = 
      amps.begin(); amp != amps.end(); ++amp ){
    
    int i = amp - amps.begin();
    
    //    cout << "adding:  " << *amp << endl;
    m_ampSelectBox->AddEntry( amp->c_str(), i );
    // default to on
    m_ampSelectBox->Select( i );
    m_generator.enableAmp( i );
  }
  m_ampSelectBox->Resize( 4*kWidth/10, 10*kHeight/15 );
  m_ampSelectBox->MapSubwindows();
  m_ampSelectBox->Layout();
  m_ampSelectBox->Associate( this );
  m_ampFrame->AddFrame( m_ampSelectBox, &reactionLayoutHints );
  
  // layout the frame that holds the sum selectors
  
  m_sumFrame = m_ampSumTab->GetTabContainer( "Coherent Sums" );
    
  m_selClrSumFrame = new TGHorizontalFrame( m_sumFrame, kWidth/2, 1, kFitHeight );
  m_selectAllSumButton = new TGTextButton( m_selClrSumFrame, "&Select All", kSelectAllSum );
  m_selectAllSumButton->Associate( this );
  m_selClrSumFrame->AddFrame( m_selectAllSumButton, &buttonLayoutHints );
  
  m_clearSumButton = new TGTextButton( m_selClrSumFrame, "&Clear", kClearSum );
  m_clearSumButton->Associate( this );
  m_selClrSumFrame->AddFrame( m_clearSumButton, &buttonLayoutHints );
  m_sumFrame->AddFrame( m_selClrSumFrame, &reactionLayoutHints );
  
  m_sumSelectBox = new TGListBox( m_sumFrame, kSumSelect );
  m_sumSelectBox->SetMultipleSelections( true );
  const std::vector<string>& sums = m_generator.uniqueSums();
  for( std::vector<string>::const_iterator sum = 
      sums.begin(); sum != sums.end(); ++sum ){
    
    int i = sum - sums.begin();
    
    //    cout << "adding:  " << *sum << endl;
    m_sumSelectBox->AddEntry( sum->c_str(), i );
    // default to on
    m_sumSelectBox->Select( i );
    m_generator.enableSum( i );
  }
  m_sumSelectBox->Resize( 4*kWidth/10, 10*kHeight/15 );
  m_sumSelectBox->MapSubwindows();
  m_sumSelectBox->Layout();
  m_sumSelectBox->Associate( this );
  m_sumFrame->AddFrame( m_sumSelectBox, &reactionLayoutHints );
  
  
  // create a frame to hold the plot select box and show data check button
  m_plotSelectFrame = new TGHorizontalFrame( m_menuFrame, kWidth, 1, kFitHeight );
  TGLayoutHints plotSelectLHints( kLHintsCenterX, 2, 2, 5, 5 );
  m_plotSelectBox = new TGComboBox( m_plotSelectFrame, kChoosePlot );
  const std::vector< string >& plots = m_factory.availablePlots();
  for( std::vector< string >::const_iterator aPlot = plots.begin();
      aPlot != plots.end(); ++aPlot ){
    
    m_plotSelectBox->AddEntry( aPlot->c_str(), 
                              aPlot - plots.begin() );
  }
  m_plotSelectBox->Resize( kWidth/3, 20 );
  m_plotSelectBox->Select( 0 );
  m_plotSelectBox->Associate( this );
  m_plotSelectFrame->AddFrame( m_plotSelectBox, &plotSelectLHints );
	
  
  // now add a plot and exit button on the bottom
  m_buttonFrame = new TGHorizontalFrame( m_menuFrame, kWidth,1, kFitHeight );
  m_plotButton = new TGTextButton( m_buttonFrame, "&Plot", kPlot );
  m_plotButton->Associate( this );
  m_buttonFrame->AddFrame( m_plotButton, &buttonLayoutHints );
	
  m_exitButton = new TGTextButton( m_buttonFrame, "&Exit", kExit );
  m_exitButton->Associate( this );
  //	m_exitButton->SetCommand( ".q" );
  m_buttonFrame->AddFrame( m_exitButton, &buttonLayoutHints );
  
  
  // create a new frame to hold the canvas layout options
  m_canvFrame = new TGGroupFrame( m_rightFrame, "Draw Options", kVerticalFrame );
  TGLayoutHints canvLayoutHints( kLHintsBottom | kLHintsExpandX, 0, 0, 5, 0 );
  
  // create a list box to select which pad to draw on
  m_padButton = new TGComboBox( m_canvFrame, kChoosePad );
  char tmp[22];
  for( int i=0; i<factory.getMaxNoPads(); i++ ){
    sprintf(tmp,"Select pad %d",i+1);
    m_padButton->AddEntry( tmp, i+1 );
  }
  m_padButton->Resize( kWidth/3, 20 );
  m_padButton->Select( 1 );
  m_padButton->Associate( this );
  m_canvFrame->AddFrame( m_padButton, &canvLayoutHints );
  
  // create a list box to select the canvas layout...
  m_canvButton = new TGComboBox( m_canvFrame, kChooseCanv );
  m_canvButton->AddEntry( "1x1 pads", 1 );
  m_canvButton->AddEntry( "2x1 pads", 2 );
  m_canvButton->AddEntry( "2x2 pads", 3 );
  m_canvButton->AddEntry( "3x2 pads", 4 );
  m_canvButton->AddEntry( "2x3 pads", 5 );
  m_canvButton->AddEntry( "3x3 pads", 6 );
  m_canvButton->AddEntry( "2x4 pads", 7 );
  m_canvButton->AddEntry( "4x2 pads", 8 );
  m_canvButton->Resize( kWidth/3, 20 );
  m_canvButton->Select( 1 );
  m_canvButton->Associate( this );
  m_canvFrame->AddFrame( m_canvButton, &canvLayoutHints );
  
  // create a button to clear the canvas...
  m_clearCanvButton = new TGTextButton( m_canvFrame, "&Clear Canvas", kclearCanv );
  m_clearCanvButton->Associate( this );
  m_clearCanvButton->Resize( kWidth/3, 20 );
  m_canvFrame->AddFrame( m_clearCanvButton, &canvLayoutHints );
  
  
  // add the sub-frames to the main frame
	
  // the menu bar
  m_menuFrame->AddFrame( m_plotSelectFrame, &mainLayoutHints );
  m_menuFrame->AddFrame( m_buttonFrame, &mainLayoutHints );
  m_menuFrame->Layout();
  
  // build the right and left frames
  m_rightFrame->AddFrame( m_reactionFrame, &mainLayoutHints ); 
  m_rightFrame->AddFrame( m_typeFrame, &mainLayoutHints );
  m_rightFrame->AddFrame( m_canvFrame, &mainLayoutHints );  
  m_leftFrame->AddFrame( m_ampSumTab, &mainLayoutHints );
	
  m_restFrame->AddFrame( m_leftFrame, &mainLayoutHints );
  m_restFrame->AddFrame( m_rightFrame, &mainLayoutHints );
  
  m_mainFrame->AddFrame( m_menuFrame, &mainLayoutHints );
  m_mainFrame->AddFrame( m_restFrame, &mainLayoutHints );
  
  m_mainFrame->Resize( m_mainFrame->GetDefaultSize() );
  m_mainFrame->MapSubwindows();
  m_mainFrame->Layout();
  m_mainFrame->MapWindow();
	
  // and stick the main frame in the main window
  TGLayoutHints mainWindowHints( kLHintsExpandX | kLHintsCenterY );
  AddFrame( m_mainFrame, &mainWindowHints );
	
  MapSubwindows();
  Layout();
	
  SetWindowName( "Amplitude Projection Plotter" );
	
  MapWindow();
}

bool
PlotterMainWindow::ProcessMessage( long mes, long p1, long p2 )
{
  
  int index;
	
  switch( GET_MSG( mes ) ){
      
    case kC_COMMAND:
			
      switch( GET_SUBMSG( mes ) ){
          
        case kCM_CHECKBUTTON:
   				
          index = kIndexMask & p1;
          switch( kButtonMask & p1 ){
              
            case kReaction:
							
              if( m_reactionButtons[index]->GetState() == kButtonDown ){
								
                m_componentManager.enableReaction( index );
              }
              else{
								
                m_componentManager.disableReaction( index );
              }
              break;
							
            case kType:
              
              switch( index ) {
                  
                case kData:
                  
                  if( m_typeButtons[kData]->GetState()  == kButtonDown ){
                    
                    m_componentManager.enableData();
                  }
                  else{
                    
                    m_componentManager.disableData();
                  }
                  break;
                case kBkgnd:
                                    
                  if( m_typeButtons[kBkgnd]->GetState()  == kButtonDown ){
                    
                    m_componentManager.enableBkgnd();
                  }
                  else{
                    
                    m_componentManager.disableBkgnd();
                  }
                  break;
                case kAccMC:
                  
                  if( m_typeButtons[kAccMC]->GetState()  == kButtonDown ){
                    
                    m_componentManager.enableAccMC();
                  }
                  else{
                    
                    m_componentManager.disableAccMC();
                  }
                  break;
                case kGenMC:
                  
                  if( m_typeButtons[kGenMC]->GetState()  == kButtonDown ){
                    
                    m_componentManager.enableGenMC();
                  }
                  else{
                    
                    m_componentManager.disableGenMC();
                  }
                  break;
              }
              break;
              
            default:
							
              cout << "Unkown button click!" << endl;
              break;
          } // end of check buttons
          break;
          
          
        case kCM_LISTBOX:
          
          switch( p1 ){
              
            case kAmpSelect:
          
              if( m_ampSelectBox->GetSelection( p2 ) ){
                m_generator.enableAmp( p2 );
              }
              else{ 
                m_generator.disableAmp( p2 );
              }
              break;

            case kSumSelect:
              
              if( m_sumSelectBox->GetSelection( p2 ) ){
                m_generator.enableSum( p2 );
              }
              else{ 
                m_generator.disableSum( p2 );
              }
              break;
              
            default:
              
              cout << "Unkown listbox!" << endl;
              break;
          } // end of listboxes
          break;
          
        case kCM_BUTTON:
          
          switch( p1 ){
              
            case kPlot:
              
              m_factory.drawPlot();
              break;
              
            case kSelectAllAmp:
              
              for( int i = 0; i < m_generator.uniqueAmplitudes().size(); ++i ){
                
                m_ampSelectBox->Select( i );
                m_generator.enableAmp( i );
              }
              break;
              
            case kClearAmp:
              
              for( int i = 0; i < m_generator.uniqueAmplitudes().size(); ++i ){
                
                m_ampSelectBox->Select( i, false );
                m_generator.disableAmp( i );
              }
              break;

            case kSelectAllSum:
              
              for( int i = 0; i < m_generator.uniqueSums().size(); ++i ){
                
                m_sumSelectBox->Select( i );
                m_generator.enableSum( i );
              }
              break;
              
            case kClearSum:
              
              for( int i = 0; i < m_generator.uniqueSums().size(); ++i ){
                
                m_sumSelectBox->Select( i, false );
                m_generator.disableSum( i );
              }
              break;
              
            case kclearCanv:
              m_factory.clearCanvas();
              break;
              
            case kExit:

              gROOT->ProcessLine( ".q" );
              break;
          }
          break;
          
        case kCM_COMBOBOX:
          
          switch( p1 ){
              
            case kChoosePlot:
              m_factory.setPlot( p2 );
              break;
              
            case kChoosePad:
              m_factory.setPad( p2 );
              break;
              
            case kChooseCanv:
              m_factory.chooseCanvLayout( p2 );
              break;
              
          }
          
          break;
      }
  }
	
  return true;
}

