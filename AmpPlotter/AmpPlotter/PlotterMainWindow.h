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

#if (!defined PLOTTERMAINWINDOW)
#define PLOTTERMAINWINDOW

#include <vector>

#include "TROOT.h"
#include "TGWindow.h"
#include "TGButton.h"
#include "TGFrame.h"
#include "TGComboBox.h"
#include "TGListBox.h"
#include "TGTab.h"

#include "AmpPlotter/PlotFactory.h"
#include "AmpPlotter/PlotComponentManager.h"
#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class PlotterMainWindow : public TGMainFrame 
{
		
 public:
	
  enum { kWidth  = 690 };
  enum { kHeight = 400 };
	
  enum { kIndexMask     = 0x000000FF };
  enum { kButtonMask    = 0xFFFFFF00 };

  enum { kReaction      = 0x000100 };
  enum { kType          = 0x000200 };
	
  enum { kExit, kPlot, kSelectAllAmp, kClearAmp, kSelectAllSum, kClearSum,
         kChoosePlot, kChoosePad , kclearCanv, kChooseCanv};
  enum { kData, kBkgnd, kAccMC, kGenMC };
  enum { kSumSelect, kAmpSelect };  
  
  
  PlotterMainWindow( const TGWindow*, PlotFactory& );
	
  ~PlotterMainWindow(){}
  
  bool ProcessMessage( long mes, long p1, long p2 );


 private:
		
  PlotterMainWindow( const PlotterMainWindow& );
  const PlotterMainWindow operator=( const PlotterMainWindow& );
		
  // fill these with empty functions -- there is no need to save
  // the state of the plotter as a stream
  void SavePrimitive(std::basic_ostream<char,
                     std::char_traits<char> >&, char const*){}
  void SavePrimitiveSubframes(std::basic_ostream<char,
                              std::char_traits<char> >&, char const*){}
  
  TGComboBox* m_padButton;
  TGComboBox* m_canvButton;
  TGTextButton* m_exitButton;
  TGTextButton* m_plotButton;
  TGTextButton* m_selectAllAmpButton;
  TGTextButton* m_clearAmpButton;
  TGTextButton* m_selectAllSumButton;
  TGTextButton* m_clearSumButton;
  TGTextButton* m_clearCanvButton;
	
  vector<TGCheckButton*> m_reactionButtons;
  vector<TGCheckButton*> m_typeButtons;
	
  TGVerticalFrame*   m_mainFrame;
  TGHorizontalFrame* m_menuFrame;
  TGHorizontalFrame* m_restFrame;

  TGVerticalFrame*   m_rightFrame;
  TGVerticalFrame*   m_leftFrame;
	
  TGTab* m_ampSumTab;
  
  TGGroupFrame* m_canvFrame;
  TGGroupFrame* m_reactionFrame;
  TGGroupFrame* m_typeFrame;
  TGCompositeFrame* m_ampFrame;
  TGCompositeFrame* m_sumFrame;
	
  TGHorizontalFrame* m_plotSelectFrame;
  TGComboBox*        m_plotSelectBox;
  TGListBox*         m_ampSelectBox;
  TGListBox*         m_sumSelectBox;
	
  TGHorizontalFrame* m_buttonFrame;
  TGHorizontalFrame* m_selClrAmpFrame;
  TGHorizontalFrame* m_selClrSumFrame;
	
  PlotFactory& m_factory;
	
  PlotComponentManager& m_componentManager;
  PlotGenerator& m_generator;
};

#endif
