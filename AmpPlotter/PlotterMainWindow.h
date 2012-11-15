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
  enum { kData, kAccMC, kGenMC };
  enum { kSumSelect, kAmpSelect };  
  
  
  PlotterMainWindow( const TGWindow*, PlotFactory& );
	
  bool ProcessMessage( long mes, long p1, long p2 );


 private:
		
  PlotterMainWindow( const PlotterMainWindow& );
  const PlotterMainWindow operator=( const PlotterMainWindow& );
		
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
