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

	virtual TH1F deliverPlot( unsigned int plotIndex, 
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
