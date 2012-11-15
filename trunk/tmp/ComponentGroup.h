#if !(defined COMPONENTGROUP)
#define COMPONENTGROUP

#include <list>

#include "AmpPlotter/PlotComponent.h"

using namespace std;

class ComponentGroup
{
	
public:
	
	ComponentGroup( string title ) : m_title( title ), m_isEnabled( false ) {}
	virtual ~ComponentGroup();
	
	void addComponent( PlotComponent* component );
	void deleteComponent( PlotComponent* component );
	
	void enable();
	void disable();
	
	bool isEnabled() const;
	string title() const;
	
	const list<PlotComponent*>& componentList() const { 
		
		return m_componentList;
	}
	
private:
		
	ComponentGroup( const ComponentGroup& );                   
	const ComponentGroup& operator=( const ComponentGroup& ); 
	
	string m_title;
	list<PlotComponent*> m_componentList;
	bool m_isEnabled;
};

#endif
