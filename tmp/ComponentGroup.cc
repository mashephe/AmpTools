
#include "AmpPlotter/ComponentGroup.h"

ComponentGroup::~ComponentGroup(){

	for( std::list<PlotComponent*>::iterator comp = m_componentList.begin();
		 comp != m_componentList.end(); ++comp ){
	
		// we don't want components with stale pointers
		(**comp).deleteGroup( this );
	}
}

void
ComponentGroup::addComponent( PlotComponent* component ){
	
	m_componentList.push_back( component );
	component->addGroup( this );
}

void
ComponentGroup::deleteComponent( PlotComponent* component ){
	
	for( std::list<PlotComponent*>::iterator comp = m_componentList.begin();
		 comp != m_componentList.end(); ++comp ){
	
		if( *comp == component ){
			
			m_componentList.erase( comp );
		}
	}
	component->deleteGroup( this );
}

void ComponentGroup::enable( void ) { m_isEnabled = true; }
void ComponentGroup::disable( void ) { m_isEnabled = false; }

bool ComponentGroup::isEnabled( void ) const { return m_isEnabled; }
std::string ComponentGroup::title( void ) const { return m_title; }

