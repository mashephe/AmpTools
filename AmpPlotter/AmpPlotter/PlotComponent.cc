
#include "AmpPlotter/PlotComponent.h"
#include "AmpPlotter/ComponentGroup.h"

PlotComponent::~PlotComponent(){
	
	for( list<ComponentGroup*>::iterator group = m_groups.begin();
		 group != m_groups.end(); ++group ){
	
		// we want to avoid groups having stale pointers
		(**group).deleteComponent( this );
	}
}

bool
PlotComponent::isAvailable( void ) const {
	
	return m_dataOK;
}

bool
PlotComponent::isEnabled( void ) const {
	
	// for a component to be enabled, _all_ of the groups it is a member
	// of must also be enabled.
	
	bool enabled = true;
	
	for( list<ComponentGroup*>::const_iterator group = m_groups.begin();
		 group != m_groups.end(); ++group ){
		
		if( ! (**group).isEnabled() ){
			
			enabled = false;
		}
	}
	
	return enabled;
}

void
PlotComponent::addGroup( ComponentGroup* group ) {

	m_groups.push_back( group );
}

void
PlotComponent::deleteGroup( ComponentGroup* group ) {

	for( list<ComponentGroup*>::iterator groupItr = m_groups.begin();
		 groupItr != m_groups.end(); ++groupItr ){
	
		if( *groupItr == group ){
			
			m_groups.erase( groupItr );
		}
	}
}
