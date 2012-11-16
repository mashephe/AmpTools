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

