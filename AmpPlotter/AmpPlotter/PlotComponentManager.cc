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
#include <sstream>
#include <algorithm>

#include "AmpPlotter/PlotComponentManager.h"
#include "AmpPlotter/DataComponent.h"
#include "AmpPlotter/BkgndComponent.h"
#include "AmpPlotter/AccMCComponent.h"
#include "AmpPlotter/GenMCComponent.h"

#include "TString.h"

using namespace std;

PlotComponentManager::PlotComponentManager( const vector< string >& reactionVec,
                                           PlotGenerator& pltGen,
                                           bool boringPlots ) :
m_generator( pltGen ),
m_boringPlots( boringPlots ),
m_reactionGroup( 0 ),
m_dataGroup( new ComponentGroup( "Data" ) ),
m_bkgndGroup( new ComponentGroup( "Background" ) ),
m_genMCGroup( new ComponentGroup( "Generated Monte Carlo" ) ),
m_accMCGroup( new ComponentGroup( "Accepted Monte Carlo" ) ),
m_reactionVec( reactionVec )
{
  for( int i = 0; i < reactionVec.size(); ++i ){
    
    m_reactionIndex[reactionVec[i]] = i;
    m_reactionGroup.push_back( new ComponentGroup( reactionVec[i] ) );
    
    DataComponent* data = new DataComponent( reactionVec[i], reactionVec[i], pltGen );
    data->setFillStyle( 0 );
    m_dataGroup->addComponent( data );
    m_reactionGroup[i]->addComponent( data );
    
    BkgndComponent* bkgnd = new BkgndComponent( reactionVec[i], reactionVec[i], pltGen );
    bkgnd->setFillStyle( 3004 );
    bkgnd->setFillColor( 47 );
    m_bkgndGroup->addComponent( bkgnd );
    m_reactionGroup[i]->addComponent( bkgnd );
    
    AccMCComponent* accMC = new AccMCComponent( reactionVec[i], reactionVec[i], pltGen );
    accMC->setFillStyle( 1001 );
    accMC->setFillColor(   32 );
    m_accMCGroup->addComponent( accMC );
    m_reactionGroup[i]->addComponent( accMC );
    
    GenMCComponent* genMC =
    new GenMCComponent( reactionVec[i], reactionVec[i], pltGen );
    genMC->setFillStyle( 1001 );
    genMC->setFillColor(   41 );
    m_genMCGroup->addComponent( genMC );
    m_reactionGroup[i]->addComponent( genMC );
  }
  
  
}

void
PlotComponentManager::enableReaction( int index ){
  
  ComponentGroup* group = m_reactionGroup[index];
  
  //	cout << "Enabling display of reaction:  " << group->title() << endl;
  group->enable();
  m_generator.enableReaction( m_reactionVec[index] );
}

void
PlotComponentManager::disableReaction( int index ){
  
  ComponentGroup* group = m_reactionGroup[index];
  
  //	cout << "Disabling display of reaction:  " << group->title() << endl;
  group->disable();
  m_generator.disableReaction( m_reactionVec[index] );
}

const ComponentGroup*
PlotComponentManager::dataGroup() const {
  
  return m_dataGroup;
}

const ComponentGroup*
PlotComponentManager::bkgndGroup() const {
  
  return m_bkgndGroup;
}

const ComponentGroup*
PlotComponentManager::genMCGroup() const {
  
  return m_genMCGroup;
}

const ComponentGroup*
PlotComponentManager::accMCGroup() const {
  
  return m_accMCGroup;
}

const std::vector<ComponentGroup*>&
PlotComponentManager::reactionGroups() const {
  
  return m_reactionGroup;
}

