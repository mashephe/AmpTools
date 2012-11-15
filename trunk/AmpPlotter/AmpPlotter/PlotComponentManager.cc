
#include <iostream>
#include <sstream>
#include <algorithm>

#include "AmpPlotter/PlotComponentManager.h"
#include "AmpPlotter/DataComponent.h"
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

