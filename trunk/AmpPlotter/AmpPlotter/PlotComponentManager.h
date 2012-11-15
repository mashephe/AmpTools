#if (!defined PLOTCOMPONENTMANAGER)
#define PLOTCOMPONENTMANAGER

#include <vector>
#include <string>
#include <map>

#include "AmpPlotter/PlotComponent.h"
#include "AmpPlotter/ComponentGroup.h"
#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class PlotComponentManager
{
	
public:
  
	PlotComponentManager( const vector< string >& reactionVec,
                        PlotGenerator& pltGen,
                        bool boringPlots = false );
	
	void toggleData( int dataIndex );
  
	const ComponentGroup* dataGroup() const;
	const ComponentGroup* genMCGroup() const;
	const ComponentGroup* accMCGroup() const;
	
	const vector<ComponentGroup*>& reactionGroups() const;
  
	void enableData()  { m_dataGroup->enable();  }
	void disableData() { m_dataGroup->disable(); }
  
  void enableGenMC()  { m_genMCGroup->enable();  }
	void disableGenMC() { m_genMCGroup->disable(); }
  
  void enableAccMC()  { m_accMCGroup->enable();  }
	void disableAccMC() { m_accMCGroup->disable(); }
  
  void enableReaction( int reactionIndex );
	void disableReaction( int reactionIndex );
  
private:
  
  PlotGenerator& m_generator;
  
	// this groups all components by reaction
	vector<ComponentGroup*> m_reactionGroup;
  
	// these groups allow for turning off and on the data, 
	// generated, and accepted MC
  ComponentGroup* m_dataGroup;
  ComponentGroup* m_genMCGroup;
  ComponentGroup* m_accMCGroup;
  
	bool m_boringPlots;
  
  map< string, int > m_reactionIndex;
  vector< string > m_reactionVec;
};

#endif
