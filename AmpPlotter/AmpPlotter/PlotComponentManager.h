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
  const ComponentGroup* bkgndGroup() const;
	const ComponentGroup* genMCGroup() const;
	const ComponentGroup* accMCGroup() const;
	
	const vector<ComponentGroup*>& reactionGroups() const;
  
	void enableData()  { m_dataGroup->enable();  }
	void disableData() { m_dataGroup->disable(); }

  void enableBkgnd()  { m_bkgndGroup->enable();  }
  void disableBkgnd() { m_bkgndGroup->disable(); }
  
  void enableGenMC()  { m_genMCGroup->enable();  }
	void disableGenMC() { m_genMCGroup->disable(); }
  
  void enableAccMC()  { m_accMCGroup->enable();  }
	void disableAccMC() { m_accMCGroup->disable(); }
  
  void enableReaction( int reactionIndex );
	void disableReaction( int reactionIndex );
  
  void weightMC( bool value ) { m_generator.setWeightMCByIntensity( value ); }
  
private:
  
  PlotGenerator& m_generator;
  
	// this groups all components by reaction
	vector<ComponentGroup*> m_reactionGroup;
  
	// these groups allow for turning off and on the data, 
	// background, generated, and accepted MC
  ComponentGroup* m_dataGroup;
  ComponentGroup* m_bkgndGroup;
  ComponentGroup* m_genMCGroup;
  ComponentGroup* m_accMCGroup;
  
	bool m_boringPlots;
  
  map< string, int > m_reactionIndex;
  vector< string > m_reactionVec;
};

#endif
