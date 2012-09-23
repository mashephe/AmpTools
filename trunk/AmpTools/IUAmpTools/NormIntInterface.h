#if !(defined NORMINTINTERFACE)
#define NORMINTINTERFACE

//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
// 
// Copyright Trustees of Indiana University 2010, all rights reserved
// 
// This software written by Matthew Shepherd, Ryan Mitchell, and 
//                  Hrayr Matevosyan at Indiana University, Bloomington
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

#include <complex>
#include <string>
#include <map>
#include <vector>

#include "GPUManager/GPUCustomTypes.h"
#include "IUAmpTools/AmpVecs.h"

class Kinematics;
class AmplitudeManager;
class DataReader;

using namespace std;

class NormIntInterface
{
	
public:
	
  NormIntInterface();
	NormIntInterface( const string& normIntFile );
	NormIntInterface( DataReader* genMCData, DataReader* accMCData, 
                   const AmplitudeManager& ampManager );
	
  istream& loadNormIntCache( istream& in );
  
	int numGenEvents() const { return m_nGenEvents; }
	int numAccEvents() const { return m_nAccEvents; }
	
  // does this interface have the ability to recalculate integrals?
  bool hasAccessToMC() const;
  
	// this integral folds in detector acceptance
	virtual complex< double > normInt( string amp, string conjAmp, bool forceUseCache = false ) const;
	bool hasNormInt( string amp, string conjAmp ) const;

	// this purely the integral of the amplitude -- perfect acceptance
	virtual complex< double > ampInt( string amp, string conjAmp, bool forceUseCache = false ) const;
	bool hasAmpInt( string amp, string conjAmp ) const;

  // needs to be virtual so parallel implementations can properly
  // override this function
  virtual void forceCacheUpdate( bool normIntOnly = false ) const;

	void exportNormIntCache( const string& fileName, bool renormalize = false ) const;
  void exportNormIntCache( ostream& output, bool renormalize = false ) const;
  
	void setGenEvents( int events ) { m_nGenEvents = events; }
	void setAccEvents( int events ) { m_nAccEvents = events; }
	
 protected:

	// protected helper functions for parallel implementations

  map< string, map< string, complex< double > > > getAmpIntegrals() const;
  map< string, map< string, complex< double > > > getNormIntegrals() const;


  void setAmpIntegral( string ampName, string cnjName,
                       complex< double > val ) const;
  void setNormIntegral( string ampName, string cnjName,
                        complex< double > val ) const;

  const AmplitudeManager* ampManager() const { return m_pAmpManager; }
  
 private:
  
	const AmplitudeManager* m_pAmpManager;
  
  DataReader* m_accMCReader;
  DataReader* m_genMCReader;
  
	int m_nGenEvents;
  int m_nAccEvents;

  mutable bool m_emptyNormIntCache;
  mutable bool m_emptyAmpIntCache;
  
  // needed to cache accepted MC data for NI recalculation
  mutable AmpVecs m_mcVecs;

	// amp -> amp* -> value
	mutable map< string, map< string, complex< double > > > m_normIntCache;
	mutable map< string, map< string, complex< double > > > m_ampIntCache;
};

inline istream& operator>>( istream& input, NormIntInterface& normInt ){
  
  return normInt.loadNormIntCache( input );
}

#endif
