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

#include <sys/time.h>
#include <cassert>
#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MinuitParameterManager.h"
#include "MinuitInterface/MinuitParameter.h"

LikelihoodCalculator::LikelihoodCalculator( const AmplitudeManager& ampManager,
                                            const NormIntInterface& normInt,
                                            DataReader& dataReader,
                                            const ParameterManager& parManager ) :
MIFunctionContribution( parManager.fitManager() ),
m_ampManager( ampManager ),
m_normInt( normInt ),
m_dataReader( dataReader ),
m_firstPass( true )
{
  
// avoid caching data here in the constructor so that MPI-based classes that
// inherit from this class can have better control of what triggers the
// caching of data
  
}

double
LikelihoodCalculator::operator()(){
  
  // split the compuation of likelihood into data and normalization
  // integral sums -- this provides parallel implementations with the
  // methods they need to compute these terms individually
  
  return ( -2 * ( dataTerm() - normIntTerm() ) );
}


double
LikelihoodCalculator::normIntTerm(){
  
  // check to be sure we can actually perform a computation of the
  // normalization integrals in case we have floating parameters
  
  if( m_ampManager.hasAmpWithFreeParam() && !m_normInt.hasAccessToMC() ){
    
    cout << "ERROR: AmplitudeManager has amplitudes with floating parameters\n"
    << "       but NormIntInterface has not been provided with MC." << endl;
    
    assert( false );
  }
  
  // the normalization integral table won't change during the course
  // of this loop, so after we have retrieved one integral (which will
  // force recalculation of the table) it is safe to use the cached table
  
  bool useCachedIntegrals = false;
  
  complex< double > normTerm( 0, 0 );
  vector< string > ampNames = m_ampManager.getAmpNames();
  for( vector< string >::iterator amp = ampNames.begin();
      amp != ampNames.end(); ++amp ){
    
    for( vector< string >::iterator conjAmp = ampNames.begin();
        conjAmp != ampNames.end(); ++conjAmp ){
      
      complex< double > thisTerm( 0, 0 );

      thisTerm = ( m_ampManager.productionAmp( *amp ) * 
                   conj( m_ampManager.productionAmp( *conjAmp ) ) *
                  m_normInt.normInt( *amp, *conjAmp, useCachedIntegrals ) );
              
      if( m_ampManager.ampsAreRenormalized() ){
        
        thisTerm /=
          sqrt( real( m_normInt.ampInt( *amp, *amp, useCachedIntegrals ) ) *
                real( m_normInt.ampInt( *conjAmp, *conjAmp, useCachedIntegrals ) ) );
      }
      
      normTerm += thisTerm;
      
      useCachedIntegrals = true;
    }
  }
  
  // normTerm should be purely real by construction
  return real( normTerm );
}

double
LikelihoodCalculator::dataTerm(){
    
  if( m_firstPass ) { 
    
    // first calculation -- need to load the data
    
    cout << "Allocating Data and Amplitude Array in LikelihoodCalculator for " 
         << m_ampManager.reactionName() << "..." << endl;
    
    m_ampVecs.loadData( &m_dataReader );
    m_ampVecs.allocateAmps( m_ampManager, true );  
    cout << "\tDone." << endl;
  }
  
  double sumLnI = m_ampManager.calcSumLogIntensity( m_ampVecs, m_firstPass );
  
  m_firstPass = false;
  
  return sumLnI;
}
