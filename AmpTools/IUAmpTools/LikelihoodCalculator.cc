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
#include <string>

#include "IUAmpTools/LikelihoodCalculator.h"
#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/NormIntInterface.h"
#include "IUAmpTools/ParameterManager.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MinuitParameterManager.h"
#include "MinuitInterface/MinuitParameter.h"

#include "IUAmpTools/report.h"
static const char* kModule = "LikelihoodCalculator";

#ifdef SCOREP
#include <scorep/SCOREP_User.h>
#endif

LikelihoodCalculator::LikelihoodCalculator( const IntensityManager& intenManager,
                                            const NormIntInterface& normInt,
                                            DataReader* dataReaderSignal,
                                            DataReader* dataReaderBkgnd,
                                            const ParameterManager& parManager ) :
MIFunctionContribution( parManager.fitManager() ),
m_intenManager( intenManager ),
m_normInt( normInt ),
m_dataReaderSignal( dataReaderSignal ),
m_dataReaderBkgnd( dataReaderBkgnd ),
m_firstDataCalc( true ),
m_firstNormIntCalc( true ),
m_sumBkgWeights( 0 ),
m_numBkgEvents( 0 ),
m_sumDataWeights( 0 ),
m_numDataEvents( 0 )
{
  
  m_hasBackground = ( dataReaderBkgnd != NULL );
  
// avoid caching data here in the constructor so that MPI-based classes that
// inherit from this class can have better control of what triggers the
// caching of data
 
  m_prodFactorArray = new double[2*intenManager.getTermNames().size()];
  m_normIntArray = normInt.normIntMatrix();
  m_ampIntArray = normInt.ampIntMatrix();
}

LikelihoodCalculator::~LikelihoodCalculator(){
  
  delete[] m_prodFactorArray;
}

double
LikelihoodCalculator::operator()(){
  
  // split the compuation of likelihood into data and normalization
  // integral sums -- this provides parallel implementations with the
  // methods they need to compute these terms individually
  
  return -2 * ( dataTerm() - normIntTerm() );
}

double
LikelihoodCalculator::numSignalEvents(){

  if( m_firstDataCalc ) dataTerm();
  return numDataEvents() - sumBkgWeights();
}


double
LikelihoodCalculator::normIntTerm(){
#ifdef SCOREP
SCOREP_USER_REGION_DEFINE( normIntTerm )                                                                                    
SCOREP_USER_REGION_BEGIN( normIntTerm, "normIntTerm", SCOREP_USER_REGION_TYPE_COMMON )
#endif
  
  // check to be sure we can actually perform a computation of the
  // normalization integrals in case we have floating parameters
  
  if( m_intenManager.hasTermWithFreeParam() && !m_normInt.hasAccessToMC() ){
    
    report( ERROR, kModule ) << "IntensityManager has terms with floating parameters\n"
         << "\tbut NormIntInterface has not been provided with MC." << endl;
    
    assert( false );
  }
  
  if( ( m_firstNormIntCalc && m_normInt.hasAccessToMC() ) ||
      ( m_intenManager.hasTermWithFreeParam() && !m_firstNormIntCalc ) ){

    m_normInt.forceCacheUpdate( true );
  }
  
  int n = m_intenManager.getTermNames().size();
  m_intenManager.prodFactorArray( m_prodFactorArray );
  
  double normTerm = 0;
  
  switch( m_intenManager.type() ){
      
    case IntensityManager::kAmplitude:
      
      for( int a = 0; a < n; ++a ){
        for( int b = 0; b <= a; ++b ){
          
          double thisTerm = 0;
          
          // we only need to compute the real part since the
          // imaginary part will sum to zero in the end
          //  want:  Re( V_a conj( V_b ) NI_a,b )
          
          double reVa = m_prodFactorArray[2*a];
          double imVa = m_prodFactorArray[2*a+1];
          double reVb = m_prodFactorArray[2*b];
          double imVb = m_prodFactorArray[2*b+1];
          double reNI = m_normIntArray[2*a*n+2*b];
          double imNI = m_normIntArray[2*a*n+2*b+1];
          
          thisTerm = ( reVa*reVb + imVa*imVb ) * reNI;
          thisTerm -= ( imVa*reVb - reVa*imVb ) * imNI;
          
          if( a != b ) thisTerm *= 2;
          
          normTerm += thisTerm;
        }
      }
      break;
      
    case IntensityManager::kMoment:
      
      break;
      
    default:
      
      report( ERROR, kModule ) << "Unknown IntensityManager type" << endl;
      assert( false );
      break;
  }
  
  m_firstNormIntCalc = false;

#ifdef SCOREP
  SCOREP_USER_REGION_END( normIntTerm )
#endif
  
  if( m_hasBackground ){
    
    // in this case let's match the number of observed events to sum of the
    // predicted signal and the provided background

    // this is the number of predicted signal and background events
    double nPred = normTerm + m_sumBkgWeights;

    normTerm = ( m_numDataEvents - m_sumBkgWeights ) * log( normTerm );
    normTerm += nPred - ( m_numDataEvents * log( nPred ) );
    
    return normTerm;
  }
  else{
    
    // this is the standard extended maximum likelihood technique that uses
    // a Poisson distribution of the predicted signal and the number of
    // observed events to normalize the fit components
    
    return normTerm;
  }

}

double
LikelihoodCalculator::dataTerm( bool suppressError ){
#ifdef SCOREP
SCOREP_USER_REGION_DEFINE( dataTerm )                                                                                    
SCOREP_USER_REGION_BEGIN( dataTerm, "dataTerm", SCOREP_USER_REGION_TYPE_COMMON )
#endif

  if( m_firstDataCalc ) {
    
    // first calculation -- need to load the data
    
    report( DEBUG, kModule ) << "Allocating Data and Amplitude Array in LikelihoodCalculator for "
      << m_intenManager.reactionName() << "..." << endl;
    
    m_ampVecsSignal.loadData( m_dataReaderSignal );
    m_ampVecsSignal.allocateTerms( m_intenManager, true );

    m_numDataEvents = m_ampVecsSignal.m_iNTrueEvents;
    m_sumDataWeights = m_ampVecsSignal.m_dSumWeights;
 
    if( m_ampVecsSignal.m_hasNonUnityWeights && m_hasBackground ){
    
      report( WARNING, kModule ) << "\n"
      << "****************************************************************\n"
      << "* WARNING: Events in the signal data sample have weights       *\n"
      << "*   that differ from 1 and a background data sample has been   *\n"
      << "*   provided.  This will produce incorrect results.  The       *\n"
      << "*   recommended solution is to set all signal weights to unity *\n"
      << "*   and put all background events in the background file with  *\n"
      << "*   weights such that the weighted sum mimics the background   *\n"
      << "*   contribution to the likelihood.                            *\n"
      << "****************************************************************\n" << endl;
    }
    
    if( m_hasBackground ){
          
      m_ampVecsBkgnd.loadData( m_dataReaderBkgnd );
      m_ampVecsBkgnd.allocateTerms( m_intenManager, true );

      if( m_ampVecsBkgnd.m_hasMixedSignWeights ){
        report( NOTICE, kModule ) << "\n"
        << "***************************************************************\n"
        << "* NOTICE:  Weights with both positive and negative signs were *\n"
        << "* detected in the background file.  This may be desirable for *\n"
        << "* some applications.  Older versions of AmpTools (v0.10.x and *\n"
        << "* prior) will not properly handle this case and will also not *\n"
        << "* print this notification to the screen.                      *\n"
        << "***************************************************************\n" << endl;
      }
      
      m_sumBkgWeights = m_ampVecsBkgnd.m_dSumWeights;
      
      // the extra boolean allows MPI jobs to suppress this check which
      // may fail on one of the follower nodes if the background sample
      // is sparse
      if( m_sumBkgWeights < 0 && !suppressError ){
        report( ERROR, kModule ) << "\n"
        << "****************************************************************\n"
        << "* ERROR: The sum of all background weights is negative.  This  *\n"
        << "*   implies a negative background in the signal region, which  *\n"
        << "*   unphysical.  The weighted sum of the background events     *\n"
        << "*   should represent the background contribution to the signal *\n"
        << "*   region.                                                    *\n"
        << "****************************************************************\n" << endl;
        assert( false );
      }
      
      m_numBkgEvents = m_ampVecsBkgnd.m_iNTrueEvents;
    }
    
    report( DEBUG, kModule ) << "\tDone." << endl;
  }
  
  double sumLnI = m_intenManager.calcSumLogIntensity( m_ampVecsSignal );
  
  // if there is a background file, we try to correct the sumLnI for background
  // by subtracting the contribution to the likelihood as derived from the
  // weighted background sample (the sign change here from v0.10.x and prior is
  // consistent with no longer forcing weights in the background sample to be negative)
  if( m_hasBackground ){

    sumLnI -= m_intenManager.calcSumLogIntensity( m_ampVecsBkgnd );
  }
  
  m_firstDataCalc = false;

#ifdef SCOREP
SCOREP_USER_REGION_END( dataTerm )
#endif
  
  return sumLnI;
}
