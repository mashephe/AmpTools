#ifndef AMPVECS
#define AMPVECS

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

#include "GPUManager/GPUCustomTypes.h"

class DataReader;
class AmplitudeManager;
class Kinematics;

/**
 * This class is intended to be a helper structure that serves as a place
 * to pack data, amplitudes, and intensities for a particular data set.
 * Internally these data are stored as flat arrays to improve speed of the code.
 * The data-packing in this structure is not the most user friendly, but the
 * design is such that the user never really needs to interact with the contents
 * of this structure.
 *
 * \ingroup IUAmpTools
 */

//A helper struct for parameter passing
struct AmpVecs
{
  
  /**
   * An integer that stores the number of events.
   */
	int m_iNEvents;
  
  /**
   * An integer that stores the true number of events.  For GPU calculations
   * it is necessary to pad iNEvents up to the next power of 2.  This integer
   * stores the actual number of unique events.
   */
	int m_iNTrueEvents;
  
  /**
   * An integer that stores the number of particles in the final state.
   */
	int m_iNParticles;
  
  /**
   * An integer that stores the number of amplitudes for a particular
   * configuration of the AmplitudeManager.
   */
  int m_iNAmps;
  
  /**
   * An integer that stores the number of unique amplitudes and permutations
   * for a particular configuration of the AmplitudeManager.  This is equivalent
   * to the number calcAmplitude calls that must be made for each event in order
   * to compute the intensity.
   *
   * \see Amplitude::calcAmplitude
   */
  int m_iNAmpFactorsAndPerms;
	
  /**
   * An array of length 4 * iNEvents * iNParticles that stores the four-vectors
   * for the entire data set.
   */
	GDouble* m_pdData;

	/**
   * An array of length iNEvents that stores the event weights for each event.
   */
  GDouble* m_pdWeights;
	
  /**
   * An array of length 2 * iNAmps * iNEvents that stores the real and imaginary
   * parts of the complete decay amplitude (product of factors) for each event.
   */
	GDouble* m_pdAmps;

  /**
   * An array of length 2 * iNAmpFactorsAndPerms * iNEvents that stores the
   * real and imaginary parts for every factor of the decay amplitude and 
   * every permutation.
   */
	GDouble* m_pdAmpFactors;	
	
  /**
   * An array of length iNEvents that stores the intensity for each event.
   */
	GDouble* m_pdIntensity;
	
  /**
   * The constructor.  All array pointers are set to zero in the constructor.
   * No memory is allocated at construction time.
   */
	AmpVecs();
  
  /**
   * This destructor frees all of the allocated memory that holds all of the
   * data and amplitude calculations.
   *
   * \see deallocAmpVecs
   */
	~AmpVecs(){ deallocAmpVecs(); }

  /**
   * This routine allocates space to store the calculated amplitudes and
   * (optionally) the calculated intensities.  It should be passed a reference
   * to the AmplitudeManager after the AmplitudeManager has been setup and
   * configured to compute the intensity as it needs to know, for example,
   * the number of permutations for each amplitude.
   *
   * \param[in] ampMan a reference to the AmplitudeManager to use as a guide
   * for allocating the arrays for amplitudes, factors, and intensities
   *
   * \param[in] bAllocIntensity if set to true this will allocate space for
   * the intensity calculation
   */
  void allocateAmps( const AmplitudeManager& ampMan, bool bAllocIntensity = false );
 
  /**
   * This routine uses the pointer to the data reader that is provided to 
   * allocate and fill the array of data and weights.  The function will
   * first reset the source and then get events from the data reader until
   * a null pointer is provided (signaling the end of the source).
   *
   * \param[in] pDataReader a pointer to a user-defined data reader
   * 
   * \see DataReader::resetSource
   * \see DataReader::getEvent
   */
  void loadData( DataReader* pDataReader );

  /**
   * This routine fills the arrays of data and weights event by event
   * rather than all at once.  If the data arrays pointers are null, then
   * memory is allocated according to the iNTrueEvents argument that is
   * passed in -- this should be set to the number of events the user plans
   * to load.  The iEvent argument is the location in the data array where
   * the event should be loaded.
   *
   * \param[in] pKinematics a pointer to a Kinematics object
   * \param[in] iEvent an event counter (increment this for every event loaded)
   * \param[in] iNTrueEvents the total number of events to be loaded
   *
   * \see loadData
   */
  void loadEvent( const Kinematics* pKinematics, int iEvent = 0, int iNTrueEvents = 1 );

  /**
   * A helper routine to get an event i from the array of data and weights.
   * This routine allows the AmpVecs class to behave in the same way that
   * a DataReader would.
   *
   * \see DataReader::getEvent
   */
  Kinematics* getEvent( int i );

  /**
   * This deallocates all allocated memory, which effective erases the entire
   * contents (data and amplitudes) of AmpVecs
   */
  void deallocAmpVecs();
};

#endif 

