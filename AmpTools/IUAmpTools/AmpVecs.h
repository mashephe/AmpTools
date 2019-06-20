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

#include <map>
#include <string>

#ifdef GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif

using namespace std;

class DataReader;
class IntensityManager;
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
  unsigned long long m_iNEvents;
  
  /**
   * An integer that stores the true number of events.  For GPU calculations
   * it is necessary to pad iNEvents up to the next power of 2.  This integer
   * stores the actual number of unique events.
   */
  unsigned long long m_iNTrueEvents;
  
  /**
   * A double that stores the absolute value of the sum of the weights.  (For
   * cases where all weights are unity, this is simply the number of true 
   * events.)
   */
  
  double m_dAbsSumWeights;
  
  /**
   * An integer that stores the number of particles in the final state.
   */
  unsigned int m_iNParticles;
  
  /**
   * An integer that stores the number of amplitudes for a particular
   * configuration of the AmplitudeManager.
   */
  unsigned int m_iNTerms;
  
  /**
   * An integer that is number of doubles required to store all factors
   * and permutations for any term for an event.
   */
  unsigned int m_maxFactPerEvent;
  
  /**
   * An integer that is the number of doubles required to store all
   * (optionally) user-calculated data per event.
   */
  unsigned int m_userVarsPerEvent;

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
   * An array of length iNEvents * m_userVarsPerEvent that is used to 
   * store user-calculated and cached data to expedite subsequent
   * amplitude calculations.
   */
  GDouble* m_pdUserVars;
  
  /**
   * An array of length 2 * iNTerms * iNTerms that holds the sums of
   * the prodcuts of A_i A_j* for all data events and all terms i,j
   * that is useful for calculating normalization integrals.
   */
  
  GDouble* m_pdIntegralMatrix;
  
  /**
   * An array of length iNEvents that stores the intensity for each event.
   */
  GDouble* m_pdIntensity;
  
  /**
   * A boolean that tracks if m_pdAmps and m_pdAmpFactors are filled
   * for the current data set in m_pdData.  This variable will be
   * managed by classes outside of AmpVecs.  It can only be invalidated
   * (set to false) by the AmpVecs class itself.
   */
  bool m_termsValid;

  /**
   * A boolean that tracks if integrals are filled for the current data
   * set in m_pdData and terms in m_pdAmps.  This variable will be
   * managed by classes outside of AmpVecs.  It can only be invalidated
   * (set to false) by the AmpVecs class itself.
   */
  bool m_integralValid;
  
  /**
   * This is a map from amplitude identifer to the location in memory
   * where user data for that amplitude exists.  It is
   * utilized by the AmplitudeManager, but the values are tied
   * to the data set so it resides in the AmpVecs struct.
   */
  map< string, unsigned long long > m_userVarsOffset;

  
#ifdef GPU_ACCELERATION
  /**
   * The GPU Manager for this data set.
   */
  GPUManager m_gpuMan;
  
#endif
  
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
   * This routine allocates space to store the calculated terms and
   * (optionally) the calculated intensities.  It should be passed a reference
   * to the IntensityManager AFTER the IntensityManager has been setup and
   * configured to compute the intensity.
   *
   * \param[in] intenMan a reference to the IntensityManager to use as a guide
   * for allocating the arrays for terms, factors, and intensities
   *
   * \param[in] bAllocIntensity if set to true this will allocate space for
   * the intensity calculation
   */
  void allocateTerms( const IntensityManager& intenMan,
                      bool bAllocIntensity = false );
  
#ifdef GPU_ACCELERATION
  /**
   * This will allocate CPU memory for storage and copying of amplitudes
   * and amplitude factors from the GPU.  For most production operations
   * this wastes CPU memory and is unnecessary since all amplitudes
   * are maintained on the GPU.  However for some debugging operations,
   * this functionality is useful.
   *
   * \param[in] intenMan a reference to the IntensityManager to use as a guide
   * for allocating the arrays for terms, factors, and intensities
   */
  void allocateCPUAmpStorage( const IntensityManager& intenMan );
#endif
  
  /**
   * This routine uses the pointer to the data reader that is provided to 
   * allocate and fill the array of data and weights.  The function will
   * first reset the source and then get events from the data reader until
   * a null pointer is provided (signaling the end of the source).
   *
   * \param[in] pDataReader a pointer to a user-defined data reader
   * \param[in] bForceNegativeWeight sets weight = -abs( weight )
   * 
   * \see DataReader::resetSource
   * \see DataReader::getEvent
   */
  void loadData( DataReader* pDataReader, bool bForceNegativeWeight = false );
  
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
   * \param[in] bForceNegativeWeight sets weight = -abs( weight )
   *
   * \see loadData
   */
  void loadEvent( const Kinematics* pKinematics, unsigned long long iEvent = 0,
                  unsigned long long iNTrueEvents = 1,
                  bool bForceNegativeWeight = false );
  
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
  
  /**
   * This clears only the four vectors from memory.  It can be used
   * in the case all amplitudes depend only on user data.
   */
  void clearFourVecs();
};

#endif 

