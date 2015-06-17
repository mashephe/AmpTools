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

#if !defined(KINEMATICS)
#define KINEMATICS

#include <vector>
#include <cassert>

#include "TLorentzVector.h"

using namespace std;

/**
 * This object contains a simple description of the particle kinematics
 * in a single event. It contains a list of four vectors and a weight
 * (optional) that can be used for a variety of different purposes.  The
 * object also has an event identifier that is automatically assigned and
 * is unique.
 *
 * Note that this class is where one of the few external dependencies in
 * AmpTools resides:  it utilizes the TLorentzVector class from CLHEP,
 * which is a convenient mechanism for manipulating four-vectors.
 *
 * \ingroup IUAmpTools
 */

class Kinematics
{
  
public:
  
  /**
   * The default constructor.  
   * 
   * This creates a Kinematics object with a unique identifier 
   * and a weight of 1.
   */
  Kinematics():
  m_eventID( Kinematics::m_globalEventID++ ), 
  m_weight( 1.0 ) {}
  
  /**
   * The constructor.
   *
   * This constructor generates a kinematics object from a vector
   * of TLorentzVector, which is a data type provided by the CLHEP
   * package.  A weight can be provided as an optional argument.
   *
   * \param[in] particleList a vector of TLorentzVector containing
   *   the four-momenta of all particles
   * \param[in] weight (optional) a weight to apply to this event
   */
  Kinematics( const vector< TLorentzVector >& particleList,
             float weight = 1.0 ) :
  m_eventID(  Kinematics::m_globalEventID++ ),
  m_particleList( particleList ),
  m_weight( weight ) { assert( particleList.size() <= kMaxParticles ); }
  
  /**
   * The destructor.
   */
  virtual ~Kinematics() {}
  
  /**
   * In many applications, fixed-size arrays are used to store kinematic
   * data.  Provided a single place to define the maximum number of four-
   * vectors that makeup an event.
   */
  
  enum { kMaxParticles = 6 };
  
  /**
   * Set the event ID.
   *
   * Event ID's are assigned at construction time and guaranteed to be unique.
   * However, the user can set the event ID to any integer using this 
   * method.  The AmpTools framework does not utilize the event ID -- it 
   * exists for historical reasons.
   *
   * \param[in] eventID the desired event ID
   *
   * \see eventID
   */
  void setEventID( int eventID ) { m_eventID = eventID; }

  /**
   * Set the list of particles.
   *
   * This method will replace the current list of four vectors with the
   * four vectors being passed in an argument to this function call.
   *
   * \param[in] particleList the vector of TLorentzVector that describes
   *   the kinematics for the event
   *
   * \see particle
   * \see particleList
   */
  void setParticleList( const vector< TLorentzVector >& particleList );
  
  /**
   * Set the weight of an event.
   *
   * This sets the weight of the event to be something other than 1 
   * (the default).  These weights are utilized in any computation of
   * the intensity or normalization integrals.
   *
   * \param[in] weight the desired weight
   *
   * \see weight
   */
  void setWeight( float weight ) { m_weight = weight; }
    
  /**
   * Get the list of four-vectors for an event.
   *
   * This method returns a vector of TLorentzVector corresponding
   * to the list of four-vectors for an event.
   *
   * \see setParticleList
   * \see particle
   */
  const vector<TLorentzVector>& particleList() const;

  /**
   * Get a particular four-vector.
   *
   * This method returns a TLorentzVector corresponding to the
   * four-momentum of the particle with the designated index.
   *
   * \param[in] index index of the particle
   *
   * \see particleList
   * \see setParticleList
   */
  const TLorentzVector& particle( unsigned int index ) const;

  /**
   * Get the event ID.
   *
   * This returns the identifier of the given event.
   *
   * \see setEventID
   */
  int eventID() const { return m_eventID; }
  
  /**
   * Get the weight.
   *
   * This returns the weight of a given event.
   *
   * \see setWeight
   */
  float weight() const { return m_weight; }
  
private:
  
  int m_eventID;
  vector<TLorentzVector> m_particleList;
  float m_weight;
  
  static int m_globalEventID;
};

#endif
