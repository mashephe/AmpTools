#if !defined( AMPLITUDE )
#define AMPLITUDE

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

#include <map>
#include <complex>
#include <vector>
#include <string>
#include <cassert>

#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/Term.h"
#include "GPUManager/GPUCustomTypes.h"

#ifdef GPU_ACCELERATION 
#include "cuda_runtime.h"
#include "GPUManager/CUDA-Complex.cuh"
class GPUManager;
#endif //GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;

/**
 * This class represents a user defined amplitude.  In its most abstract
 * sense it is a mechanism to turn a set of four vectors describing an
 * event into a complex number that represents the amplitude for 
 * some physics process.  These amplitudes can be added to the AmplitudeManager
 * as factors that build-up a decay amplitude for a particular final state.
 * 
 * It is recommended that the user factorize the physics amplitudes as much as
 * possible.  For example, a user may have an Amplitude class that defines
 * how to compute a Breit-Wigner and another Amplitude class that describes
 * some two-body angular distribution.  The complete decay amplitude for a
 * two-body resonance is given by a product of these two.  The factorization
 * ultimately improves optimization potential since only factors with
 * free parameters need to be recomputed at each fit iteration.
 *
 * In addition users should typically not permute identical particles in 
 * the construction of the amplitude.  Doing so "breaks" the factorization 
 * feature built into the AmplitudeManager since ( a1 + a2 ) * ( b1 + b2 ) is 
 * not equal to ( a1 b1 + a2 b2 ) for permutations 1 and 2 and amplitudes a and b.
 * Instead the AmplitudeManager provides a mechanism to permute idnetical
 * particles or charge conjugates for a particular amplitude.
 *
 * Finally the Amplitude class is intended to define methods only for
 * calculating amplitudes and should not be used as a data storage class.
 * The same instance of the Amplitude class will be used for calculations
 * over different sets of data.
 *
 * \ingroup IUAmpTools
 */

class Amplitude : public Term
{
  
public:
  
  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  Amplitude() : Term() { }
  
  /**
   * This constructor takes a list of arguments to inititialize an
   * amplitude and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  Amplitude( const vector< string >& args ) : Term( args ) { }
  
  /**
   * This is the destructor.
   */
  virtual ~Amplitude(){}
  
  /**
   * Must be overridden by the user to provide the name of the amplitude.
   * This is necessary to connect ConfigurationInfo to specific class
   * instances when the AmplitudeManager is being setup.
   */
  virtual string name() const = 0;
  

  /**
   * This must be overriden by the user and indicates how to convert a list
   * of strings (arguments) into a pointer to a new instance of the
   * users defined amplitude.
   *
   * The user can avoid writing this method by inheriting from the
   * UserAmplitude class (which derives from this class).
   *
   * \param[in] args a list of string arguments that may, for example, be
   * specified in a configuration file
   *
   *  \see UserAmplitude
   */
  virtual Amplitude* newAmplitude( const vector< string >& args ) const = 0;
  
  /**
   * A function that the user must write that indicates how an amplitude
   * can duplicate itself.  It returns a pointer to a new instance of the
   * Amplitude that behaves in exactly the same way as the original.
   *
   * The user can avoid writing this method by inheriting from the
   * UserAmplitude class (which derives from this class).
   *
   *  \see UserAmplitude
   */
  virtual Amplitude* clone() const = 0;
  
  /**
   * This loops over all events and calculates the amplitude for each event.
   * It does so by calling the user-defined calcAmplitude routine for each
   * permutation of particles.  The function is virtual since, in principle,
   * the user may choose to override it and perform the loop and computation
   * directly and more efficiently than using a separate function call
   * for each event.
   *
   * \param[in] pdData a pointer to the array of data.  This is a long list
   * of GDoubles E, px, py, pz repeated sequentially for each particle in 
   * the event and then the series repeated for each event in the data set.
   *
   * \param[out] pdAmps a pointer to a block of memory where the calculated
   * amplitudes should go.  The values are stored as real and imaginary parts
   * repeated for each permutation and then for each event.  
   * Total length of memory needed is 2 * sizeof( GDDouble ) * N_evts * N_perms
   *
   * \param[in] iNEvents the number of events in the data set
   *
   * \param[in] pvPermutations a pointer to the vector that contains the various
   * permutations.  The size of this vector is the number of permutations and
   * the size of one element of this vector is the number of particles.
   *
   * \param[in] pdUserVars is a pointer to the block of user data that is
   *  optionally filled by the user with precalculated variables in advance
   *  of the computation
   *
   * \see calcAmplitudeAll
   * \see AmplitudeManager::addAmpPermutation
   */
  virtual void calcAmplitudeAll( GDouble* pdData, GDouble* pdAmps, int iNEvents,
                                const vector< vector< int > >* pvPermutations,
                                GDouble* pdUserVars = 0 ) const;
  
  
  /**
   * This is the user-defined function that computes a single complex amplitude
   * for a set of four-vectors that describe the event kinematics.  As discussed
   * above this function should be factorized as much as possible and not
   * include permutations of particles.  The user must override this function
   * in his or her amplitude class.
   *
   * \param[in] pKin a pointer to a single event.  pKin[0][0-3] define E, px,
   * py, pz for the first particle, pKin[1][0-3] for the second, and so on
   *
   */
  virtual complex< GDouble > calcAmplitude( GDouble** pKin ) const ;
  
  /**
   * This is the user-defined function that computes a single complex amplitude
   * for a set of four-vectos that describe the event kinematics.  
   * 
   * For the user to utilize user-defiend data in the amplitude calculation,
   * this function must be overridden by the derived class.  Either this function
   * or the function above must be defined for any Amplitude class.
   *
   * \param[in] pKin a pointer to a single event.  pKin[0][0-3] define E, px,
   * py, pz for the first particle, pKin[1][0-3] for the second, and so on
   *
   * \param[in] userVars is an optional pointer to the user data block associated
   * with this event and this permutation of particles.  It can be used to store
   * intermediate portions of the calculation in the case that calcAmplitude
   * must be called multiple times during the course of a fit.  The userVars
   * memory block is filled in calcUserVars.
   */
  virtual complex< GDouble > calcAmplitude( GDouble** pKin,
                                            GDouble* userVars ) const;
  
  /**
   * \overload
   *
   * This is a user-friendly interface to calcAmplitude used for diagnostics.
   * It takes a pointer to a Kinematics object, converts it into a GDouble** format,
   * then calls the standard calcAmplitude method.
   *
   * \param[in] pKin a pointer to a Kinematics object
   *
   * \param[in] userVars (optional) address of userVars memory block that will
   * get passed onto the standard calcAmplitude method
   *
   * \see calcAmplitude
   */
  
  complex< GDouble > calcAmplitude( const Kinematics* pKin,
                                    GDouble* userVars = 0 ) const;
  
  
  /**
   * \overload
   *
   * This is a user-friendly interface to calcAmplitude used for diagnostics.
   * It takes a pointer to a Kinematics object, rearranges it according to the
   * permutation specified, converts it into a GDouble** format,
   * then calls the standard calcAmplitude method.
   *
   * \param[in] pKin a pointer to a Kinematics object
   * \param[in] permutation a vector of permuted particle indices
   *
   * \param[in] userVars (optional) address of userVars memory block that will
   * get passed onto the standard calcAmplitude method
   *
   * \see calcAmplitude
   * \see AmplitudeManager::addAmpPermutation
   */
  
  complex< GDouble > calcAmplitude( const Kinematics* pKin,
                                    const vector < int >& permutation,
                                    GDouble* userVars = 0 ) const;
  
 
#ifdef GPU_ACCELERATION 
  
  /**
   * If GPU_ACCELERATION flag is set this is the member function that sets
   * the current permutation and then calls the user-defined routine
   * to launch the GPU kernel.  It is the GPU analog of calcAmplitudeAll.
   */
  virtual void calcAmplitudeGPU( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO,
                                const vector< int >& perm ) const;
  
  /**
   * The user override this route and use it pass any parameters to a global
   * C function that actually launches the GPU kernel
   */
  virtual void launchGPUKernel( dim3 dimGrid, dim3 dimBlock, GPU_AMP_PROTO ) const {
    
    cout << "\nNo GPU function for calculating " << name() << " is defined." << endl;
    assert( false );
  }
  
#endif //GPU_ACCELERATION
  
protected:
  
  /**
   * A helper function so that the user's calcAmplitude routine can access
   * which permutation has been passed to it.  In general the user does not
   * need to know this as the calcAllAmplitudes routine permutes the particles
   * before calling the user-defined calcAmplitude.  However in some case
   * the user code may want to include additional factors, e.g., isospin 
   * Clebsch-Gordan coefficients that can only be computed if the the
   * permutation of the particles is known.
   *
   * \see calcAllAmplitudes
   * \see calcAmplitude
   */
  inline const vector< int >& getCurrentPermutation() const { return m_currentPermutation; }
  
private:
  
  mutable vector< int > m_currentPermutation;
};


#endif
