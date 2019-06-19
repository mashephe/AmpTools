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
#include <set>

#include "IUAmpTools/Kinematics.h"
#include "GPUManager/GPUCustomTypes.h"

#ifdef GPU_ACCELERATION 
#include "cuda_runtime.h"
#include "GPUManager/CUDA-Complex.cuh"
class GPUManager;
#endif //GPU_ACCELERATION

using std::complex;
using namespace std;

class Kinematics;
class AmpParameter;

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

class Amplitude
{
  
public:
  
  /**
   * The default constructor.  The user's derived class should contain a
   * default constructor that calls this constructor.
   */
  Amplitude( ) : m_isDefault(true) { }
  
  /**
   * This constructor takes a list of arguments to inititialize an 
   * amplitude and then stores them.  The user's derived class should
   * contain a similar constructor that calls this one.
   */
  Amplitude( const vector< string >& args ) : m_isDefault(false),
  m_args(args){ }
  
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
   * This method returns a string that uniquely identifies the instance
   * of the amplitude.  It is the amplitude name and all the arguments
   * concatenated together.  Every amplitude with the same identifier
   * should behave the same way.
   */
  string identifier() const;
  
  /**
   * Returns a boolean to indicate if this amplitude contains at least one
   * floating parameter.
   */
  bool containsFreeParameters() const;
  
  /**
   * If the user intendends to store intermediate calculations that are
   * static but associated with each event and each permutation of this
   * particles, then this method should be overriden with a function that
   * returns the number of user variables that will be stored.  It is 
   * recommended these are indexed with an enum.  The user must also define
   * the calcUserVars method.
   */
  virtual unsigned int numUserVars() const { return 0; }
  
  /**
   * If the user can calculate the amplitude from only the user-computed
   * data block and there is no need for the four-vectors, then the
   * user should override this function and return true.  If all amplitudes
   * in a fit can be calculated from user data then the memory consumption
   * in GPU fits can be optimizes as the raw four-vectors will not
   * be copied to the GPU.
   */
  
  virtual bool needsUserVarsOnly() const { return false; }
  
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
   * A function that indicates if this amplitude was created with the
   * default constructor.
   */
  bool isDefault() const { return ( m_isDefault == true ); }
  
  /**
   * Returns the list of arguments that was passed to the constructor.
   */
  vector<string> arguments() const { return m_args; }
  
  /**
   * The user can override this to do specific one-time tasks that
   * before the fit begins, but after the parameters have been initialized.
   * In general these tasks should go into the this init routine
   * rather than the consructor.
   */
  virtual void init(){}
  
  /**
   * This tells the AmpParameter with the indicated name to obtain its
   * value from some external location in memory.  This functionality is
   * is useful when fitting as the ParameterManager can automatically
   * updated these external values.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this amplitude, false otherwise.  The function calls
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] ptr the location in memory to find the value of the AmpParameter
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParPtr
   */
  bool setParPtr( const string& name, const double* ptr ) const;
  
  /**
   * This tells the AmpParameter with the indicated name to set its
   * value to the specified value.  The function returns true if the
   * a parameter of the indicated name was found in the list of registered
   * parameters for this amplitude, false otherwise. Function calls 
   * updatePar after setting the value.
   *
   * \param[in] name the name of the AmpParameter
   * \param[in] val the value to set AmpParameter to
   *
   * \see updatePar
   * \see AmplitudeManager::setAmpParValue
   */
  bool setParValue( const string& name, double val ) const;
  
  /**
   * The user may override this function to recalculate member data in the
   * amplitude class whenever a parameter changes.  For example, the user's
   * amplitude calculation may have expensive integration or other functions
   * that only need to be computed whenever a parameter changes, rather than
   * being computed on the fly for every event.
   *
   * \param[in] par a const reference to the parameter that has been updated
   */
  virtual void updatePar( const AmpParameter& par ) {}
  
  /**
   * \overload
   * 
   * A function used by the framework to signal that a parameter with a
   * particular name has changed.  This function searches the registered
   * parameters to see if this amplitude contains that parameter.  If so
   * it calls the virtual function updatePar with the AmpParameter.  The
   * function returns true of a parameter with the specified name is found.
   *
   * \param[in] name the name of the updated parameter
   */
  
  bool updatePar( const string& name ) const;
  
  // speed may be enhanced if the two functions below are combined
  // as this avoids an extra function call, but puts more complicated
  // calcAllAmplitudes loops in user code
  
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
   * for a set of four-vectos that describe the event kinematics.  As discussed
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
  
 
  /**
   * This loops over all events and calculates an optionally-defined
   * function specified by the user, calcUserVars, that allows the 
   * user to compute and store particular quantities associated with 
   * each event and each permutation of particles.  These may be expensive
   * quantities computed from kinematics that remain fixed for subsequent
   * amplitude calculations and hence caching them will expedite the fit.  
   *
   * \param[in] pdData a pointer to the array of data.  This is a long list
   * of GDoubles E, px, py, pz repeated sequentially for each particle in
   * the event and then the series repeated for each event in the data set.
   *
   * \param[out] pdUser data is a pointer to a block of memory to store the
   * data.  It is expected that the function will write to locations in this 
   * given by a length of iNEvents * pvPermutations.size() * numUserVars.
   *
   * \param[in] iNEvents the number of events in the data set
   *
   * \param[in] pvPermutations a pointer to the vector that contains the various
   * permutations.  The size of this vector is the number of permutations and
   * the size of one element of this vector is the number of particles.
   *
   * \see calcUserVars
   */
 
   void calcUserVarsAll( GDouble* pdData, GDouble* pdUserVars, int iNEvents,
                         const vector< vector< int > >* pvPermutations ) const;
  
  /**
   * The user should override this function in order to calculate data
   * that can be cached for each event and each permutation of particles.
   *
   * \param[in] pKin a pointer to a single event.  pKin[0][0-3] define E, px,
   * py, pz for the first particle, pKin[1][0-3] for the second, and so on
   *
   * \param[in] userVars is an optional pointer to the user data block associated
   * with this event and this permutation of particles.  It can be used to store
   * intermediate portions of the calculation in the case that calcAmplitude
   * must be called multiple times during the course of a fit.  It is strongly
   * recommended that this data block be indexed by an enum as there is
   * no ability to check that the user writes within bounds of the block.
   */
  
  virtual void calcUserVars( GDouble** pKin, GDouble* userVars ) const {}
  
  /**
   * If the calculated user data is the same for all instances of the
   * amplitude (i.e. independent of the amplitude arguments) then the
   * user can override this function and return 'true' to reduce
   * memory consumption.
   */
  virtual bool areUserVarsStatic() const { return false; }
  
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
   * Any user-defined class derived from Amplitude that has a parameter
   * in the amplitude should register the parameter using this routine.  This
   * shoudl be done for each parameter in the constructor of the user's 
   * Amplitude class.
   *
   * \param[in] par a reference to the AmpParmeter object
   */
  void registerParameter( AmpParameter& par );
  
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
  
  bool m_isDefault;
  
  vector<string> m_args;
  
  vector< AmpParameter* > m_registeredParams;
  
  mutable vector< int > m_currentPermutation;
};


#endif
