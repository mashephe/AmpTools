#if !defined( AMPLITUDEMANAGER )
#define AMPLITUDEMANAGER

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
#include <iostream>
#include <vector>
#include <string>
#include <complex>

#include "IUAmpTools/IntensityManager.h"
#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/AmpVecs.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "GPUManager/GPUCustomTypes.h"

#ifdef GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif //GPU_ACCELERATION

class Kinematics;
class NormIntInterface;

using namespace std;

/** 
 * A class for managing the construction of a sum of amplitudes.
 *
 * This class is used to construct the intensity for each event.
 * It manages multiple interfering amplitudes, each of which can be
 * composed of multiple amplitude factors.  It also maintains a link
 * to the production parameters (typically fit parameters) stored by
 * the parameter manager.
 *
 * Tyipcally the user will setup this class once at the start of some job
 * which could be either fitting data or generating Monte Carlo.  Once
 * the class is setup, then the const member functions can be called to
 * calculate intensities for blocks of events.
 *
 * There should be one instance of the AmplitudeManager for each reaction or
 * set of unique final state particles in the fit.
 *
 * \ingroup IUAmpTools
 */

class AmplitudeManager : public IntensityManager
{
  
public:
  
  /** Constructor.
   * Constructs an AmplitudeManager
   *
   * \param[in] finalState a vector of strings, one to identify each final state
   * argument.  Particles with identical string identifiers will be automatically
   * permuted when the intesnity is calculated.
   * 
   * \param[in] reactionName an optional name for the reaction
   *
   * \see getPermutations
   * \see addAmpPermutation
   */
  AmplitudeManager( const vector< string >& finalState, 
                    const string& reactionName = ""); 
  
  /** Destructor.
   */
  ~AmplitudeManager();
  
  /**
   * This returns an enum to identify the type of intensity calculation
   * being done.
   */
  Type type() const { return kAmplitude; }
  
  /**
   * This function returns the number of doubles needed to store all factors
   * for all amplitudes for each permutation in each event.  It is the maximum
   * of number of factors * number of permutations for any term.
   */
  
  unsigned int maxFactorStoragePerEvent() const;
  
  /**
   * This function returns the number of doubles needed to store all complete
   * complex decay amplitudes for each event.  It is just 2 * nAmps
   * the size of a double.
   */
  
  unsigned int termStoragePerEvent() const;
  
  /**
   * This function returns the number of doubles needed to store optional
   * user data for all factors and permutations.
   */
  
  unsigned int userVarsPerEvent() const;
  
  
  /**
   * This function triggers the calculation of optional user data that
   * can be stored with each amplitude to expedite future calculations.
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure, four vectors
   * will be read from this class and caculations will be written to
   * this class
   *
   */
  
  void calcUserVars( AmpVecs& ampVecs ) const;
  
  /**
   * This function caculates the amplitudes for each data event and stores
   * them in the AmpVecs structure.  It returns true if it alters the
   * terms in the structure (helpful to see if an update actually
   * happenend when recalculating).
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure, four vectors
   * will be read from this class and amplitude caculations will be written to
   * this class
   *
   * \see calcIntensities
   * \see calcSumLogIntensity
   * \see calcIntegrals
   */
  bool calcTerms( AmpVecs& ampVecs ) const;
  
  /**
   * This function calculates the intensity (one number) for all events and
   * stores it in the AmpVecs structure.  It returns the maximum intensity
   * in the data set, which is useful for accept/reject MC generation
   * algorithms.  The routine will call calcAmplitudes, so the amplitudes
   * do not need to be explicitly calculated first.
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure,
   * four vectors will be read from this class and intensities written to it.
   *
   * \see calcAmplitudes
   * \see calcSumLogIntensity
   * \see calcIntegrals
   */
  double calcIntensities( AmpVecs& ampVecs ) const;

  /**
   * This function calculates and returns the sum of the log of the intensities
   * stored in the AmpVecs structure.  The routine will call calcIntensities
   * so the intensities do not need to be explicitly calculated first.
   *
   * \param[in,out] ampVecs a reference to the ampVecs structure from which the
   * intensities are to be read.  (This structure will be modified and updated
   * by underlying calls to calcIntensities and calcAmplitudes.)
   *
   * \see calcAmplitudes
   * \see calcIntensities
   * \see calcIntegrals
   */
  double calcSumLogIntensity( AmpVecs& ampVecs ) const;
  
  /**
   * This routine calculates a square matrix with dimension equal to the number
   * of amplitudes and each element (index i,j) set to
   * \f$ \sum_{k=1}^{N} w_k A_i A_j^* / N \f$, where N is the number of events
   * in the ampVecs structure and \f$ w_k \f$ is the weight of each event.
   * In general the matrix is block diagonal since the elements i,j where
   * amplitudes i and j appear in different coherent sums are zero.  
   * The returned matrix is packed in a flat array of doubles with size 2*n^2.
   * Access to the real and imaginary parts of element i,j is by 2*i*n+2*j and
   * 2*i*n+2*j+1, respectively.
   *
   * \param[in,out] ampVecs a reference to the ampVecs structure from which the
   * amplitudes are to be read.  (This structure will be modified and updated
   * by underlying calls to calcAmplitudes.)
   *
   * \param[in] iNGenEvents the number of genereated events (N in the above
   * computation)
   *
   * \param[in] bIsFirstPass an optional boolean argument to aid in optimizaiton
   * of calculations if the amplitudes have not changed.  Set to true if it
   * is the first computation on this data set, and set to false otherwise.
   *
   * \see calcAmplitudes
   * \see calcIntensities
   */
  void calcIntegrals( AmpVecs& ampVecs, int iNGenEvents ) const;

  /**
   * The function returns a list of permutations that will be performed on
   * a particular amplitude when then amplitude is computed.  Each permutation
   * is a list of vectors that specifies how to reorder the particles in the
   * amplitude calculation.  For example, in a four-particle final state
   * the default permutation (the one that exists for every amplitude) is
   * 0 1 2 3.  If particles 2 and 3 are to be permuted in computing the
   * amplitude then the vector 0 1 3 2 will also be in the list of permutations.
   *
   * \param[in] name the name of the amplitude to get the permutations for
   *
   * \see addAmpPermutation
   */
  const vector< vector< int > >& getPermutations( const string& name ) const;
  
  /**
   * This function returns a vector of const points to the Amplitude classes
   * that make up the factors that are multipled together to get a single
   * amplitude.
   *
   * \param[in] name the name of the amplitude to get the factors of
   *
   * \see addAmpFactor
   * \see registerAmplitudeFactor
   */
  const vector< const Amplitude* >& getFactors( const string& name ) const;

  /**
   * This function returns a boolean indicating if any amplitude in the
   * AmplitudeManager contains a free parameter than might be changing
   * with each fit iteration.  It is used to trigger the recalculation of
   * normalization integrals at each fit iteration.
   *
   * \see setParPtr
   * \see setParValue
   */
  bool hasTermWithFreeParam() const;
  
  /**
   * This function will return true if every amplitude factor can be
   * calculated from user-defined data variables.  In some instances
   * this flag is used to optimize memory consumption as it means
   * that the raw four-vectors are not needed for amplitude
   * computations.
   */
  
  bool needsUserVarsOnly() const { return m_needsUserVarsOnly &&
    m_flushFourVecsIfPossible && !m_forceUserVarRecalculation; }

  //
  // The functions below modify the state of the AmplitudeManager
  //
  
  /**
   * This function sets up the AmplitudeManager based on information
   * provided by a ConfigurationInfo object.  It is intended to be the
   * most user-friendly way of configuring the amplitude manager for
   * use.
   *
   * \param[in] configInfo a pointer to a ConfigurationInfo object
   */
  void setupFromConfigurationInfo( const ConfigurationInfo* configInfo );
 
  /**
   * This routine will create an instance of the amplitude specified in
   * factorName using the arguments specified in args and add it to the
   * list of factors that define the amplitude ampName in the (optional)
   * sum.  Before addAmpFactor can be called, the specific type of
   * amplitude (factorName) that user wants to add needs to be registered
   * with the AmplitudeManager.
   *
   * \param[in] ampName the name of the amplitude inside of the AmplitudeManager
   * to which the factor is to be added.  If an amplitude with this name
   * does not exist already, then it will be created.
   *
   * \param[in] factorName the type of amplitude to be added.  This is typically
   * the name of a user-defined class that inherits from the Amplitude base class.
   * It must be registered in advance and factorName should match the return
   * value of the member function name() in the amplitude class.
   *
   * \param[in] args a list of arguments to passed in the creation of the new
   * amplitude factor.  This vector gets passed on to the newAmplitude function
   * defined in the users amplitude class.
   *
   * \param[in] sum (optional) the name of a coherent sum that the amplitude
   * should be added to (for use in cases where multiple coherent sums are
   * needed to construct the intensity).
   *
   * \param[in] scale (optional) the name of a parameter to use for the
   * the scale factor for the full amplitude (not a factor).  This must be
   * passed in the with the first factor of the amplitude.
   * The setupFromConfigurationInfo method does this properly.
   *
   * \see getFactors
   * \see registerAmplitudeFactor
   * \see Amplitude::newAmplitude
   */
  void addAmpFactor( const string& ampName, const string& factorName,
                     const vector< string >& args, const string& sum = "",
                     const string& scale = "1.0" );

  /**
   * This adds an additional permutation to any amplitude that is has not
   * already been added by the constructor.  Particles with identical names
   * are permuted automatically (this is setup in the constructor), but there
   * are other cases where one may want to permute particles for a particular
   * amplitude.  For example if one is analyzing the final state with particle
   * names "K+" "K-" "pi0", the default permutation 0 1 2 will exist for all
   * amplitudes.  It is possible that one wants to define a K* amplitude that
   * incorporates CP symmetry and therefore the permutation 1 0 2 (swapping the
   * two charged kaons) should also be added to the amplitude.  The routine
   * checks to be sure that all permutations are unique to avoid
   * double-counting.
   *
   * \param[in] ampName the name of the amplitude to which to add the permutation
   *
   * \param[in] permutation the vector of integers describing how to rearrange
   * the particle indices.  The four momentum of the particle in the i^th position
   * of the vector will be replaced by the four momentum of the particle indexed
   * by the value in the i^th position of the vector.
   *
   * \see getPermutations
   */
  void addAmpPermutation( const string& ampName, const vector< int >& permutation );

  /**
   * This function is used to register a user-defined amplitude, which inherits
   * from the amplitude base class.  Amplitudes must be registered before they
   * can be added as factors to specifc, named decay amplitudes in the
   * AmplitudeManager.
   *
   * \param[in] defaultAmplitude a reference to user-defined amplitude object
   * that is simply constructed using the default construction.  For example
   * if the user amplitude class is called MyAmp, the call
   * \code
   * ampManager.registerAmplitudeFactor( MyAmp() )
   * \endcode
   * is sufficent to register the amplitude.
   *
   * \see addAmpFactor
   */
  void registerAmplitudeFactor( const Amplitude& defaultAmplitude );
  
  /**
   * This tells an amplitude to use an external pointer to resolve the value
   * of some parameter inside of an amplitude.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the amplitude factors and check each one for the parameter
   * and then use the functionality of the Amplitude class to set the parameter
   * pointer.
   *
   * \param[in] ampName the name of the amplitude that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParPtr a pointer to external memory that tells the amplitude
   * where to find the value of the parameter
   *
   * \see setAmpParValue
   * \see Amplitude::setParPtr
   */
  void setParPtr( const string& termName, const string& parName,
                  const double* ampParPtr );
  
  /**
   * This tells an amplitude to use a particular value for
   * some parameter inside of an amplitude.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the amplitude factors and check each one for the parameter
   * and then use the functionality of the Amplitude class to set the parameter
   * value.
   *
   * \param[in] ampName the name of the amplitude that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParValue the value to which to set the parameter
   *
   * \see setAmpParValue
   * \see Amplitude::setParValue
   */
  void setParValue( const string& termName, const string& parName,
                    double ampParValue );
  
  /**
   * This function will be called whenever a parameter is updated.
   *
   * \see Amplitude::updatePar
   */
  void updatePar( const string& parName ) const;
  
  /**
   * This function can be used to set a flag to optimize subsequent
   * calls to calcTerms in the case that amplitudes have free parameters.
   * It uses functionality in the updatePar function to only force a
   * recalculation of the amplitude for a particular AmpVecs class if
   * one of the AmpParameters for that amplitude has changed since
   * the last calculation.  This should signficiantly enhance the speed
   * for fits with parameters floating in amplitudes, since there will
   * only be an expensive recomputation when MINUIT changes the parameter.
   *
   * \param[in] flag set to true to enable the optimization
   */
  void setOptimizeParIteration( bool flag ) { m_optimizeParIteration = flag; }
  
  /**
   * This function will allow the AmplitudeManager to clear four vectors
   * from memory if all amplitudes can be calculated from user variables.
   *
   * \param[in] flag set to true to enable the optimization
   */
  void setFlushFourVecsIfPossible( bool flag ) { m_flushFourVecsIfPossible = flag; }
  
  /**
   * If set to true, this flag will force the recalculation of the user
   * data every time that calcTerms is called.  This is needed in cases
   * where one reuses a common memory block multiple times with different
   * kinematics, like what happens in MC generation.
   *
   * \param[in] flag set to true to enable the optimization
   */
  void setForceUserVarRecalculation( bool flag ) {
    m_forceUserVarRecalculation = flag;
    m_flushFourVecsIfPossible = false;
  }


private:
  
  // recursive routine to symmetrize final state
  void generateSymmetricCombos( const vector< pair< int, int > >& prevSwaps,
                               vector< vector< pair< int, int > > > remainingSwaps,
                               const vector< int >& defaultOrder );
  
  // amplitude name -> vector of amplitude factors
  map< string, vector< const Amplitude* > > m_mapNameToAmps;

  // amplitude name -> vector of particle permutations
  // by default this starts as m_symmCombos for each amp
  map< string, vector< vector< int > > > m_ampPermutations;
  
  // a vector to hold all of the symmetric combinations of final
  // state particles
  vector< vector< int > > m_symmCombos;

  // check to see if amplitudes have already been symmetrized so user can
  // be warned if additional amplitudes are added after symmetrization is done
  bool m_symmetrizeCalled;

  // this holds "default" amplitudes for all registered amplitudes
  map< string, Amplitude* > m_registeredFactors;

  // the sum that each amplitude belongs to indexed on amplitude index
  vector< string > m_ampSum;
  
  // a 2-D vector to track whether or not the two particular amplitudes interfere
  vector< vector< bool > > m_sumCoherently;
  
  // vector to short-cut recomputation of terms with all fixed factors
  vector< bool > m_vbIsAmpFixed;
    
  // some internal members to optimize amplitude recalculation
  bool m_needsUserVarsOnly;
  bool m_optimizeParIteration;
  bool m_flushFourVecsIfPossible;
  bool m_forceUserVarRecalculation;
  
  mutable map< const Amplitude*, int > m_ampIteration;
  mutable map< AmpVecs*, map< const Amplitude*, int > > m_dataAmpIteration;
  mutable map< string, unsigned long long > m_staticUserVarsOffset;
};

#endif
