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

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/AmpVecs.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "GPUManager/GPUCustomTypes.h"

#ifdef	GPU_ACCELERATION
#include "GPUManager/GPUManager.h"
#endif	//GPU_ACCELERATION

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

class AmplitudeManager
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
		
  /*
   * These functions perform computations based on the current state
   * of the AmplitudeManager.  The helper class AmpVecs is used to
   * store computations.  These functions do not change the state
   * of the AmplitudeManager.
   */
  
  /**
   * This function caculates the amplitudes for each data event and stores
   * them in the AmpVecs structure.
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure, four vectors
   * will be read from this class and amplitude caculations will be written to
   * this class
   *
   * \param[in] bIsFirstPass an optional boolean argument to aid in optimizaiton 
   * of calculations if the amplitudes have not changed.  Set to true if it
   * is the first computation on this data set, and set to false otherwise.
   *
   * \param[in] useMC an optional argument to indicate amplitudes are being
   * computed for MC.  This is only relevant when using GPU accelerated code
   * since the AmplitudeManager has two GPUManager objects, one configured for
   * computations on data and one for Monte Carlo.
   *
   * \see calcIntensities
   * \see calcSumLogIntensity
   * \see calcIntegrals
   */
  void calcAmplitudes( AmpVecs& ampVecs, bool bIsFirstPass = true, 
                       bool useMC = false ) const;
  
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
   * \param[in] bIsFirstPass an optional boolean argument to aid in optimizaiton 
   * of calculations if the amplitudes have not changed.  Set to true if it
   * is the first computation on this data set, and set to false otherwise.
   *
   * \see calcAmplitudes
   * \see calcSumLogIntensity
   * \see calcIntegrals
   */
	double calcIntensities( AmpVecs& ampVecs, bool bIsFirstPass = true ) const;

  /**
   * This function calculates the intensity for one event using a Kinematics
   * object.  This is only useful for diagnostics on a small number of events.  
   * The calcIntensities method should be used to calculate intensities
   * for many events.
   *
   * \param[in] kinematics a pointer to a Kinematics object
   */
  double calcIntensity( const Kinematics* kinematics ) const;

  /**
   * This function calculates and returns the sum of the log of the intensities
   * stored in the AmpVecs structure.  The routine will call calcIntensities
   * so the intensities do not need to be explicitly calculated first.
   *
   * \param[in,out] ampVecs a reference to the ampVecs structure from which the
   * intensities are to be read.  (This structure will be modified and updated
   * by underlying calls to calcIntensities and calcAmplitudes.)
   *
   * \param[in] bIsFirstPass an optional boolean argument to aid in optimizaiton 
   * of calculations if the amplitudes have not changed.  Set to true if it
   * is the first computation on this data set, and set to false otherwise.
   *
   * \see calcAmplitudes
   * \see calcIntensities
   * \see calcIntegrals
   */
	double calcSumLogIntensity( AmpVecs& ampVecs, bool bIsFirstPass = true ) const;
  
  /**
   * This routine calculates a square matrix with dimension equal to the number
   * of amplitudes and each element (index i,j) set to 
   * \f$ \sum_{k=1}^{N} w_k A_i A_j^* / N \f$, where N is the number of events  
   * in the ampVecs structure and \f$ w_k \f$ is the weight of each event.
   * In general the matrix is block diagonal since the elements i,j where
   * amplitudes i and j appear in different coherent sums are zero.  The returned
   * map is indexed on the string names of the amplitudes i and j.
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
  map< string, map< string, complex< double > > >
     calcIntegrals( AmpVecs& ampVecs, int iNGenEvents, 
                    bool bIsFirstPass = true ) const;

  //
  // The functions below all return information that describes the 
  // state and configuration of the AmplitudeManager.
  //
  
  /**
   * This function returns a vector of strings.  Each string corresponds to
   * the name of a unique amplitude in the amplitude manager.  The location
   * of each amplitude in the vector represents the ampIndex, an integer
   * corresponding to each amplitude that is largely used internally.
   *
   * \see ampIndex
   */
	const vector< string >& getAmpNames() const;
  
  /**
   * The function returns a reference to the parameter that is acting as the
   * scale factor for an amplitude.  Even if the user does not set a scale
   * factor, by default a fixed parameter is created and its value is
   * inititalized to 1.
   */
  const AmpParameter& getScale( const string& name ) const;
  
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
   * \see setAmpParPtr
   * \see setAmpParValue
   */
  bool hasAmpWithFreeParam() const;
  
  /**
   * This returns the internal index of an amplitude.  It is useful for users
   * who may want to index amplitude data in an array.
   *
   * \param[in] ampName the name of the amplitude to get the index of
   *
   * \see getAmpNames
   */
  int ampIndex( const string& ampName ) const;
  
  /**
   * This returns the name of the final state that was optionally passed
   * to the constructor of the AmplitudeManager.
   */
  string reactionName() const {return m_reactionName;}

  /**
   * This return returns the complex production amplitude (\f$ V_i\f$) for a
   * a single amplitude.
   *
   * \param[in] ampName the string identifying the amplitude
   *
   * \see hasProductionAmp
   * \see setExternalProductionAmplitude
   * \see setDefaultProductionAmplitude
   * \see resetProductionAmplitudes
   */
	complex< double > productionAmp( const string& ampName ) const;

	/**
   * \overload
   *
   * \param[in] ampIndex the index of the amplitude
   */
  complex< double > productionAmp( int ampIndex ) const;
  
  /**
   * This class returns a boolean indicating whether or not a production
   * amplitude exists for a particular amplitude.
   *
   * \param[in] ampName the string identifying the amplitude
   * 
   * \see productionAmp
   * \see setExternalProductionAmplitude
   * \see setDefaultProductionAmplitude
   * \see resetProductionAmplitudes
   */  
	bool hasProductionAmp( const string& ampName ) const;
  
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
   * This is used to set the default production amplitude (\f$ V_i \f$) for
   * a specific amplitude.  The production amplitude is the complex number
   * that multiplies the decay amplitude (a collection of amplitude factors)
   * in the construction of the intensity.  When fitting, this production
   * amplitude is typically a free parameter or is constrained to some other
   * free parameter in the fit.  To generate Monte Carlo, it is useful to
   * be able to set this value to any arbitrary complex number
   *
   * \param[in] ampName the name of the amplitude that is having its
   * production amplitude set
   *
   * \param[in] prodAmp the complex value of the production amplitude to set
   *
   * \see productionAmp
   * \see setExternalProductionAmplitude
   * \see resetProductionAmplitudes
   */
	void setDefaultProductionAmplitude( const string& ampName,
                                     complex< double > prodAmp );
	
  /**
   * This is used to set the memory location to which the AmplitudeManager
   * should look to find the production amplitude for a particular amplitude.
   * Note that the AmplitudeManager does not have control over this memory.
   * This relies on the user to create and preserve the memory that holds the
   * acutal values.  This is particularly useful when fitting where the
   * ParameterManager manages and updates the production amplitudes and
   * the AmplitudeManager can just look at this memory to recalculate the
   * intensity.
   *
   * \param[in] ampName the name of the amplitude to set the production 
   * parameter for
   *
   * \param[in] prodAmpPtr a pointer to memory where the AmplitudeManager
   * can look for the production amplitude
   *
   * \see productionAmp
   * \see setDefaultProductionAmplitude
   * \see resetProductionAmplitudes
   */
	void setExternalProductionAmplitude( const string& ampName,
                                       const complex< double >* prodAmpPtr );
  
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
  void setAmpParPtr( const string& ampName, const string& parName,
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
  void setAmpParValue( const string& ampName, const string& parName,
                       double ampParValue );
	
  /**
   * This loops over all amplitudes and notifies them that the parameter passed
   * in as a string has been updated.
   *
   * \param[in] parName the name of the parameter that has been updated
   *
   * \see Amplitude::updatePar
   */
  void updateAmpPar( const string& parName ) const;
    
  /**
   * This resets all of the production amplitudes to their default values.
   * If production amplitudes were referenced from external pointers, this
   * external referencing is broken by this call and the amplitudes are instead
   * derived by internally stored default values.
   */
	void resetProductionAmplitudes();

  /**
   * This will cause amplitude manager to renormalize amplitudes <b>in
   * computations of the intensity only</b>.  Each amplitude will be
   * multipled by a factor of \f$ 1 / \lang A_i(x) A_i^*(x) \rang \f$
   * that comes from the normalization integral interface.
   *
   * This scaling is useful when fitting data as it results in a trivial
   * relation between the fit parameters (production amplitudes) and
   * the number of events associated with each amplitude.
   *
   * This feature should be used with caution when coupled with constraints
   * on the production parameters.  If two amplitudes are constrained
   * to have the same production coefficient and have the same
   * scale, turning on renormalization will essentially for the two
   * amplitudes to have the same number of events in the fit.
   *
   * The default behavior is to have renormalization turned off.
   *
   * Results of calcAmplitudes or calcIntegrals will remain unchanged.
   * Note that this is necessary to avoid a circular dependency.
   *
   * \param[in] normInt the interface to provide normalization integrals
   *
   * \see disableRenormalization
   * \see ampsAreRenormalized
   */
  void renormalizeAmps( const NormIntInterface* normInt );

  /**
   * This function disables renormalization of amplitudes in the
   * calculation of the intensity.
   *
   * \see renormalizeAmps
   * \see ampsAreRenormalized
   */
  void disableRenormalization();
  
  /**
   * This function checks to see whether amplitudes are renormalized
   * in the intensity calculation.
   *
   * \see renormalizeAmps
   * \see disableRenormalization
   */
  bool ampsAreRenormalized() const { return m_renormalizeAmps; }
    
private:
	
	// recursive routine to symmetrize final state
	void generateSymmetricCombos( const vector< pair< int, int > >& prevSwaps,
                               vector< vector< pair< int, int > > > remainingSwaps,
                               const vector< int >& defaultOrder );
  
  string m_reactionName;
  
	// check to see if amplitudes have already been symmetrized so user can
	// be warned if additional amplitudes are added after symmetrization is done
	bool m_symmetrizeCalled; 
  
	// amplitude name -> vector of amplitude factors
	map< string, vector< const Amplitude* > > m_mapNameToAmps;
	
	//Vector to short-cut recomputation of amps with all fixed factors
	vector< bool > m_vbIsAmpFixed;
  
  // amplitude name -> vector of particle permutations
  // by default this starts as m_symmCombos for each amp
  map< string, vector< vector< int > > > m_ampPermutations;
	
	// amplitude name -> production amplitude
	map< string, const complex< double >* > m_prodAmp;
  
  // amplitude index -> production amplitude
  vector< const complex< double >* > m_prodAmpVec;
	
	// a vector of amplitude names -- keep also a set of amplitude indices 
  // these can be useful speeding up intesity calculations by removing
  // slower map element look-ups
	vector< string > m_ampNames;
  map< string, int > m_ampIndex;
		
  // the sum that each amplitude belongs to indexed on amplitude index
  vector< string > m_ampSum;
  
  // a 2-D vector to track whether or not the two particular amplitudes interfere
  vector< vector< bool > > m_sumCoherently;
  
	// a map to hold a set of default production amplitudes
	map< string, complex< double > > m_defaultProdAmp;
	
	// this holds "default" amplitudes for all registered amplitudes
	map< string, Amplitude* > m_registeredFactors;
	
	// a vector to hold all of the symmetric combinations of final
	// state particles
	vector< vector< int > > m_symmCombos;
    
  // a flag to track if we are renormalizing the amplitudes
  bool m_renormalizeAmps;
  const NormIntInterface* m_normInt;
  
  vector< AmpParameter > m_ampScaleVec;
    
#ifdef GPU_ACCELERATION
	
  void initGPU( const AmpVecs& a, bool useMC ) const;

  // for amplitudes with floating parameters we need a GPU manager
  // to manage both the data and the accepted Monte Carlo
  
  // cheat for now -- GPUManager class still isn't quite "const correct"
  mutable GPUManager m_dataGPUManGTX;	
  mutable GPUManager m_mcGPUManGTX;
  	
#endif //GPU_ACCELERATION
  
};

#endif
