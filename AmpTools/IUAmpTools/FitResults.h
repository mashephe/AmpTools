#ifndef FITRESULTS
#define FITRESULTS

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


#include <vector>
#include <string>
#include <map>
#include <complex>

class LikelihoodCalculator;
class ParameterManager;
class MinuitMinimizationManager;
class IntensityManager;
class NormIntInterface;
class ConfigurationInfo;

using namespace std;

/**
 * A class for storing the information from a fit and providing an
 * interface to stored information so that subsequent analysis 
 * can be performed.
 *
 * This class captures a snapshot of the relevant information from the
 * AmplitudeManager, LikelihoodCalculator, MinuitMinimizationManager
 * and ParameterManager after a fit.  In addition it has the capability
 * to store and redeliver NormIntInterface objects to provide the
 * normalization integrals for the amplitude.  It also stores a working
 * copy of the configuration file used to perfrom the fit.  The
 * stored config file is not a verbatim copy of the original file, but
 * it functions identically.
 */

class FitResults
{
  
public:

  /**
   * Constructor for recording results of fits.
   *
   * This constructor should be used when the FitResults object is
   * going to capture the results of a fit.  It needs to have pointers
   * to the various classes that contain critical information. FitResults d
   * does not delete these pointers in its destructor so the calling function 
   * is responsible for managing the memory allocated by these objects. The
   * AmpToolsInterface class utilizes this constructor.  It is typically
   * not the constructor that a user would use.
   *
   * \param[in] cfgInfo pointer to ConfigurationInfo class
   * \param[in] intenManVec vector of IntensityManager pointers, one for each reaction
   * \param[in] likCalcMap map from reaction name to LikelihoodCalculator pointer
   * \param[in] normIntMap map from reaction name to NormIntInterface pointer
   * \param[in] parManager a pointer to the ParameterManager used in the fit
   */

  FitResults( ConfigurationInfo* cfgInfo,
              vector< IntensityManager* > intenManVec,
              map< string, LikelihoodCalculator* > likCalcMap,
              map< string, NormIntInterface* > normIntMap,
              MinuitMinimizationManager* minManager,
              ParameterManager* parManager );

  /**
   * Constructor for accessing the results of fits.
   *
   * This constructor will reconstitute a FitResults object from 
   * a file.  This is the constructor that is most useful for users
   * to open and access saved results of the fit.
   *
   * \param[in] inFile the name of the input file
   */
  FitResults( const string& inFile );
  
  /**
   * The destructor.
   *
   * The FitResults object only cleans up memory that it allocated itself.
   * It will not delete the AmplitudeManagers, LikelihoodCalculator,
   * etc. pointers that were passed in at construction time.
   */
  ~FitResults();

  /**
   * Empty constructor.
   *
   * Used only for PlotGenerator to contain histograms during
   * event generation.
   */
  FitResults() {};

  /**
   * A function to test if the FitResults object is valid.  This
   * returns true if the object was successfully constructed from
   * a file
   */
   
  bool valid() const { return m_isValid; }

  /**
   * Returns a const pointer to the NormIntInterface that is being
   * used in the fit or was reconstructed from the input file.
   *
   * \param[in] reactionName the name reaction to return the interface for
   */
  const NormIntInterface* normInt( const string& reactionName ) const;

  /**
   * Returns a const pointer to the configuration info that is being
   * used in the fit or was read in from the input file.
   */
  const ConfigurationInfo* configInfo() const { return m_cfgInfo; }
  
  /**
   * Return the global likelihood.
   */
  double likelihood() const { return m_likelihoodTotal; }

  /**
   * Return the contribution the global likelihood from a particular
   * reaction.
   *
   * \param[in] reactionName the name of the reaction
   */
  double likelihood( const string& reactionName ) const;
  
  /**
   * Return the intensity (first) and its statistical error (second)
   * for all amplitudes.  
   *
   * \param[in] accCorrected an optional boolean argument that when true
   * (default) will correct the intenisty for acceptance and when false 
   * will not perform the correction
   */
  pair< double, double > intensity( bool accCorrected = true ) const;

  /**
   * Return the intensity (first) and its statistical error (second)
   * for a subset of amplitudes.  By default this is the acceptenace corrected
   * intensity.  To obtain the product of intensity weighted by the
   * efficiency then pass in 'false' to the function.
   *
   * \param[in] amplitudes a vector of amplitudes names to sum over when
   * calculating the intensity
   *
   * \param[in] accCorrected an optional boolean argument that when true
   * (default) will correct the intenisty for acceptance and when false
   * will not perform the correction
   */
  pair< double, double > intensity( const vector< string >& amplitudes,
                                    bool accCorrected = true ) const;
  
  /**
   * Return the pahse difference (first) and its error (second) in radians
   * for two different amplitudes.
   *
   * \param[in] amp1 the first amplitude
   *
   * \param[in] amp2 the second amplitude
   */
  pair< double, double > phaseDiff( const string& amp1, const string& amp2 ) const;
  
  /**
   * Return the name of the parameter for the real value of the production  
   * coefficienct of an amplitude.
   *
   * \param[in] amplitude the amplitude name
   */
  string realProdParName( const string& amplitude ) const;
  
  /**
   * Return the name of the parameter for the imaginary value of the production
   * coefficienct of an amplitude.
   *
   * \param[in] amplitude the amplitude name
   */
  string imagProdParName( const string& amplitude ) const;
  
  /**
   * Return the name of the parameter for the scale factor for an amplitude
   *
   * \param[in] amplitude the amplitude name
   */
  string ampScaleName( const string& amplitude ) const;

  /**
   * Return the value of a parameter.
   *
   * \param[in] parName the name of the parameter
   */
  double parValue( const string& parName ) const;

  /**
   * Return the error on a parameter.
   *
   * \param[in] parName the name of the parameter
   */
  double parError( const string& parName ) const;

  /**
   * Return the covaraince of two parameters
   *
   * \param[in] par1 the first parameter
   * \param[in] par2 the second parameter
   */
  double covariance( const string& par1, const string& par2 ) const;

  /**
   * Return a map of amplitude names (strings) to values for the
   * complex production amplitudes for each amplitude.
   */
  map< string, complex< double > > ampProdParMap() const;

  /**
   * Return a map of amplitude names (strings) to values for the
   * scale factors for each amplitude.
   */
  map< string, double > ampScaleParMap() const;

  /**
   * Return a map of parameter names (strings) to the values of the
   * parameters for all amplitude parameters used in the fit.
   */
  map< string, double > ampParMap() const;
  
  /**
   * Return the value of the complex production parameter for an
   * amplitude.
   *
   * \param[in] ampName the name of the amplitude
   */
  complex< double > productionParameter( const string& ampName ) const;

  /**
   * Return the value of the complex production parameter scaled by
   * the scale factor.  This is the actual production amplitude used
   * in the fit.
   *
   * \param[in] ampName the name of the amplitude.
   */
  complex< double > scaledProductionParameter( const string& ampName ) const;
  
  /**
   * Return the real-valued scale for the amplitude.  The default for
   * this is unity, but some fits may scale amplitudes, e.g., for
   * isospin relations within a fit, or float a common scale for
   * several amplitudes.
   *
   * \param[in] ampName the name of the amplitude.
   */
  double ampScale( const string& ampName ) const;

  /**
   * Returns the status of the error matrix as reported by MINUIT
   * at the end of the of the fit.  The return value has the following
   * meanings:
   *
   *       0 = not calculated at all
   *       1 = approximation only, not accurate
   *       2 = full matrix, but forced positive-definite
   *       3 = full accurate covariance matrix
   *
   * \see MinuitMinimizationManager::EMatrixStatus
   */
  int eMatrixStatus() const { return m_eMatrixStatus; }
  
  /**
   * Returns the status of the last command executed by MINUIT.
   * It may be useful to use this to verify that minimization did not
   * terminate abnormally.  Many of the other potential return values are
   * not possible because the user does not interact with MINUIT directly.
   * The return value has the following meanings:
   *
   *       -1 = undefined status
   *        0 = normal
   *        1 = blank command
   *        2 = unreadable command
   *        3 = unkown command
   *        4 = abnormal termination (e.g., MIGRAD not converged)
   *
   * \see MinuitMinimizationManager::MinuitStatus
   */
  int lastMinuitCommandStatus() const { return m_lastCommandStatus; }

   /**
    * An integer corresponding to the last command that was executed
    * by the MinuitMinimizationManager.  The return value has the
    * following meaning:
    *
    *       0 = unkown
    *       1 = migrad
    *       2 = minos
    *       3 = hesse
    *
    * \see MinuitMinimizationManager::Commands
    */
  int lastMinuitCommand() const { return m_lastCommand; }
  
  /**
   * The value that MINUIT is assuming for the machine precision
   * when performing calculations of the likelihood.  This is typically
   * automatically determined by MINUIT, however some fits,
   * especially those that use GPU hardware, may need to specify
   * this precision so MINUIT doesn't try to minimize what are
   * fluctions in the likelihood due to machine precision.
   *
   */
  double minuitPrecision() const { return m_precision; }
  
  /**
   * The minimization strategy used by MINUIT.  See MINUIT documenation
   * for details.  Strategy ranges from 0-2 with 0 being the fastest
   * and 2 being the most robust.  Any final fits should be done
   * with strategy 2 to guarantee a proper minium.
   */
  int minuitStrategy() const { return m_strategy; }
  
  /**
   * The estimated distance to minimum as reported by MINUIT.
   */
  double estDistToMinimum() const { return m_estDistToMin; }
  
  /**
   * The value of the best minimum found by MINUIT.
   */
  double bestMinimum() const { return m_bestMin; }
  
  /**
   * A vector of all the names of all of the parameters that used in
   * the fit.  These are the actual parameters that MINUIT is 
   * manipulating.  The values of the parameters can be obtained
   * by the parValueList function and appear in the same order as
   * the names.
   * 
   * \see FitResults::parValueList
   */
  const vector< string >& parNameList() const { return m_parNames; }

  /**
   * A vector of all of the values of the parameters that are used in
   * the fit.  The order is as they appear in the parNameList.
   *
   * \see FitResults::parNameList
   */
  const vector< double >& parValueList() const { return m_parValues; }
  
  /**
   * The covariance matrix of the fit.  The order of the rows and columns
   * is as the parameters appear in the parNameList.
   *
   * \see FitResults::parNameList
   */
  const vector< vector< double > >& errorMatrix() const { return m_covMatrix; }
  
  /**
   * A vector of all of the reaction names.
   */
  const vector< string >& reactionList() const { return m_reactionNames; }

  /**
   * A single vector of all the amplitude names for all reactions.
   */
  vector< string > ampList() const;
  
  /**
   * A vector of amplitude names for a particular reaction.
   *
   * \param[in] reaction the name of the reaction
   */
  vector< string > ampList( const string& reaction ) const;
  
  /**
   * This function fetches the likelihood, parameters, MINUIT status, from
   * all of the various classes and captures this info as member data of
   * this class.  It is the function used after a fit to take a snapshot
   * of the results.
   */
  void saveResults();
  
  /**
   * Writes the current fit results object to a file.  It sould typically
   * be called after saveResults.
   *
   * \param[in] fileName the name of the output file
   */
  void writeResults( const string& fileName ) const;
  
  /**
   * Writes a file that is useful for seeding the results of subsequent fits.
   * The file contains commands to initalize the production parameters of
   * the amplitudes to the current values stored in the FitResults object.
   *
   * \param[in] fileName the name of the output file
   */
  void writeSeed( const string& fileName ) const;
  
  /**
   * Reconstitutues a FitResults object from a file.  Users will typically
   * not use this function directly, but instead use the constructor
   * that calls this function. It remains public for convenience.
   *
   * \param[in] fileName loadResults
   */
  void loadResults( const string& fileName );

  /**
   * Rotates the results of a fit to have the following phase convention:
   *
   * 1. In a loop over the amplitudes, the first one which is constrained to 
   *      be real must have a positive value
   * 2. The total amplitude must have a positive phase
   *
   */
  void rotateResults();
  void writeFitResults( const string& fileName );

protected:
  
  /**
   * Copy constructor.
   *
   * This is disabled to avoid memory handling complications since
   * FitResults may hold pointers to many objects.
   */
  FitResults( const FitResults& fitResults );

  /**
   * The assignment operator is disabled due to complications with
   * memory handling.
   */
  FitResults& operator=( const FitResults& results );

private:
  
  vector< string > stringSplit(const string& str, const string& delimiters ) const;

  void recordAmpSetup();
  void recordLikelihood();
  void recordFitStats();
  void recordParameters();
  
  // amplitude manager info
  
  int m_numReactions;
  vector< string > m_reactionNames;
  vector< int > m_numAmps;
  vector< vector< string > > m_ampNames;
  vector< vector< string > > m_ampScaleNames;
  vector< vector< double > > m_ampScaleValues;
  
  map< string, int > m_reacIndex;
  vector< map< string, int > > m_ampIndex;
  
  // likelihood info
  
  map< string, double > m_likelihoodMap;
  double m_likelihoodTotal;
  
  // info from the fitter
  
  int m_eMatrixStatus;
  int m_lastCommandStatus;
  int m_lastCommand;
  double m_precision;
  int m_strategy;
  double m_estDistToMin;
  double m_bestMin;
  
  // info about the parameters
  
  vector< string > m_parNames;
  vector< double > m_parValues;
  vector< vector< double > > m_covMatrix;
  
  map< string, int > m_parIndex;
  
  // hold pointers to the classes that are needed to fetch the results
  
  ConfigurationInfo* m_cfgInfo;
  vector< IntensityManager* > m_intenManVec;
  map< string, LikelihoodCalculator* > m_likCalcMap;
  map< string, NormIntInterface* > m_normIntMap;
  MinuitMinimizationManager* m_minManager;
  ParameterManager* m_parManager;
  
  bool m_createdFromFile;
  mutable bool m_warnedAboutFreeParams;
  bool m_isValid;
};


#endif
