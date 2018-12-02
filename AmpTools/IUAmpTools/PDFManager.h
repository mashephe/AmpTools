#if !defined( PDFMANAGER )
#define PDFMANAGER

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
#include "IUAmpTools/PDF.h"
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
 * A class for managing the construction of a sum of PDFs.
 *
 * This class is used to construct the intensity for each event.
 * It manages multiple PDFs, each of which can be
 * composed of multiple PDF factors.  It also maintains a link
 * to the production parameters (typically fit parameters) stored by
 * the parameter manager.
 *
 * Tyipcally the user will setup this class once at the start of some job
 * which could be either fitting data or generating Monte Carlo.  Once
 * the class is setup, then the const member functions can be called to
 * calculate intensities for blocks of events.
 *
 * There should be one instance of the PDFManager for each reaction or
 * set of unique final state particles in the fit.
 *
 * \ingroup IUAmpTools
 */

class PDFManager : public IntensityManager
{
  
public:
  
  /** Constructor.
   * Constructs an PDFManager
   *
   * \param[in] finalState a vector of strings, one to identify each final state
   * argument.
   *
   * \param[in] reactionName an optional name for the reaction
   *
   */
  PDFManager( const vector< string >& finalState,
             const string& reactionName = "");
  
  /** Destructor.
   */
  ~PDFManager();
  
  /**
   * This returns an enum to identify the type of intensity calculation
   * being done.
   */
  Type type() const { return kPDF; }
  
  /**
   * This function returns the number of doubles needed to store all factors
   * for all PDFs for each permutation in each event.  It is the maximum
   * of number of factors * number of permutations for any term.
   */
  
  unsigned int maxFactorStoragePerEvent() const;
  
  /**
   * This function returns the number of doubles needed to store all complete
   * complex decay PDFs for each event.  It is just 2 * nAmps
   * the size of a double.
   */
  
  unsigned int termStoragePerEvent() const;
  
  /**
   * This function caculates the PDFs for each data event and stores
   * them in the AmpVecs structure.  It returns true if it alters the
   * terms in the structure (helpful to see if an update actually
   * happenend when recalculating).
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure, four vectors
   * will be read from this class and PDF caculations will be written to
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
   * algorithms.  The routine will call calcPDFs, so the PDFs
   * do not need to be explicitly calculated first.
   *
   * \param[in,out] ampVecs a reference to the AmpVecs storage structure,
   * four vectors will be read from this class and intensities written to it.
   *
   * \see calcPDFs
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
   * by underlying calls to calcIntensities and calcPDFs.)
   *
   * \see calcPDFs
   * \see calcIntensities
   * \see calcIntegrals
   */
  double calcSumLogIntensity( AmpVecs& ampVecs ) const;
  
  /**
   * This routine calculates a square matrix with dimension equal to the number
   * of PDFs and each element (index i,j) set to
   * \f$ \sum_{k=1}^{N} w_k A_i A_j^* / N \f$, where N is the number of events
   * in the ampVecs structure and \f$ w_k \f$ is the weight of each event.
   * In general the matrix is block diagonal since the elements i,j where
   * PDFs i and j appear in different coherent sums are zero.
   * The returned matrix is packed in a flat array of doubles with size 2*n^2.
   * Access to the real and imaginary parts of element i,j is by 2*i*n+2*j and
   * 2*i*n+2*j+1, respectively.
   *
   * \param[in,out] ampVecs a reference to the ampVecs structure from which the
   * PDFs are to be read.  (This structure will be modified and updated
   * by underlying calls to calcPDFs.)
   *
   * \param[in] iNGenEvents the number of genereated events (N in the above
   * computation)
   *
   * \param[in] bIsFirstPass an optional boolean argument to aid in optimizaiton
   * of calculations if the PDFs have not changed.  Set to true if it
   * is the first computation on this data set, and set to false otherwise.
   *
   * \see calcPDFs
   * \see calcIntensities
   */
  void calcIntegrals( AmpVecs& ampVecs, int iNGenEvents ) const;
  
  /**
   * This function returns a vector of const points to the PDF classes
   * that make up the factors that are multipled together to get a single
   * PDF.
   *
   * \param[in] name the name of the PDF to get the factors of
   *
   * \see addAmpFactor
   * \see registerPDFFactor
   */
  const vector< const PDF* >& getFactors( const string& name ) const;
  
  /**
   * This function returns a boolean indicating if any PDF in the
   * PDFManager contains a free parameter than might be changing
   * with each fit iteration.  It is used to trigger the recalculation of
   * normalization integrals at each fit iteration.
   *
   * \see setParPtr
   * \see setParValue
   */
  bool hasTermWithFreeParam() const;
  
  //
  // The functions below modify the state of the PDFManager
  //
  
  /**
   * This function sets up the PDFManager based on information
   * provided by a ConfigurationInfo object.  It is intended to be the
   * most user-friendly way of configuring the PDF manager for
   * use.
   *
   * \param[in] configInfo a pointer to a ConfigurationInfo object
   */
  void setupFromConfigurationInfo( const ConfigurationInfo* configInfo );
  
  /**
   * This routine will create an instance of the PDF specified in
   * factorName using the arguments specified in args and add it to the
   * list of factors that define the PDF ampName in the (optional)
   * sum.  Before addAmpFactor can be called, the specific type of
   * PDF (factorName) that user wants to add needs to be registered
   * with the PDFManager.
   *
   * \param[in] ampName the name of the PDF inside of the PDFManager
   * to which the factor is to be added.  If an PDF with this name
   * does not exist already, then it will be created.
   *
   * \param[in] factorName the type of PDF to be added.  This is typically
   * the name of a user-defined class that inherits from the PDF base class.
   * It must be registered in advance and factorName should match the return
   * value of the member function name() in the PDF class.
   *
   * \param[in] args a list of arguments to passed in the creation of the new
   * PDF factor.  This vector gets passed on to the newPDF function
   * defined in the users PDF class.
   *
   * \param[in] scale (optional) the name of a parameter to use for the
   * the scale factor for the full PDF (not a factor).  This must be
   * passed in the with the first factor of the PDF.
   * The setupFromConfigurationInfo method does this properly.
   *
   * \see getFactors
   * \see registerPDFFactor
   * \see PDF::newPDF
   */
  void addPDFFactor( const string& termName, const string& factorName,
                    const vector< string >& args,
                    const string& scale = "1.0" );
  
  /**
   * This function is used to register a user-defined PDF, which inherits
   * from the PDF base class.  PDFs must be registered before they
   * can be added as factors to specifc, named decay PDFs in the
   * PDFManager.
   *
   * \param[in] defaultPDF a reference to user-defined PDF object
   * that is simply constructed using the default construction.  For example
   * if the user PDF class is called MyAmp, the call
   * \code
   * ampManager.registerPDFFactor( MyAmp() )
   * \endcode
   * is sufficent to register the PDF.
   *
   * \see addAmpFactor
   */
  void registerPDFFactor( const PDF& defaultPDF );
  
  /**
   * This tells an PDF to use an external pointer to resolve the value
   * of some parameter inside of an PDF.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the PDF factors and check each one for the parameter
   * and then use the functionality of the PDF class to set the parameter
   * pointer.
   *
   * \param[in] ampName the name of the PDF that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParPtr a pointer to external memory that tells the PDF
   * where to find the value of the parameter
   *
   * \see setAmpParValue
   * \see PDF::setParPtr
   */
  void setParPtr( const string& termName, const string& parName,
                 const double* ampParPtr );
  
  /**
   * This tells an PDF to use a particular value for
   * some parameter inside of an PDF.  The parameter is not a
   * production paramater as discussed above, but instead a parameter that
   * may be floating in the description of the decay, e.g., the mass or width
   * of some state.  This routine is here for convenience, one could also
   * get a list of all the PDF factors and check each one for the parameter
   * and then use the functionality of the PDF class to set the parameter
   * value.
   *
   * \param[in] ampName the name of the PDF that contains a factor that
   * contains the paramter of interest
   *
   * \param[in] parName the name of the parameter
   *
   * \param[in] ampParValue the value to which to set the parameter
   *
   * \see setAmpParValue
   * \see PDF::setParValue
   */
  void setParValue( const string& termName, const string& parName,
                   double ampParValue );
  
  /**
   * This function will be called whenever a parameter is updated.
   *
   * \see PDF::updatePar
   */
  void updatePar( const string& parName ) const;
  
  /**
   * This function can be used to set a flag to optimize subsequent
   * calls to calcTerms in the case that PDFs have free parameters.
   * It uses functionality in the updatePar function to only force a
   * recalculation of the PDF for a particular AmpVecs class if
   * one of the AmpParameters for that PDF has changed since
   * the last calculation.  This should signficiantly enhance the speed
   * for fits with parameters floating in PDFs, since there will
   * only be an expensive recomputation when MINUIT changes the parameter.
   *
   * \param[in] flag set to true to enable the optimization
   */
  void setOptimizeParIteration( bool flag ) { m_optimizeParIteration = flag; }
  
private:
  
  // PDF name -> vector of PDF factors
  map< string, vector< const PDF* > > m_mapNameToPDFs;
  
  // this holds "default" PDFs for all registered PDFs
  map< string, PDF* > m_registeredFactors;
  
  // vector to short-cut recomputation of terms with all fixed factors
  vector< bool > m_vbIsPDFFixed;
  
  // some internal members to optimize PDF recalculation
  bool m_optimizeParIteration;
  mutable map< const PDF*, int > m_pdfIteration;
  mutable map< AmpVecs*, map< const PDF*, int > > m_dataPDFIteration;
};

#endif
