//******************************************************************************
// This file is part of AmpTools, a package for performing PDFlitude Analysis
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
// EXPRESS OR IMPLIED.  By way of expdfle, but not limitation, INDIANA
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS,
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be
// held liable for any liability with respect to any claim by the user or
// any other party arising from use of the program.
//******************************************************************************

#include <iostream>
#include <sstream>
#include <fstream>

#include <string.h>

#include "IUAmpTools/PDFManager.h"
#include "IUAmpTools/NormIntInterface.h"

#ifdef VTRACE
#include "vt_user.h"
#endif

PDFManager::PDFManager( const vector< string >& reaction,
                       const string& reactionName) :
IntensityManager( reaction, reactionName ),
m_needsUserVarsOnly( true ),
m_optimizeParIteration( false )
{
  cout << endl << "## PDF MANAGER INITIALIZATION ##" << endl;
  cout << " Creating PDF manager for reaction:  " << reactionName << endl;
  
  
  // setup the vector of symmetric combinations -- first initialize
  // it with numberOfCombos copies of the default ordering
  // then go in and make the swaps
  vector< int > defaultOrder( reaction.size() );
  for( unsigned int i = 0; i < reaction.size(); ++i ){
    
    defaultOrder[i] = i;
  }

  // for a sum of PDFs the notion of "symmetrization" doesn't
  // really make sense -- in this case there is only one
  // ordering of the particles, the deafult one
    
  m_defaultOrderVec.push_back( defaultOrder );
}

PDFManager::~PDFManager() {
  
  for( map< string, PDF* >::iterator mapItr = m_registeredFactors.begin();
      mapItr != m_registeredFactors.end(); ++mapItr ){
    
    delete mapItr->second;
  }
  
  for( map< string, vector< const PDF* > >::iterator mapItr = m_mapNameToPDFs.begin();
      mapItr != m_mapNameToPDFs.end(); ++mapItr ){
    
    for( vector< const PDF* >::iterator vecItr = mapItr->second.begin();
        vecItr != mapItr->second.end(); ++vecItr ){
      
      delete *vecItr;
    }
  }
}

unsigned int
PDFManager::maxFactorStoragePerEvent() const {
  
  vector< string > pdfNames = getTermNames();
  
  unsigned int nPDFMaxFactor = 0;
  
  for( int i = 0; i < getTermNames().size(); i++ ) {
    
    unsigned int iNFactors = getFactors( pdfNames[i] ).size();
    
    if( iNFactors > nPDFMaxFactor )
      nPDFMaxFactor =  iNFactors;
  }
  
  // for each factor we need to store a real number
  return nPDFMaxFactor;
}

unsigned int
PDFManager::termStoragePerEvent() const {
  
  // for each PDF we need to store a real number
  return getTermNames().size();
}

bool
PDFManager::hasTermWithFreeParam() const {
  
  for( vector< bool >::const_iterator isFixed = m_vbIsPDFFixed.begin();
      isFixed != m_vbIsPDFFixed.end();
      ++isFixed ){
    
    if( !(*isFixed) ) return true;
  }
  
  return false;
}

bool
PDFManager::calcTerms( AmpVecs& a ) const
{
  
#ifdef VTRACE
  VT_TRACER( "PDFManager::calcTerms" );
#endif
  
  // on the first pass through this data set be sure to calculate
  // the user data first, if needed, before doing term calculations
  if( !a.m_termsValid && a.m_userVarsPerEvent > 0 ){
    
    calcUserVars( a );
    if( m_needsUserVarsOnly ) a.clearFourVecs();
  }
  
  bool modifiedTerm = false;
  
  const vector< string >& pdfNames = getTermNames();
  
  int iNPDFs = pdfNames.size();
  
  assert( iNPDFs && a.m_iNEvents && a.m_iNTrueEvents );
#ifndef GPU_ACCELERATION
  assert( a.m_pdAmps && a.m_pdAmpFactors);
#endif
  
  int iPDFIndex;
  for( iPDFIndex = 0; iPDFIndex < iNPDFs; iPDFIndex++ )
  {
    
    vector< const PDF* > vPDFs =
    m_mapNameToPDFs.find(pdfNames.at(iPDFIndex))->second;
    
    int iFactor, iNFactors = vPDFs.size();
    
    // if it is not the first pass through and this particular
    // PDF is fixed, then we can skip to the next PDF
    // and avoid all of the symmetrization computation as well
    // as checking each individual factor for free parameters.
    if( a.m_termsValid && m_vbIsPDFFixed[iPDFIndex] ) continue;
    
    // now figure out if we need to recalculate a factors for an PDF
    // in the case we have optimization turned on, we may not have changed
    // a parameter in the last iteration that is related to a factor in
    // this PDF
    
    bool recalculateFactors = false;
    
    const PDF* pCurrPDF = 0;
    for( iFactor=0; iFactor < iNFactors; iFactor++ ){
      
      pCurrPDF = vPDFs.at( iFactor );
      
      // now check to see if the value of this PDFs parameters
      // are the same as they were the last time they were evaluated
      // for this particular dataset -- if not, recalculate
      
      if( !( a.m_termsValid && m_optimizeParIteration &&
            m_dataPDFIteration[&a][pCurrPDF] == m_pdfIteration[pCurrPDF] ) ){
        
        recalculateFactors = true;
        m_dataPDFIteration[&a][pCurrPDF] = m_pdfIteration[pCurrPDF];
      }
    }
    
    if( !recalculateFactors ) continue;
    
    // if we get to here, we are changing the stored factors of the
    // PDF
    
    modifiedTerm = true;
    
    // calculate all the factors that make up an PDF for
    // for all events serially on CPU or in parallel on GPU
    int iLocalOffset = 0;
    for( iFactor=0; iFactor < iNFactors;
        iFactor++, iLocalOffset += a.m_iNEvents ){
      
      pCurrPDF = vPDFs.at( iFactor );

      // if we have static user data, look up the location in the data array
      // if not, then look up by identifier
      unsigned long long uOffset =
      ( pCurrPDF->areUserVarsStatic() ?
       a.m_staticUserVarsOffset[pCurrPDF->name()] :
       a.m_userVarsOffset[pCurrPDF->identifier()] );

#ifndef GPU_ACCELERATION
      pCurrPDF->
      calcPDFAll( a.m_pdData,
                 a.m_pdAmpFactors + iLocalOffset,
                 a.m_iNEvents, a.m_iNParticles,
                 a.m_pdUserVars + uOffset);
#else
      a.m_gpuMan.calcPDFAll( pCurrPDF, iLocalOffset,
                             uOffset );
#endif//GPU_ACCELERATION
    }
    
    
    // now assemble all the factors in an PDF into a single
    // symmetrized PDF for each event
    
#ifndef GPU_ACCELERATION
    
    int iEvent;
    unsigned long long iOffsetA, iOffsetF;
    
    // re-ordering of data will be useful to not fall out of (CPU) memory cache!!!
    
    // zeroing out the entire range
    memset( (void*)( a.m_pdAmps + a.m_iNEvents * iPDFIndex ), 0,
           a.m_iNEvents * sizeof(GDouble) );
    
    // only sum over the true events from data and skip paddings
    for( iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ )
    {
      iOffsetA = a.m_iNEvents * iPDFIndex + iEvent;
      
      a.m_pdAmps[iOffsetA] = a.m_pdAmpFactors[iEvent];
      
      for( iFactor = 1; iFactor < iNFactors; iFactor++ )
      {
        iOffsetF = iEvent + a.m_iNEvents * iFactor;
        a.m_pdAmps[iOffsetA] *= a.m_pdAmpFactors[iOffsetF];
      }
    }
    
#else
    // on the GPU the terms are assembled and never copied out
    // of GPU memory
    
    a.m_gpuMan.assembleTerms( iPDFIndex, iNFactors, iNPermutations );
#endif
    
  }
  
  a.m_termsValid = true;
  
  return modifiedTerm;
}

double
PDFManager::calcIntensities( AmpVecs& a ) const
{
  
#ifdef VTRACE
  VT_TRACER( "PDFManager::calcIntensities" );
#endif
  
  // check to be sure destination memory has been allocated
  assert( a.m_pdIntensity );
  
  double maxInten = 0;
  
  // first update the PDFs
  calcTerms( a );
  
  // In GPU running mode PDFs are maintained on the GPU and
  // the sum of the logs of intensities are calculated directly.
  // This is a CPU calculation that was likely called from the
  // parent IntensityManager with a single Kinematics object or
  // perhaps during MC generation. Copy the PDFs out of the
  // GPU and compute the intensities (on the CPU).  There is no
  // GPU accelerated intensity calculation, just a GPU accelerated
  // log( intensity ) calculation.
  
#ifdef GPU_ACCELERATION
  a.allocateCPUAmpStorage( *this );
  a.m_gpuMan.copyAmpsFromGPU( a );
#endif
  
  const vector< string >& pdfNames = getTermNames();
  
  int iNPDFs = pdfNames.size();
  
  double* coef = new double[iNPDFs];
  
  for( int i = 0; i < iNPDFs; ++i ){
    
    coef[i] = productionFactor( i ).real();
    if( termsAreRenormalized() ){

      coef[i] /= normInt()->ampInt( pdfNames[i], pdfNames[i] ).real();
    }
  }
  
  double dIntensity;
  double cAiAjRe,cAiAjIm;
  
  //Re-ordering of data will be useful to not fall out of (CPU) memory cache!!!
  //Only sum over the true events from data and skip paddings
  int iEvent;
  for( iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ )
  {
    dIntensity = 0;
    
    for( int i = 0; i < iNPDFs; i++ ){
      
      dIntensity += coef[i]*a.m_pdAmps[a.m_iNEvents*i+iEvent];
    }
    
    dIntensity *= a.m_pdWeights[iEvent];
    
    a.m_pdIntensity[iEvent] = dIntensity;
    if( dIntensity > maxInten ) maxInten = dIntensity;
  }
  
  delete[] coef;
  
  return maxInten;
}


double
PDFManager::calcSumLogIntensity( AmpVecs& a ) const
{
  
#ifdef VTRACE
  VT_TRACER( "PDFManager::calcSumLogIntensity" );
#endif
  
  // this may be inefficienct since there are two
  // loops over events, one here and one in the
  // calculation of intensities -- however, this
  // streamlines the code a little
  // this may be a place for optimization later
  
  double dSumLogI = 0;
  
#ifndef GPU_ACCELERATION
  
  calcIntensities( a );
  
  for( int iEvent=0; iEvent < a.m_iNTrueEvents; iEvent++ ){
    
    if( a.m_pdIntensity[iEvent] / a.m_pdWeights[iEvent] <= 0 ){
      
      cout << "WARNING:  intensity for event " << iEvent << " is negative.\n"
           << "          Skipping this event in likelihood calculation because\n"
           << "          one cannot take the log of a negative number." << endl;

      continue;
    }
    
    // here divide out the weight that was put into the intensity calculation
    // and weight the log -- in practice this just contributes an extra constant
    // term in the likelihood equal to sum -w_i * log( w_i ), but the division
    // helps avoid problems with negative weights, which may be used
    // in background subtraction
    dSumLogI += a.m_pdWeights[iEvent] *
    G_LOG( a.m_pdIntensity[iEvent] / a.m_pdWeights[iEvent] );
  }
  
#else
  
  // need to compute the production coefficients with all scale factors
  // taken into account
  
  vector< string > pdfNames = getTermNames();
  
  vector< complex< double > > gpuProdPars( pdfNames.size() );
  
  for( int i = 0; i < pdfNames.size(); ++i ){
    
    gpuProdPars[i] = productionFactor( pdfNames[i] );
    
    if( termsAreRenormalized() ){
      
      gpuProdPars[i] /= sqrt( normInt()->ampInt( pdfNames[i], pdfNames[i] ) );
    }
  }
  
  // need to explicitly do PDF calculation
  // since intensity and sum is done directly on GPU
  
  if( !a.m_termsValid || hasTermWithFreeParam() ){
    
    calcTerms( a );
  }
  
  dSumLogI = a.m_gpuMan.calcSumLogIntensity( gpuProdPars, m_sumCoherently );
  
#endif
  
  return( dSumLogI );
}


void
PDFManager::calcIntegrals( AmpVecs& a, int iNGenEvents ) const
{
  
#ifdef VTRACE
  VT_TRACER( "PDFManager::calcIntegrals" );
#endif
  
  GDouble* integralMatrix = a.m_pdIntegralMatrix;
  
  assert( iNGenEvents );
  bool termChanged = calcTerms( a );
  
  // if nothing changed and it isn't the first pass, return
  if( !termChanged && a.m_integralValid ) return;
  
  int iNPDFs = a.m_iNTerms;
  
  // the allocated integralMatrix object is the same size as
  // that used for Amplitude based fits where it is 2*NAmp^2
  // numbers to hold each complex-valued integral
  //
  // for a PDF fit we only need NAmp doubles, so only fill
  // the diagnonal elements
  
  int i, iEvent;
  for( i = 0; i < iNPDFs;i++ ) {
    
    
    // if the PDF isn't floating and it isn't the first pass
    // through these data, then its integral can't change
    if( a.m_integralValid && m_vbIsPDFFixed[i] ){
      
      // if the PDF isn't floating and it isn't the first pass
      // through these data, then its integral can't change
      
      continue;
    }
    else{
      
      // otherwise zero it out and recalculate it
      integralMatrix[2*i*iNPDFs+2*i] = 0;
    }
    
#ifndef GPU_ACCELERATION
    
    for( iEvent = 0; iEvent < a.m_iNTrueEvents; iEvent++ )
    {
      integralMatrix[2*i*iNPDFs+2*i] += a.m_pdWeights[iEvent] *
      a.m_pdAmps[a.m_iNEvents*i+iEvent];
    }
    // normalize
    integralMatrix[2*i*iNPDFs+2*i] /= static_cast< GDouble >( iNGenEvents );
    
#else
    a.m_gpuMan.calcIntegral( &(integralMatrix[2*i*iNPDFs+2*i]), i, i, iNGenEvents );
#endif
  }
  
  a.m_integralValid = true;
}

vector< const Term* >
PDFManager::getFactors( const string& name ) const {
  
  vector< const Term* > vTermPtr;
  
  map< string, vector< const PDF* > >::const_iterator mapItr =
  m_mapNameToPDFs.find( name );
  
  // check to be sure the PDF is there:
  assert( mapItr != m_mapNameToPDFs.end() );
  
  for( vector< const PDF* >::const_iterator vecItr = mapItr->second.begin();
      vecItr != mapItr->second.end(); ++vecItr )
    vTermPtr.push_back( *vecItr );
  
  return vTermPtr;
}

void
PDFManager::addPDFFactor( const string& name,
                         const string& factorName,
                         const vector< string >& args,
                         const string& scale ){
  
  map< string, PDF* >::iterator defaultPDF =
  m_registeredFactors.find( factorName );
  
  if( defaultPDF == m_registeredFactors.end() ){
    
    cout << "ERROR: PDF factor with name " << factorName
    << " has not been registered." << endl;
    assert( false );
  }
  
  PDF* newPDF = defaultPDF->second->newPDF( args );
  
  // check to see if this is a new term and do some setup if it is
  if( !hasTerm( name ) ){
    
    addTerm( name, scale );
    
    m_vbIsPDFFixed.push_back( true );
    m_mapNameToPDFs[name] = vector< const PDF* >( 0 );
  }
  
  m_mapNameToPDFs[name].push_back( newPDF );

  m_needsUserVarsOnly = m_needsUserVarsOnly && newPDF->needsUserVarsOnly();

  //Enable a short-cut if no factors are variable in the PDF
  m_vbIsPDFFixed[termIndex(name)] =
  m_vbIsPDFFixed[termIndex(name)] && !newPDF->containsFreeParameters();
}

void
PDFManager::setupFromConfigurationInfo( const ConfigurationInfo* configInfo ){
  
  vector< string > sumName;
  
  // loop over PDFs in the ConfigurationInfo
  vector<AmplitudeInfo*> pdfInfoVector = configInfo->amplitudeList(reactionName());
  for (unsigned int i = 0; i < pdfInfoVector.size(); i++){
    
    string pdfName = pdfInfoVector[i]->fullName();
    string scale   = pdfInfoVector[i]->scale();
    
    // add PDFs
    vector< vector<string> > pdfFactors = pdfInfoVector[i]->factors();
    for (unsigned int j = 0; j < pdfFactors.size(); j++){
      string factorName = pdfFactors[j][0];
      vector<string> pdfParameters = pdfFactors[j];
      pdfParameters.erase(pdfParameters.begin());
      addPDFFactor( pdfName, factorName, pdfParameters, scale );
    }
    
    // add production PDFs
    setDefaultProductionFactor(pdfName, pdfInfoVector[i]->value());
    
    // if the PDF has parameters we should go ahead and set their
    // values, otherwise they will be set to the default value as defined
    // by the PDFParameter class -- in a fit, the ParameterManager will
    // later reset these pointers to point to floating parameters in MINUIT
    vector< ParameterInfo* > pars = pdfInfoVector[i]->parameters();
    for( vector< ParameterInfo* >::iterator parItr = pars.begin();
        parItr != pars.end();
        ++parItr ){
      
      setParValue( pdfName, (**parItr).parName(), (**parItr).value() );
    }
    
    // finally initialize the PDFs
    vector< const PDF* > pdfVec = m_mapNameToPDFs[pdfName];
    for( vector< const PDF* >::iterator pdf = pdfVec.begin();
        pdf != pdfVec.end();
        ++pdf ) {
      
      // init needs to be non-const or else the user has to deal with
      // mutable data -- in reality we really want some aspects of the
      // PDF like the name, arguments, etc. to never change and
      // other aspects to be mutable, but this seems to put an extra
      // burden on the user
      
      const_cast< PDF* >(*pdf)->init();
    }
  }
}

void
PDFManager::setParPtr( const string& name, const string& parName,
                      const double* pdfParPtr ){
  
  IntensityManager::setParPtr( name, parName, pdfParPtr );
  
  // now look for the parameter as part of the PDF factors
  
  for( vector< const PDF* >::iterator factorItr = m_mapNameToPDFs[name].begin();
      factorItr != m_mapNameToPDFs[name].end();
      ++factorItr ){
    
    if( (**factorItr).setParPtr( parName, pdfParPtr ) ){
      
      m_vbIsPDFFixed[termIndex(name)] = false;
    }
  }
}

void
PDFManager::setParValue( const string& name, const string& parName,
                        double val ){
  
  IntensityManager::setParValue( name, parName, val );
  
  // we will redetermine the status of this variable in the loop below
  
  m_vbIsPDFFixed[termIndex(name)] = true;
  
  // now loop through the PDF factors looking for the parameter
  
  for( vector< const PDF* >::iterator factorItr = m_mapNameToPDFs[name].begin();
      factorItr != m_mapNameToPDFs[name].end();
      ++factorItr ){
    
    (**factorItr).setParValue( parName, val );
    
    m_vbIsPDFFixed[termIndex(name)] = m_vbIsPDFFixed[termIndex(name)] &&
    !(**factorItr).containsFreeParameters();
  }
}

void
PDFManager::updatePar( const string& parName ) const {
  
  for( map< string, vector< const PDF* > >::const_iterator mapItr = m_mapNameToPDFs.begin();
      mapItr != m_mapNameToPDFs.end();
      ++mapItr ){
    
    for( vector< const PDF* >::const_iterator pdfItr = mapItr->second.begin();
        pdfItr != mapItr->second.end();
        ++pdfItr ){
      
      // if we find an PDF with the parameter update the iteration
      // counter; this may result in multiple increments over one fuction
      // call but that is OK -- iteration numbers just need to be unique
      if( (**pdfItr).updatePar( parName ) ){
        
        ++m_pdfIteration[*pdfItr];
      }
    }
  }
}

void
PDFManager::registerPDFFactor( const PDF& PDF ){
  
  m_registeredFactors[PDF.name()] = PDF.clone();
}

