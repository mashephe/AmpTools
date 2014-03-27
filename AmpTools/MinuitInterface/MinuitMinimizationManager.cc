
// This file is a part of MinuitInterface - a front end for the Minuit minimization
//       package (Minuit itself was authored by Fred James, of CERN)
// 
// 
// Copyright Cornell University 1993, 1996, All Rights Reserved.
// 
// This software written by Lawrence Gibbons, Cornell University.
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
// obtained from Cornell University.
// 
// CORNELL MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  By way
// of example, but not limitation, CORNELL MAKES NO REPRESENTATIONS OR
// WARRANTIES OF MERCANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
// THE USE OF THIS SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS,
// COPYRIGHTS, TRADEMARKS, OR OTHER RIGHTS.  Cornell University shall not be
// held liable for any liability with respect to any claim by the user or any
// other party arising from use of the program.
//

#include <sstream>
#include "UpRootMinuit/URMinuit.h"

#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MIFunctionContribution.h"

#include <sys/time.h>

using namespace std;

MinuitMinimizationManager::MinuitMinimizationManager( int maxParameters ) :
   MISubject(),
   URFcn(),
   m_parameterManager(*this),
   m_fitter( maxParameters ),
   m_derivativesEnabled( false ),
   m_status( kUndefinedStatus ),
   m_newFlagFunction( 0 ),
   m_lastMinuitFlag( -1 ),
   m_bestMin( 0 ),
   m_estDistToMin( 0 ),
   m_eMatrixStat( 0 ),
   m_strategy( 1 ),
   m_precision( 0 ) // MINUIT determines automatically
{
//   m_parameterManager = new MinuitParameterManager( *this );
   // tell the fitter that this object is the function to evaluate
   m_fitter.SetFCN( this );
  
  setStrategy( m_strategy );
}

MinuitMinimizationManager::~MinuitMinimizationManager()
{}

double
MinuitMinimizationManager::evaluateFunction() {

   m_parameterManager.update();
   notify();

  double totalContribution = 0;
   MISubject::ObserverList& contributors = observerList();
   for ( MISubject::ObserverList::iterator iter = contributors.begin();
         iter != contributors.end();
         ++iter ) {
      MIFunctionContribution* contributor = dynamic_cast<MIFunctionContribution*>(*iter);

      // if the dynamic cast fails, i.e., the MIObserver is
      // not an MIFunctionContribution then don't ask for 
      // its contribution
     
      if( contributor != NULL )
        totalContribution += contributor->contribution();
   }
  
  ++m_functionCallCounter;
   
   return totalContribution;
}

void
MinuitMinimizationManager::computeDerivatives( double* grad )
{ 
    m_parameterManager.update();
    notify();
  
    int gradIndex = 0;
    for( MinuitParameterManager::const_iterator par = m_parameterManager.begin();
         par != m_parameterManager.end();
         ++par ){
     
        double thisDerivative = 0;
        for ( MISubject::ObserverList::iterator iter = observerList().begin();
              iter != observerList().end();
              ++iter ) {
            MIFunctionContribution* contributor = dynamic_cast<MIFunctionContribution*>(*iter);
          
          // if the dynamic cast fails, i.e., the MIObserver is
          // not an MIFunctionContribution then don't ask for 
          // its contribution

          if( contributor != NULL )
            thisDerivative += contributor->derivative( **par );
        }
        
        grad[gradIndex++] = thisDerivative;
    }
}

void 
MinuitMinimizationManager::setPrecision( double precision ){
 
  m_precision = precision;
  
  std::ostringstream command;
  command << "SET EPS " << precision;

  m_fitter.Command( command.str().c_str() );
}

double
MinuitMinimizationManager::precision() const {
  
  return m_precision;
}

void
MinuitMinimizationManager::setMaxIterations( int maxIter ){
  
  m_fitter.SetMaxIterations( maxIter );
}

int
MinuitMinimizationManager::maxIterations() const {
  
  return m_fitter.GetMaxIterations();
}

void 
MinuitMinimizationManager::setStrategy( int strat ){
  
  m_strategy = strat;
  
  std::ostringstream command;
  command << "SET STRATEGY " << strat;
  
  m_fitter.Command( command.str().c_str() );
}

int
MinuitMinimizationManager::strategy() const {
  
  return m_strategy;
}

void
MinuitMinimizationManager::enableDerivatives( Option optionArg )
{
    // set status flag first since subsequent commands
    // may indireclty rely on it
    m_derivativesEnabled = true;

    if( optionArg == kCheckDerivativeCalc ){
        
        m_fitter.Command( "SET GRADIENT" );
    }
    else{
        
        m_fitter.Command( "SET GRADIENT 1" );
    }
}

void
MinuitMinimizationManager::disableDerivatives()
{
    // set status flag first since subsequent commands
    // may indireclty rely on it
    m_derivativesEnabled = false;

    m_fitter.Command( "SET NOGRADIENT" );
}

void
MinuitMinimizationManager::migradMinimization() {
  
  m_functionCallCounter = 0;
  
  timeval tStart,tStop,tSpan;
  double dTime;
  
  gettimeofday( &(tStart), NULL );
  
   int dummyI;
   double dummyD;
   m_lastMinuitFlag = -1;
   m_lastCommand = kMigrad;
   m_parameterManager.synchronizeMinuit();
   m_parameterManager.fitInProgress();
   m_status = m_fitter.Migrad();
   m_fitter.mnstat( m_bestMin, m_estDistToMin, dummyD, dummyI, dummyI, m_eMatrixStat );
   // minuit doesn't make a final function call after evaluation of
   // errors.  This ensures that final parameter and function values
   // correspond to the minimimum found
   evaluateFunction();
   m_parameterManager.noFitInProgress();
  
  gettimeofday( &(tStop), NULL );
  timersub( &(tStop), &(tStart), &tSpan );
  dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
  
  cout << "\n MIGRAD evaluation total wall time:  " << dTime << " s." << endl;
  cout << "    average time per function call:  " << dTime * 1000 / m_functionCallCounter
  << " ms." << endl << endl;
}

void
MinuitMinimizationManager::minosMinimization() {

  m_functionCallCounter = 0;
  
  timeval tStart,tStop,tSpan;
  double dTime;
  
  gettimeofday( &(tStart), NULL );

  int dummyI;
   double dummyD;
   m_lastMinuitFlag = -1;
   m_lastCommand = kMinos;
   m_parameterManager.synchronizeMinuit();
   m_parameterManager.fitInProgress();
   m_status = m_fitter.Minos();
   m_fitter.mnstat( m_bestMin, m_estDistToMin, dummyD, dummyI, dummyI, m_eMatrixStat );
   // minuit doesn't make a final function call after evaluation of
   // errors.  This ensures that final parameter and function values
   // correspond to the minimimum found
   evaluateFunction();
   m_parameterManager.noFitInProgress();
  
  gettimeofday( &(tStop), NULL );
  timersub( &(tStop), &(tStart), &tSpan );
  dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
  
  cout << "\n MINOS evaluation total wall time:  " << dTime << " s." << endl;
  cout << "   average time per function call:  " << dTime * 1000 / m_functionCallCounter
       << " ms." << endl << endl;

}

vector< vector< double > >
MinuitMinimizationManager::hesseEvaluation() {
  
  m_functionCallCounter = 0;
  
  timeval tStart,tStop,tSpan;
  double dTime;

  gettimeofday( &(tStart), NULL );
  
   int dummyI;
   double dummyD;
   m_lastMinuitFlag = -1;
   m_lastCommand = kHesse;
   m_parameterManager.synchronizeMinuit();
   m_parameterManager.fitInProgress();
   m_status = m_fitter.Hesse();
  
  // copy out the Hessian matrix, which is filled internally,
  // and return it to the user
  int nPar = parameterManager().numFloatingPars();
  
  vector< vector< double > > secDerivMatrix;

  // first allocate space for the matrix
  for( int i = 0; i < nPar; ++i )
    secDerivMatrix.push_back( vector< double >( nPar ) );
  
  // then fill from the stored internal half matrix
  for( int i = 0; i < nPar; ++i ){
  
    for( int j = 0; j <= i; ++j ){
      
      int index = i*(i+1)/2 + j;
      secDerivMatrix[i][j] = m_fitter.fSecDer[index];
      secDerivMatrix[j][i] = m_fitter.fSecDer[index];
      
    }
  }

   m_fitter.mnstat( m_bestMin, m_estDistToMin, dummyD, dummyI, dummyI, m_eMatrixStat );
   // minuit doesn't make a final function call after evaluation of
   // errors.  This ensures that final parameter and function values
   // correspond to the minimimum found
   evaluateFunction();
   m_parameterManager.noFitInProgress();

  gettimeofday( &(tStop), NULL );
  timersub( &(tStop), &(tStart), &tSpan );
  dTime = tSpan.tv_sec + tSpan.tv_usec/1000000.0; // 10^6 uSec per second
  
  cout << "\n HESSE evaluation total wall time:  " << dTime << " s." << endl;
  cout << "   average time per function call:  " << dTime * 1000 / m_functionCallCounter
       << " ms." << endl << endl;
  
  return secDerivMatrix;
}

void
MinuitMinimizationManager::setLogStream( std::ostream& logStream ) {
   m_fitter.SetLogStream( logStream );
}

void
MinuitMinimizationManager::operator()( int &npar, double *grad, double &fval, const std::vector<double>& par, int flag) {
    
    if ( m_newFlagFunction ) {
      if ( flag != m_lastMinuitFlag ) {
         (*m_newFlagFunction)(flag);
         m_lastMinuitFlag = flag;
      }
   }
    
    if( flag == kComputeDerivatives ){
        
        // fill the grad array with the derivatives
        computeDerivatives( grad );
    }
    
   fval = evaluateFunction();
}

int
MinuitMinimizationManager::status() const {

  return m_status;
}

int
MinuitMinimizationManager::lastCommand() const {
  
  return m_lastCommand;
}

int
MinuitMinimizationManager::eMatrixStatus() const {
  
  return m_eMatrixStat;
}

double
MinuitMinimizationManager::bestMinimum() const {
  
  return m_bestMin;
}

double
MinuitMinimizationManager::estDistToMinimum() const {
  
  return m_estDistToMin;
}

const vector<double>& 
MinuitMinimizationManager::minuitWorkingValues() const {
   return m_fitter.GetParameterList();
}

URMinuit&
MinuitMinimizationManager::minuitMinimizer()  {
   return m_fitter;
}
