
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

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <ios>
#include <iomanip>
#include <cmath>

#include "MinuitInterface/MinuitParameterManager.h"
#include "MinuitInterface/MinuitParameter.h"
#include "MinuitInterface/MinuitMinimizationManager.h"
#include "MinuitInterface/MIPointerListIterator.h"

using namespace std;

typedef list<MinuitParameter*>::iterator ListIter;
typedef list<MinuitParameter*>::const_iterator ConstListIter;
typedef MIPointerListIterator<ListIter,MinuitParameter,MinuitParameter> ParamIter;

MinuitParameterManager::MinuitParameterManager( MinuitMinimizationManager& aManager ) :
   list<MinuitParameter*>(),
   m_minimizationManager( aManager ),
   m_fitInProgress( false )
{}

MinuitParameterManager::~MinuitParameterManager() 
{}

void
MinuitParameterManager::update() 
{
   for ( ParamIter parameter = begin(); parameter != end(); ++parameter ) 
   {
      int minuitId = parameter->minuitId();
      const vector<double>& minuitValues = m_minimizationManager.minuitWorkingValues();
//       cout << "Updating param " << parameter->name() 
// 	   << " with new value " << minuitValues[minuitId] << endl;
      parameter->setValue( minuitValues[minuitId] );
      parameter->invalidateErrors();
   }
}

void
MinuitParameterManager::updateErrors() 
{
   if ( m_fitInProgress || 
        m_minimizationManager.status() != MinuitMinimizationManager::kNormal) 
   { 
      return; 
   }
   
   URMinuit& theMinimizer = m_minimizationManager.minuitMinimizer();
   for ( ParamIter parameter = begin(); parameter != end(); ++parameter ) 
   {
      int minuitId = parameter->minuitId();
      double upwardsError;
      double downwardsError;
      double parabolicError;
      double correlationCoefficient;
//       cout << "querying minuit about error for param: " << parameter->name() 
//          << " MinuitId " << minuitId << endl;
      theMinimizer.mnerrs( minuitId, upwardsError, downwardsError, 
			   parabolicError, correlationCoefficient );
      parameter->setAsymmetricErrors( std::pair<double,double>(downwardsError,
							       upwardsError) );
      parameter->setGlobalCorrelation( correlationCoefficient );
      parameter->validateErrors();
      // hmm design problem?  Must have validateErrors set before we can do
      // the updateError with a notify...
      parameter->setError( parabolicError, true );
   }
}

void
MinuitParameterManager::synchronizeMinuit() 
{
  // we'll need the minuitMinimizer to check/modify parameter status
  URMinuit& theMinimizer = m_minimizationManager.minuitMinimizer();
  
  m_numFloatingPars = 0;
  
  // check that minuit has each parameter's status correct
  for ( ParamIter parameter = begin(); parameter != end(); ++parameter ) 
  {
    int minuitId = parameter->minuitId();
    
    // get the current values from minuit
    string name;
    double value; 
    double error; 
    double lowerLimit;
    double upperLimit; 
    int minuitVariableId;
    theMinimizer.mnpout( minuitId, name, value, error, 
                        lowerLimit, upperLimit, minuitVariableId );
    
    // force the minuit value to the current value of the parameter
    if ( value != parameter->value() ) 
    {
      double commandParameters[]  = { static_cast<double>(minuitId), parameter->value() };
      int status;
      theMinimizer.mnexcm( "SET PAR", commandParameters, 2, status );
    }
    
    // make sure that minuit has all of its parameters held fixed as necessary
    // Note: a positive minuitVariableId => minuit considers the parameter floating
    if ( ! parameter->floating() ) {
//      cout << "Need to fix parameter " << minuitId << " MinuitVariableId = "
//      << minuitVariableId << endl;
      if ( minuitVariableId > 0 ) {
        theMinimizer.FixParameter( minuitId );
      }
    }
    else{
      
      // we'll need this to dimension the error matrix
      ++m_numFloatingPars;
    }
    
    // check any bounds that are set on the parameter, or unbound if necessary
    bool updateLimits =  
    (lowerLimit != parameter->lowerBound()) || 
    (upperLimit != parameter->upperBound());
    if ( updateLimits ) 
    {
      int status;
      string limCommand = "SET LIM";
      if ( ! parameter->bounded() ) 
      {
        double commandParameter = minuitId;
        theMinimizer.mnexcm( limCommand, &commandParameter, 1, status );
      } 
      else 
      {
        double commandParameters[] = {static_cast<double>(minuitId), parameter->lowerBound(), 
          parameter->upperBound()};
        theMinimizer.mnexcm( limCommand, commandParameters, 3,  status );
      }
    }
    
  } // end of loop over parameters
}

bool
MinuitParameterManager::registerParameter( MinuitParameter* newParameter ) 
{

   
   // the minuit user Id for this parameter
   unsigned int newParameterMinuitId = size(); 
   ++newParameterMinuitId;
   
   // the place to insert the parameter into the list
   ListIter insertPosition = end();

   // check to see if the list is fully populated.  If not, look for the
   // first unused id, and the parameter before which the new parameter
   // will be inserted
   if ( size() && ! (back()->minuitId() == size()) ) 
   {
      newParameterMinuitId = 1;
      insertPosition = begin();
      for ( ParamIter iter = begin(); iter != end(); ++iter ) 
      {
         if ( newParameterMinuitId == iter->minuitId() ) 
	 { 
            ++newParameterMinuitId; 
            ++insertPosition;
         }  
      }
   }
   
   insert( insertPosition, newParameter );
   newParameter->setMinuitId( newParameterMinuitId );
   
   // inform minuit of the parameter, with bounds as defined 
   // in MinuitParameter
   double lowerBound = 0;
   double upperBound = 0;
   if( newParameter->bounded() )
   {
      lowerBound = newParameter->lowerBound();
      upperBound = newParameter->upperBound();
   }

   double initialStepSize = 0.1 * fabs(newParameter->value());
   // Make sure the initial step size is reasonable
   if ( initialStepSize == 0 ) 
   { 
      initialStepSize = 0.1; 
   }
   if ( ( (newParameter->value()-initialStepSize) <= lowerBound ) ||
	( (newParameter->value()+initialStepSize) >= upperBound ) )
   {
      initialStepSize *= 0.1;
   }
   URMinuit& theMinimizer = m_minimizationManager.minuitMinimizer();
   int error = theMinimizer.DefineParameter( newParameter->minuitId(), 
                                             newParameter->name(), 
                                             newParameter->value(),
                                             initialStepSize,
                                             lowerBound,
                                             upperBound );
   
   return error == 0;
}

vector< vector< double > > 
MinuitParameterManager::covarianceMatrix(){
	
	// allocate space for a square matrix of doubles with 
	// dimension equal to the number of parameters
	double* emat = (double*)malloc( m_numFloatingPars * m_numFloatingPars * 
                                  sizeof( double ) );
	
	m_minimizationManager.minuitMinimizer().mnemat( emat, m_numFloatingPars );
	
	vector< vector< double > > errorMatrix;
	
	for( unsigned int i = 0; i < m_numFloatingPars; ++i ){
		
		errorMatrix.push_back( vector< double >( 0 ) );
		for( unsigned int j = 0; j < m_numFloatingPars; ++j ){
			
			errorMatrix[i].push_back( emat[i+j*m_numFloatingPars] );
		}
	}
	
	free( emat );
	
	return errorMatrix;
}

vector< vector< double > > 
MinuitParameterManager::correlationMatrix(){


	vector< vector< double > > covMatrix = covarianceMatrix();
	vector< vector< double > > corrMatrix = covMatrix;

	for( unsigned int i = 0; i < corrMatrix.size(); ++i ){
		for( unsigned int j = 0; j < corrMatrix[i].size(); ++j ){

			corrMatrix[i][j] /= sqrt( covMatrix[i][i] * covMatrix[j][j] );
		}
	}
	
	return corrMatrix;
}

void 
MinuitParameterManager::removeParameter( MinuitParameter* aParameter ) {
   remove( aParameter );
   cout << " Warning: UpRootMinuit class does not yet support parameter removal!" << endl;
}

std::ostream& 
MinuitParameterManager::dump( std::ostream& aStream ) const 
{
 
   unsigned int maxLength = 0;
   for ( ConstListIter itParam = begin(); itParam != end(); ++itParam ) 
   {
      maxLength = (maxLength < (*itParam)->name().size()) ? 
	 (*itParam)->name().size() :  maxLength;
   }
   aStream << "\n\n" << setw(maxLength) << "name" << setw(9) 
	   << "value" << setw(9) << "error"
           << setw(9) << "-error" << setw(9) << "+error" << '\n';
   for ( ConstListIter itParam = begin(); itParam != end(); ++itParam ) 
   {
      aStream << **itParam << '\n';
   }

   return aStream;
}
