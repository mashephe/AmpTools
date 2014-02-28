#if !defined(MINUITINTERFACE_MINUITPARAMETERMANAGER_H)
#define MINUITINTERFACE_MINUITPARAMETERMANAGER_H

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

#include <iostream>

#include <list>
#include <string>
#include <vector>

#include "MinuitInterface/MIObserver.h"

// forward declarations
class MinuitParameter;
class MinuitMinimizationManager;

using namespace std;

class MinuitParameterManager : public std::list<MinuitParameter*> 
{
   friend class MinuitParameter;
   friend class MinuitMinimizationManager;
   
public:

   // con/destructors
   MinuitParameterManager( MinuitMinimizationManager& aManager );
   ~MinuitParameterManager();

   // update parameter-related values
   void update(); // central values only
   void updateErrors();  // will only update if a fit is not in progress
      
   // synchronize minuit to current parameter states (fixed, floating, etc.)
   void synchronizeMinuit();
   
   // to be implemented
   // vector<ValuePair> scanParameter(int numberPointsInScan);
   // vector<ParameterPair> scanParameterPair( int numberPointsInScan, 
   //                                          double deltaFunction=1 );

   // could use some external matrix class... this works for now
   vector< vector< double > > covarianceMatrix();
   vector< vector< double > > correlationMatrix();

   // dump the list of parameters to an ostream.  It will generally be more
   // convenient to use the << operator defined as an inline below
   std::ostream& dump( std::ostream& s ) const;
  
   int numFloatingPars() const { return m_numFloatingPars; }
  
   // ------------ static member functions ---------------------
   
protected:
   // allow a MinuitParameter to register itself
   bool registerParameter( MinuitParameter* newParameter );
   
   // allow a MinuitParameter to remove itself from management
   void removeParameter( MinuitParameter* aParameter );
   
   // to maintain efficiency, MinuitMinimizationManager can inform when a fit is in progress
   void fitInProgress() { m_fitInProgress = true;}
   void noFitInProgress() { m_fitInProgress = false;}
   
private:
   // make default constructor private: must pass in a MinuitMinimizationManager object
   MinuitParameterManager();

   // no defaulted copy constructor or assignment operators
   MinuitParameterManager( const MinuitParameterManager& );
   
   // ----------------------- member items --------------------------
   
   MinuitMinimizationManager& m_minimizationManager; // the controlling minimization manager   
   bool m_fitInProgress;
  
   int m_numFloatingPars;
   
   // ------------ static data members -------------------------
};

inline  std::ostream& operator<<( std::ostream& aStream, 
                                  const MinuitParameterManager& aManager ) 
{
   return aManager.dump( aStream );
}
#endif
