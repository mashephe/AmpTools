#if !defined(MINUITINTERFACE_MINUITMINIMIZATIONMANAGER_H)
#define MINUITINTERFACE_MINUITMINIMIZATIONMANAGER_H

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
#include <functional>
#include <string>

#include "MinuitInterface/MISubject.h"
#include "MinuitInterface/MinuitParameterManager.h"
#include "UpRootMinuit/URFcn.h"
#include "UpRootMinuit/URMinuit.h"

// forward declarations
class MIFunctionContribution;

class MinuitMinimizationManager : public MISubject, public URFcn
{
public:

   enum MinuitStatus { kUndefinedStatus = -1,
                       kNormal = 0,
                       kBlankCommand = 1,
                       kUnreadableCommand = 2,
                       kUnknownCommand = 3,
                       kAbnormalTermination = 4
   };
   
   enum EMatrixStatus { kNotCalculated = 0,
                        kApproxNotAccurate = 1,
                        kFullForcedPosDef = 2,
                        kFullAccurate = 3 };
  
   enum Option { kCheckDerivativeCalc, kTrustDerivativeCalc };
 
   enum FitFlag { kComputeDerivatives = 2 };
  
   enum Commands { kUnknown = 0, kMigrad = 1, kMinos = 2, kHesse = 3 };
   
   // friends
   friend class MinuitParameterManager;
   
   // con/destructors
   MinuitMinimizationManager( int maxParameters = 50 );
   ~MinuitMinimizationManager();
   
   // calculate the current value of the function to be minimized
   double evaluateFunction(); // central values
   
   // status of the last command attempted
   //   = 0: command executed normally
   //     1: command is blank, ignored
   //     2: command line unreadable, ignored
   //     3: unknown command, ignored
   //     4: abnormal termination (e.g., MIGRAD not converged)
   //     5: command is a request to read PARAMETER definitions
   //     6: 'SET INPUT' command
   //     7: 'SET TITLE' command
   //     8: 'SET COVAR' command
   //     9: reserved
   //    10: END command
   //    11: EXIT or STOP command
   //    12: RETURN command  
   int status() const;
 
   // gets a flag that indicates the last command executed
   int lastCommand() const;

   // error matrix status:
   //   0= not calculated at all
   //   1= approximation only, not accurate
   //   2= full matrix, but forced positive-definite
   //   3= full accurate covariance matrix
   int eMatrixStatus() const; 
  
   // return some information about the best minimum so far
   double bestMinimum() const;
   double estDistToMinimum() const;
  
   // turn off/on derivative calculations -- when enabled, minuit will check
   // to ensure that derivatives are being properly computed
   // if minuit rejects user-defined derivative calculation, it can be forced
   // to accept them by setting the flag
   void enableDerivatives( Option optionArg = kCheckDerivativeCalc );
   void disableDerivatives();

   bool derivativesEnabled() const { return m_derivativesEnabled; }
   
   // standard minimization procedures
   void migradMinimization();
   void minosMinimization();
   vector< vector< double > > hesseEvaluation();
   
   // change the output stream for the minuit logging
   void setLogStream( std::ostream& logStream );
   
   // change the internal precision of MINUIT (calls SET EPS)
   void setPrecision( double precision );
   double precision() const;
 
   // change the maximum number of iterations for certain MINUIT operations (default 5000)
   void setMaxIterations( int maxIter );
 
   // retrieve the setting for maximum number of iterations
   int maxIterations() const;
  
   // change the minimization strategy -- from the MINUIT manual:
   //   In the current release, this parameter can take on three integer values (0, 1, 2), 
   //   and the default value is 1. Value 0 indicates to Minuit that it should economize 
   //   function calls; it is intended for cases where there are many variable parameters 
   //   and/or the function takes a long time to calculate and/or the user is not 
   //   interested in very precise values for parameter errors. On the other hand, 
   //   the value 2 indicates that Minuit is allowed to waste function calls in order 
   //   to be sure that all values are precise; it is intended for cases where the function 
   //   is evaluated in a very short time and/or where the parameter errors must be 
   //   calculated reliably

   void setStrategy( int strategy );
   int strategy() const;
  
   // serve up the parameterManager for the user to define/delete minuit parameters
   MinuitParameterManager& parameterManager() { return m_parameterManager; }

   // but sometimes we just want to peek at the parameter manager too
   const MinuitParameterManager& parameterManager() const { return m_parameterManager; }
   
   // allow the user to specify an action everytime operator() is called
   // with a new type of flag
   void setUserFlagFunction( void (*newFlagFunction)(int) ) {
      m_newFlagFunction = newFlagFunction;
   }
   
   // the "fcn" that URMinuit needs
   void operator()( int &npar, double *grad, double &fval, const std::vector<double>& par, int flag);
   
   // ------------ static member functions ---------------------
   
protected:
   // access to minuit information
      const std::vector<double>& minuitWorkingValues() const;  // the current parameter values
      URMinuit& minuitMinimizer(); // the URMinuit object itself
   
   
private:
   // no defaulted copy constructor or assignment operators
   MinuitMinimizationManager( const MinuitMinimizationManager& );
   
   // ----------------------- member items --------------------------

   void computeDerivatives( double* grad );
   
   MinuitParameterManager m_parameterManager;
   URMinuit m_fitter;
   
   bool m_derivativesEnabled;
   void (*m_newFlagFunction)(int);

   int m_lastMinuitFlag;
  
   int m_lastCommand;
   int m_status;
                   
   int m_strategy;
  
   double m_precision;
  
   double m_bestMin;
   double m_estDistToMin;
 
   int m_eMatrixStat;
  
  int m_functionCallCounter;
};
#endif
