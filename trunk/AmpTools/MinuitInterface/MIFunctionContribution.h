#if !defined(MINUITINTERFACE_MIFUNCTIONCONTRIBUTION_H)
#define MINUITINTERFACE_MIFUNCTIONCONTRIBUTION_H

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

#include "MinuitInterface/MISubject.h"
#include "MinuitInterface/MIObserver.h"

class MinuitMinimizationManager;
class MinuitParameter;

class MIFunctionContribution : public MIObserver
{
public :
    
    // this flag is 
    static const double kUnknownDerivative;
    
   MIFunctionContribution( MinuitMinimizationManager* );
   virtual ~MIFunctionContribution();
   
   virtual double operator()() = 0;
   
   // function contributors can override this method if they wish to provide
   // derivatives of the function with respect to fit parameters
   virtual double derivative( const MinuitParameter& par );
   
   virtual void update( const MISubject* );
   virtual double contribution();
   
   // turn off or turn on whether this function contributes.  Check status.
   void stopContributing();
   void restartContributing();
   bool contributing() const {return m_contributing;}
   
private:
   
   MinuitMinimizationManager* m_manager;
   
   bool m_functionEvaluated;
   double m_contribution;
   bool m_contributing;
   
};

#endif
