#if !defined(MINUITINTERFACE_PARAMETERLIST_H)
#define MINUITINTERFACE_PARAMETERLIST_H

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

#include <vector>
#include <algorithm>
#include <functional>

#include "MinuitInterface/Parameter.h"
#include "MinuitInterface/MIPointerListIterator.h"
#include "MinuitInterface/MIObserver.h"

class ParameterList : public std::vector<Parameter*>
{
public:
   ParameterList() : std::vector<Parameter*>() {}
   ParameterList( int size ) : std::vector<Parameter*>(size) {}
   ParameterList( int size, Parameter*& param ) : std::vector<Parameter*>(size,param) {}
   ~ParameterList() {}
   
   typedef std::vector<Parameter*> lType;

   typedef MIPointerListIterator< lType::iterator, Parameter, Parameter > iterator;
   typedef MIPointerListIterator< lType::const_iterator, const Parameter, Parameter > const_iterator;
   
   void addParameter( Parameter& aParameter ) { push_back(&aParameter); }
   void addParameter( Parameter* aParameter ) { push_back(aParameter); }
   
   void attach( MIObserver* observer) { 
	   for( iterator par = begin(); par != end(); ++par ) par->attach( observer ); 
   }
   void detach( MIObserver* observer) { 
	   for( iterator par = begin(); par != end(); ++par ) par->detach( observer );
   }
};
#endif