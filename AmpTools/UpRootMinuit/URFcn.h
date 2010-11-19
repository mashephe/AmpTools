/*
 *  URFcn.h
 *  UpRootMinuit
 *
 *  Base class allowing user to specify the minuit-style function object to be
 *  invoked in a minimization procedure.   Replaces the pointer-to-function
 *  object in Root's original implementation to allow greater user flexibility
 *
 *  Created by Lawrence Gibbons on Tue Oct 21 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef UPROOT_URFcn
#define UPROOT_URFcn

#include <vector>

#include "UpRootMinuit/URtypes.h"

class URFcn
{
   
public:
   URFcn() {}
   virtual ~URFcn() {}
   
   virtual void operator()( Int_urt &npar, Double_urt *grad, Double_urt &fval, const std::vector<Double_urt>& par, Int_urt flag) = 0;
};

#endif
