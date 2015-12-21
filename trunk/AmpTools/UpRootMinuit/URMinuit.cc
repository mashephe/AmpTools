// @(#)root/minuit:$Name:  $:$Id: URMinuit.cc,v 1.5 2010/07/14 20:53:57 mashephe Exp $
// Author: Rene Brun, Frederick James   12/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//#include <stdlib.h>
//#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <functional>

#include "UpRootMinuit/URMinuit.h"
#include "UpRootMinuit/URMath.h"

//#include "TError.h"
//#include "TPluginManager.h"
//#include "TClass.h"
//#include "Api.h"

using namespace std;

typedef unsigned int uint;

const char charal[29] = " .ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const Int_urt kDefaultMaximumInternalParameters = 200; // necessary until class rewritten completely w/ STL
const Int_urt kDefaultMaximumExternalParameters = 2 * kDefaultMaximumInternalParameters;
const string kUndefinedParameterName( ")UNDEFINED" );

//______________________________________________________________________________
URMinuit::URMinuit() :
   m_userParameterIdToInternalId(kDefaultMaximumExternalParameters+1,0),
   m_userParameterValue(kDefaultMaximumExternalParameters+1,0),
   m_userParameterName(kDefaultMaximumExternalParameters+1,kUndefinedParameterName),
   m_userParameterFlag(kDefaultMaximumExternalParameters+1,-1),
   m_logStream( &cout )
{
//*-*-*-*-*-*-*-*-*-*-*Minuit normal constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ========================

   zeroPointers();
   BuildArrays(kDefaultMaximumInternalParameters);

   fStatus     = 0;
   fEmpty      = 0;
   SetMaxIterations();
   mninit(5,6,7);
}

//______________________________________________________________________________
URMinuit::URMinuit(Int_urt maxpar) :
   m_userParameterIdToInternalId(2*maxpar+1),
   m_userParameterValue(2*maxpar+1,0),
   m_userParameterName(2*maxpar+1,kUndefinedParameterName),
   m_userParameterFlag(2*maxpar+1,-1),
   m_logStream( &cout )
{
//*-*-*-*-*-*-*-*-*-*-*Minuit normal constructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ========================
//
//  maxpar is the maximum number of parameters used with this URMinuit object.

   zeroPointers();
   BuildArrays(maxpar);

   fStatus     = 0;
   fEmpty      = 0;
   SetMaxIterations();

   mninit(5,6,7);
}

//______________________________________________________________________________
URMinuit::~URMinuit()
{
//*-*-*-*-*-*-*-*-*-*-*Minuit default destructor*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  =========================

   DeleteArrays();
}

//______________________________________________________________________________
void
URMinuit::zeroPointers()
{
   //preset all pointers to null
//   fCpnam     = 0;
//   fU         = 0;
   fAlim      = 0; 
   fBlim      = 0; 
   fPstar     = 0;
   fGin       = 0;
//   fNvarl     = 0; 
//   fNiofex    = 0; 
   
   fNexofi    = 0; 
   fIpfix     = 0; 
   fErp       = 0; 
   fErn       = 0; 
   fWerr      = 0; 
   fGlobcc    = 0; 
   fX         = 0;
   fXt        = 0; 
   fDirin     = 0; 
   fXs        = 0; 
   fXts       = 0; 
   fDirins    = 0; 
   fGrd       = 0; 
   fG2        = 0; 
   fGstep     = 0; 
   fDgrd      = 0; 
   fGrds      = 0;
   fG2s       = 0; 
   fGsteps    = 0; 
   fPstst     = 0; 
   fPbar      = 0; 
   fPrho      = 0; 
   fWord7     = 0; 
   fVhmat     = 0;
   fSecDer    = 0;
   fVthmat    = 0;
   fP         = 0; 
   fXpt       = 0; 
   fYpt       = 0; 
   fChpt      = 0; 
   fCONTgcc   = 0; 
   fCONTw     = 0; 
   fFIXPyy    = 0; 
   fGRADgf    = 0; 
   fHESSyy    = 0; 
   fIMPRdsav  = 0; 
   fIMPRy     = 0; 
   fMATUvline = 0; 
   fMIGRflnu  = 0; 
   fMIGRstep  = 0; 
   fMIGRgs    = 0; 
   fMIGRvg    = 0; 
   fMIGRxxs   = 0; 
   fMNOTxdev  = 0; 
   fMNOTw     = 0; 
   fMNOTgcc   = 0; 
   fPSDFs     = 0; 
   fSEEKxmid  = 0; 
   fSEEKxbest = 0; 
   fSIMPy     = 0; 
   fVERTq     = 0; 
   fVERTs     = 0; 
   fVERTpp    = 0; 
   fCOMDplist = 0; 
   fPARSplist = 0;   
}

//______________________________________________________________________________
void URMinuit::BuildArrays(Int_urt maxpar)
{
//*-*-*-*-*-*-*Create internal Minuit arrays for the maxpar parameters*-*-*-*
//*-*          =======================================================

//   fMaxpar = kDefaultMaximumInternalParameters;
//   fMaxpar2= 2*fMaxpar;
//   if (maxpar > kDefaultMaximumInternalParameters) {
      fMaxpar = maxpar;
      fMaxpar2= 2*fMaxpar;
//      m_userParameterIdToInternalId.resize( fMaxpar2+1, 0 );
//      m_userParameterValue.resize( fMaxpar2+1, 0 );
//      m_userParameterName.resize( fMaxpar2+1, kUndefinedParameterName );
//      m_userParameterFlag.resize( fMaxpar2+1, -1 );
//   }
   fMaxpar1= fMaxpar*(fMaxpar+1);
   fMaxpar5= fMaxpar1/2;
   fMaxcpt = 101;
//   fCpnam  = new string[fMaxpar2];
//   fU      = new Double_urt[fMaxpar2];
   fAlim   = new Double_urt[fMaxpar2];
   fBlim   = new Double_urt[fMaxpar2];
   fPstar  = new Double_urt[fMaxpar2];
   fGin    = new Double_urt[fMaxpar2];
//   fNvarl  = new Int_urt[fMaxpar2];
//   fNiofex = new Int_urt[fMaxpar2];
   
   fNexofi = new Int_urt[fMaxpar];
   fIpfix  = new Int_urt[fMaxpar];
   fErp    = new Double_urt[fMaxpar];
   fErn    = new Double_urt[fMaxpar];
   fWerr   = new Double_urt[fMaxpar];
   fGlobcc = new Double_urt[fMaxpar];
   fX      = new Double_urt[fMaxpar];
   fXt     = new Double_urt[fMaxpar];
   fDirin  = new Double_urt[fMaxpar];
   fXs     = new Double_urt[fMaxpar];
   fXts    = new Double_urt[fMaxpar];
   fDirins = new Double_urt[fMaxpar];
   fGrd    = new Double_urt[fMaxpar];
   fG2     = new Double_urt[fMaxpar];
   fGstep  = new Double_urt[fMaxpar];
   fDgrd   = new Double_urt[fMaxpar];
   fGrds   = new Double_urt[fMaxpar];
   fG2s    = new Double_urt[fMaxpar];
   fGsteps = new Double_urt[fMaxpar];
   fPstst  = new Double_urt[fMaxpar];
   fPbar   = new Double_urt[fMaxpar];
   fPrho   = new Double_urt[fMaxpar];
   fWord7  = new Double_urt[fMaxpar];
   fVhmat  = new Double_urt[fMaxpar5];
   fSecDer = new Double_urt[fMaxpar5];
   fVthmat = new Double_urt[fMaxpar5];
   fP      = new Double_urt[fMaxpar1];
   fXpt    = new Double_urt[fMaxcpt];
   fYpt    = new Double_urt[fMaxcpt];
   fChpt   = new char[fMaxcpt+1];
   // initialisation of dynamic arrays used internally in some functions
   // these arrays had a fix dimension in Minuit
   fCONTgcc   = new Double_urt[fMaxpar];
   fCONTw     = new Double_urt[fMaxpar];
   fFIXPyy    = new Double_urt[fMaxpar];
   fGRADgf    = new Double_urt[fMaxpar];
   fHESSyy    = new Double_urt[fMaxpar];
   fIMPRdsav  = new Double_urt[fMaxpar];
   fIMPRy     = new Double_urt[fMaxpar];
   fMATUvline = new Double_urt[fMaxpar];
   fMIGRflnu  = new Double_urt[fMaxpar];
   fMIGRstep  = new Double_urt[fMaxpar];
   fMIGRgs    = new Double_urt[fMaxpar];
   fMIGRvg    = new Double_urt[fMaxpar];
   fMIGRxxs   = new Double_urt[fMaxpar];
   fMNOTxdev  = new Double_urt[fMaxpar];
   fMNOTw     = new Double_urt[fMaxpar];
   fMNOTgcc   = new Double_urt[fMaxpar];
   fPSDFs     = new Double_urt[fMaxpar];
   fSEEKxmid  = new Double_urt[fMaxpar];
   fSEEKxbest = new Double_urt[fMaxpar];
   fSIMPy     = new Double_urt[fMaxpar];
   fVERTq     = new Double_urt[fMaxpar];
   fVERTs     = new Double_urt[fMaxpar];
   fVERTpp    = new Double_urt[fMaxpar];
   fCOMDplist = new Double_urt[fMaxpar];
   fPARSplist = new Double_urt[fMaxpar];

   for (int i = 0; i < fMaxpar; i++) {
      fErp[i] = 0;
      fErn[i] = 0;
   }
}


//______________________________________________________________________________
Int_urt URMinuit::Command(const char *command)
{
// execute a Minuit command
//     Equivalent to MNEXCM except that the command is given as a
//     character string.
// See URMinuit::mnhelp for the full list of available commands
//
// Returns the status of the execution:
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
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   Int_urt status = 0;
   mncomd(command,status);
   return status;
}


//______________________________________________________________________________
Int_urt URMinuit::DefineParameter( Int_urt parNo, const string& name, Double_urt initVal, Double_urt initErr, Double_urt lowerLimit, Double_urt upperLimit )
{
// Define a parameter

     Int_urt err;

     const string& sname = name;
     mnparm( parNo, sname, initVal, initErr, lowerLimit, upperLimit, err);

     return err;
}

//______________________________________________________________________________
void URMinuit::DeleteArrays()
{
//*-*-*-*-*-*-*-*-*-*-*-*Delete internal Minuit arrays*-*-*-*-*-*-*-*-*
//*-*                    =============================
   if (fEmpty) return;
//   delete [] fCpnam;
//   delete [] fU;
   delete [] fAlim;
   delete [] fBlim;
   delete [] fErp;
   delete [] fErn;
   delete [] fWerr;
   delete [] fGlobcc;
//   delete [] fNvarl;
//   delete [] fNiofex;
   delete [] fNexofi;
   delete [] fX;
   delete [] fXt;
   delete [] fDirin;
   delete [] fXs;
   delete [] fXts;
   delete [] fDirins;
   delete [] fGrd;
   delete [] fG2;
   delete [] fGstep;
   delete [] fGin;
   delete [] fDgrd;
   delete [] fGrds;
   delete [] fG2s;
   delete [] fGsteps;
   delete [] fIpfix;
   delete [] fVhmat;
   delete [] fVthmat;
   delete [] fP;
   delete [] fPstar;
   delete [] fPstst;
   delete [] fPbar;
   delete [] fPrho;
   delete [] fWord7;
   delete [] fXpt;
   delete [] fYpt;
   delete [] fChpt;
   
   delete [] fCONTgcc;
   delete [] fCONTw;
   delete [] fFIXPyy;
   delete [] fGRADgf;
   delete [] fHESSyy;
   delete [] fIMPRdsav;
   delete [] fIMPRy;
   delete [] fMATUvline;
   delete [] fMIGRflnu;
   delete [] fMIGRstep;
   delete [] fMIGRgs;
   delete [] fMIGRvg;
   delete [] fMIGRxxs;
   delete [] fMNOTxdev;
   delete [] fMNOTw;
   delete [] fMNOTgcc;
   delete [] fPSDFs;
   delete [] fSEEKxmid;
   delete [] fSEEKxbest;
   delete [] fSIMPy;
   delete [] fVERTq;
   delete [] fVERTs;
   delete [] fVERTpp;
   delete [] fCOMDplist;
   delete [] fPARSplist;
   
   fEmpty = 1;
}

//______________________________________________________________________________
Int_urt URMinuit::Eval(Int_urt npar, Double_urt *grad, Double_urt &fval, const vector<Double_urt>& par, Int_urt flag)
{
// Evaluate the minimisation function
//  Input parameters:
//    npar:    number of currently variable parameters
//    par:     array of (constant and variable) parameters
//    flag:    Indicates what is to be calculated (see example below)
//    grad:    array of gradients
//  Output parameters:
//    fval:    The calculated function value. 
//    grad:    The (optional) vector of first derivatives).
// 
// The meaning of the parameters par is of course defined by the user, 
// who uses the values of those parameters to calculate his function value. 
// The starting values must be specified by the user.
// Later values are determined by Minuit as it searches for the minimum 
// or performs whatever analysis is requested by the user.
//
// Note that this virtual function may be redefined in a class derived from URMinuit.
// The default function calls the function specified in SetFCN
//
// Example of Minimisation function:
/*
   if (flag == 1) {
      read input data, 
      calculate any necessary constants, etc. 
   }
   if (flag == 2) {
      calculate GRAD, the first derivatives of FVAL 
     (this is optional) 
   }
   Always calculate the value of the function, FVAL, 
   which is usually a chisquare or log likelihood. 
   if (iflag == 3) {
      will come here only after the fit is finished. 
      Perform any final calculations, output fitted data, etc. 
   }
*/
//  See concrete examples in TH1::H1FitChisquare, H1FitLikelihood
   
   if (fFCN) (*fFCN)(npar,grad,fval,par,flag);
   return 0;
}

//______________________________________________________________________________
bool URMinuit::parameterFixed( Int_urt externalParameterId )
{

   // obtain internal Minuit id
   Int_urt kint = m_userParameterIdToInternalId[externalParameterId];
   
   //*-*-            see if parameter is floating
   if (kint > 0) {
      return false;
   }
   return true;
}

//______________________________________________________________________________
Int_urt URMinuit::FixParameter( Int_urt parNo)
{
// fix a parameter

     Int_urt err;
//     Double_urt tmp = parNo+1; //set internal Minuit nuumbering
     Double_urt tmp = parNo; 

     mnexcm( "FIX", &tmp,  1,  err );

     return err;
}


//______________________________________________________________________________
Int_urt URMinuit::GetParameter( Int_urt parNo, Double_urt &currentValue, Double_urt &currentError ) const
{
// return parameter value and error
     Int_urt    err;
     string   name; // ignored
     Double_urt bnd1, bnd2; // ignored

     mnpout( parNo, name, currentValue, currentError, bnd1, bnd2, err );

     return err;
}

//______________________________________________________________________________
Int_urt URMinuit::GetNumFixedPars() const
{
// returns the number of currently fixed parameters

   return fNpfix;
}

//______________________________________________________________________________
Int_urt URMinuit::GetNumFreePars() const
{
// returns the number of currently free parameters

   return fNpar;
}

//______________________________________________________________________________
Int_urt URMinuit::GetNumPars() const
{
// returns the total number of parameters that have been defined.
// (fixed and free)

   return fNpar + fNpfix;
}

//______________________________________________________________________________
Int_urt URMinuit::Migrad()
{
// invokes the MIGRAD minimizer
     Int_urt err;
  
// configure command arguments here
     Double_urt plist[] = { static_cast<Double_urt>(fMaxIterations) };
     Int_urt npar = 1;

     mnexcm( "MIGRAD", plist, npar, err );

     return err;
}

//______________________________________________________________________________
Int_urt URMinuit::Hesse()
{
   // invokes the HESSE error evaluation
   Int_urt err;

  // configure command arguments here
   Double_urt plist[] = { static_cast<Double_urt>(fMaxIterations) };
   Int_urt npar = 1;
  
   mnexcm( "HESSE", plist, npar, err );
   
   return err;
}

//______________________________________________________________________________
Int_urt URMinuit::Minos()
{
   // invokes the MINOS minimizer
   Int_urt err;
  
  // configure command arguments here
   Double_urt plist[] = { static_cast<Double_urt>(fMaxIterations) };
   Int_urt npar = 1;
   
   mnexcm( "MINOS", plist, npar, err );
   
   return err;
}

//______________________________________________________________________________
Int_urt URMinuit::Release( Int_urt parNo)
{
// release a parameter

     Int_urt err;
//     Double_urt tmp = parNo+1; //set internal Minuit numbering
     Double_urt tmp = parNo; //set internal Minuit numbering

     mnexcm( "RELEASE", &tmp, 1, err );

     return err;
}

//______________________________________________________________________________
Int_urt URMinuit::SetErrorDef( Double_urt up )
{
     Int_urt err;

     mnexcm( "SET ERRDEF", &up, 1, err );

     return err;
}

//______________________________________________________________________________
void URMinuit::SetFCN( URFcn* fcn )
{
//*-*-*-*-*-*-*To set the address of the minimization function*-*-*-*-*-*-*-*
//*-*          ===============================================
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   fFCN = fcn;
}

//______________________________________________________________________________
Int_urt URMinuit::SetPrintLevel( Int_urt printLevel )
{
     Int_urt    err;
     Double_urt tmp = printLevel;

     mnexcm( "SET PRINT", &tmp, 1, err );

     return err;
}

//______________________________________________________________________________
void URMinuit::mnamin()
{
//*-*-*-*-*-*-*-*-*-*-*-*-*Initialize AMIN*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                      ===============
//*-*C        Called  from many places.  Initializes the value of AMIN by
//*-*C        calling the user function. Prints out the function value and
//*-*C        parameter values if Print Flag value is high enough.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt fnew;
    Int_urt nparx;

    nparx = fNpar;
    if (fISW[4] >= 1) {
       *m_logStream << " FIRST CALL TO USER FUNCTION AT NEW START POINT, WITH IFLAG=4." << endl;;
    }
    mnexin(fX);
//    Eval(nparx, fGin, fnew, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, fnew, m_userParameterValue, 4);    ++fNfcn;
    fAmin = fnew;
    fEDM  = fBigedm;
} /* mnamin_ */

//______________________________________________________________________________
void URMinuit::mnbins(Double_urt a1, Double_urt a2, Int_urt naa, Double_urt &bl, Double_urt &bh, Int_urt &nb, Double_urt &bwid)
{
//*-*-*-*-*-*-*-*-*-*-*Compute reasonable histogram intervals*-*-*-*-*-*-*-*-*
//*-*                  ======================================
//*-*        Function TO DETERMINE REASONABLE HISTOGRAM INTERVALS
//*-*        GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2
//*-*        AND DESIRED MAXIMUM NUMBER OF BINS NAA
//*-*        PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID
//*-*        F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt awid,ah, al, sigfig, sigrnd, alb;
    Int_urt kwid, lwid, na=0, log_;

    al = URMath::Min(a1,a2);
    ah = URMath::Max(a1,a2);
    if (al == ah) ah = al + 1;

//*-*-       IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE
    if (naa == -1) goto L150;
L10:
    na = naa - 1;
    if (na < 1) na = 1;

//*-*-        GET NOMINAL BIN WIDTH IN EXPON FORM
L20:
    awid = (ah-al) / Double_urt(na);
    log_ = Int_urt(URMath::Log10(awid));
    if (awid <= 1) --log_;
    sigfig = awid*URMath::Power(10, -log_);
//*-*-       ROUND MANTISSA UP TO 2, 2.5, 5, OR 10
    if (sigfig > 2) goto L40;
    sigrnd = 2;
    goto L100;
L40:
    if (sigfig > 2.5) goto L50;
    sigrnd = 2.5;
    goto L100;
L50:
    if (sigfig > 5) goto L60;
    sigrnd = 5;
    goto L100;
L60:
    sigrnd = 1;
    ++log_;
L100:
    bwid = sigrnd*URMath::Power(10, log_);
    goto L200;
//*-*-       GET NEW BOUNDS FROM NEW WIDTH BWID
L150:
    if (bwid <= 0) goto L10;
L200:
    alb  = al / bwid;
    lwid = Int_urt(alb);
    if (alb < 0) --lwid;
    bl   = bwid*Double_urt(lwid);
    alb  = ah / bwid + 1;
    kwid = Int_urt(alb);
    if (alb < 0) --kwid;
    bh = bwid*Double_urt(kwid);
    nb = kwid - lwid;
    if (naa > 5) goto L240;
    if (naa == -1) return;
//*-*-        REQUEST FOR ONE BIN IS DIFFICULT CASE
    if (naa > 1 || nb == 1) return;
    bwid *= 2;
    nb = 1;
    return;
L240:
    if (nb << 1 != naa) return;
    ++na;
    goto L20;
} /* mnbins_ */

//______________________________________________________________________________
void URMinuit::mncalf(Double_urt *pvec, Double_urt &ycalf)
{
//*-*-*-*-*-*-*-*-*-*Transform FCN to find further minima*-*-*-*-*-*-*-*-*-*
//*-*                ====================================
//*-*        Called only from MNIMPR.  Transforms the function FCN
//*-*        by dividing out the quadratic part in order to find further
//*-*        minima.    Calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Int_urt ndex, i, j, m, n, nparx;
    Double_urt denom, f;

    nparx = fNpar;
    mninex(&pvec[0]);
//    Eval(nparx, fGin, f, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, f, m_userParameterValue, 4);    ++fNfcn;
    for (i = 1; i <= fNpar; ++i) {
	fGrd[i-1] = 0;
	for (j = 1; j <= fNpar; ++j) {
	    m = URMath::Max(i,j);
	    n = URMath::Min(i,j);
	    ndex = m*(m-1) / 2 + n;
	    fGrd[i-1] += fVthmat[ndex-1]*(fXt[j-1] - pvec[j-1]);
	}
    }
    denom = 0;
    for (i = 1; i <= fNpar; ++i) {denom += fGrd[i-1]*(fXt[i-1] - pvec[i-1]); }
    if (denom <= 0) {
	fDcovar = 1;
	fISW[1] = 0;
	denom   = 1;
    }
    ycalf = (f - fApsi) / denom;
} /* mncalf_ */

//______________________________________________________________________________
void URMinuit::mncler()
{
//*-*-*-*-*-*-*-*-*-*-*Resets the parameter list to UNDEFINED*-*-*-*-*-*-*-*
//*-*                  ======================================
//*-*        Called from MINUIT and by option from MNEXCM
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt i;

    fNpfix = 0;
    fNu = 0;
    fNpar = 0;
    fNfcn = 0;
    fNwrmes[0] = 0;
    fNwrmes[1] = 0;
    for (i = 1; i <= fMaxext; ++i) {
//	fU[i-1]      = 0;
       m_userParameterValue[i] = 0;
//	fCpnam[i-1]  = fCundef;
       m_userParameterName[i] = fCundef;
//	fNvarl[i-1]  = -1;
       m_userParameterFlag[i] = -1;
//	fNiofex[i-1] = 0;
       m_userParameterIdToInternalId[i] = 0;
    }
    mnrset(1);
    fCfrom  = "CLEAR   ";
    fNfcnfr = fNfcn;
    fCstatu = "UNDEFINED ";
    fLnolim = kurTRUE;
    fLphead = kurTRUE;
} /* mncler_ */

//______________________________________________________________________________
void URMinuit::mncntr(Int_urt ike1, Int_urt ike2, Int_urt &ierrf)
{
//*-*-*-*-*Print function contours in two variables, on line printer*-*-*-*-*
//*-*      =========================================================
//*-*
//*-*                input arguments: parx, pary, devs, ngrid
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    static string clabel = "0123456789ABCDEFGHIJ";

    /* Local variables */
    Double_urt d__1, d__2;
    Double_urt fcna[115], fcnb[115], contur[20];
    Double_urt  ylabel, fmn, fmx, xlo, ylo, xup, yup;
    Double_urt devs, xsav, ysav,  bwidx,  bwidy, unext, ff, xb4;
    Int_urt i,  ngrid, ixmid, nparx, ix, nx, ny, ki1, ki2, ixzero, iy, ics;
    string chmid, chln, chzero;

    // modified ke1, ke2 to not increment ike1,ike2
    Int_urt ke1 = ike1;
    Int_urt ke2 = ike2;
    if (ke1 <= 0 || ke2 <= 0) goto L1350;
    if (ke1 > fNu || ke2 > fNu) goto L1350;
//    ki1 = fNiofex[ke1-1];
//    ki2 = fNiofex[ke2-1];
    ki1 = m_userParameterIdToInternalId[ke1];
    ki2 = m_userParameterIdToInternalId[ke2];
    if (ki1 <= 0 || ki2 <= 0) goto L1350;
    if (ki1 == ki2) goto L1350;

    if (fISW[1] < 1) {
	mnhess();
	mnwerr();
    }
    nparx = fNpar;
//    xsav = fU[ke1-1];
//    ysav = fU[ke2-1];
    xsav = m_userParameterValue[ke1];
    ysav = m_userParameterValue[ke2];
    devs = fWord7[2];
    if (devs <= 0) devs = 2;
//    xlo = fU[ke1-1] - devs*fWerr[ki1-1];
//    xup = fU[ke1-1] + devs*fWerr[ki1-1];
//    ylo = fU[ke2-1] - devs*fWerr[ki2-1];
//    yup = fU[ke2-1] + devs*fWerr[ki2-1];
    xlo = m_userParameterValue[ke1] - devs*fWerr[ki1-1];
    xup = m_userParameterValue[ke1] + devs*fWerr[ki1-1];
    ylo = m_userParameterValue[ke2] - devs*fWerr[ki2-1];
    yup = m_userParameterValue[ke2] + devs*fWerr[ki2-1];
    ngrid = Int_urt(fWord7[3]);
    if (ngrid <= 0) {
	ngrid = 25;
//*-*  Computing MIN
	nx = URMath::Min(fNpagwd - 15,ngrid);
//*-*  Computing MIN
	ny = URMath::Min(fNpagln - 7,ngrid);
    } else {
	nx = ngrid;
	ny = ngrid;
    }
    if (nx < 11)   nx = 11;
    if (ny < 11)   ny = 11;
    if (nx >= 115) nx = 114;

//*-*-        ask if parameter outside limits
//    if (fNvarl[ke1-1] > 1) {
    if (m_userParameterFlag[ke1] > 1) {
	if (xlo < fAlim[ke1-1]) xlo = fAlim[ke1-1];
	if (xup > fBlim[ke1-1]) xup = fBlim[ke1-1];
    }
//    if (fNvarl[ke2-1] > 1) {
    if (m_userParameterFlag[ke2] > 1) {
	if (ylo < fAlim[ke2-1]) ylo = fAlim[ke2-1];
	if (yup > fBlim[ke2-1]) yup = fBlim[ke2-1];
    }
    bwidx = (xup - xlo) / Double_urt(nx);
    bwidy = (yup - ylo) / Double_urt(ny);
    ixmid = Int_urt(((xsav - xlo)*Double_urt(nx) / (xup - xlo)) + 1);
    if (fAmin == fUndefi) mnamin();

    for (i = 1; i <= 20; ++i) {	contur[i-1] = fAmin + fUp*(i-1)*(i-1); }
    contur[0] += fUp*.01;
//*-*-               fill FCNB to prepare first row, and find column zero/
//    fU[ke2-1] = yup;
    m_userParameterValue[ke2] = yup;
    ixzero = 0;
    xb4 = 1;
//TH
    chmid.resize(nx+1, ' ');
    chzero.resize(nx+1, ' ');
    chln.resize(nx+1, ' ');
    for (ix = 1; ix <= nx + 1; ++ix) {
//       fU[ke1-1] = xlo + Double_urt(ix-1)*bwidx;
	m_userParameterValue[ke1] = xlo + Double_urt(ix-1)*bwidx;
//       Eval(nparx, fGin, ff, fU, 4);
	Eval(nparx, fGin, ff, m_userParameterValue, 4);
	fcnb[ix-1] = ff;
//        if (xb4 < 0 && fU[ke1-1] > 0) ixzero = ix - 1;
	if (xb4 < 0 && m_userParameterValue[ke1] > 0) ixzero = ix - 1;
//        xb4          = fU[ke1-1];
	xb4          = m_userParameterValue[ke1];
	chmid[ix-1]  = '*';
	chzero[ix-1] = '-';
    }
//    std::printf(" Y-AXIS: PARAMETER %3d: %s\n",ke2, fCpnam[ke2-1].c_str() );
//    std::printf(" Y-AXIS: PARAMETER %3d: %s\n",ke2, m_userParameterName[ke2].c_str() );
    *m_logStream << " Y-AXIS: PARAMETER " << setw(3) << ke2 << ": " << m_userParameterName[ke2] << '\n';
    if (ixzero > 0) {
	chzero[ixzero-1] = '+';
	chln = " ";
        *m_logStream << "             X=0\n";
	//std::printf("             X=0\n");
    }
//*-*-                loop over rows
    for (iy = 1; iy <= ny; ++iy) {
//       unext = fU[ke2-1] - bwidy;
	unext = m_userParameterValue[ke2] - bwidy;
//*-*-                prepare this line background pattern for contour
	chln = " ";
// TH	
	chln.resize(nx+1, ' ');
	chln[ixmid-1] = '*';
	if (ixzero != 0) chln[ixzero-1] = ':';
//        if (fU[ke2-1] > ysav && unext < ysav) chln = chmid;
//        if (fU[ke2-1] > 0 && unext < 0)       chln = chzero;
	if (m_userParameterValue[ke2] > ysav && unext < ysav) chln = chmid;
	if (m_userParameterValue[ke2] > 0 && unext < 0)       chln = chzero;
//        fU[ke2-1] = unext;
//        ylabel = fU[ke2-1] + bwidy*.5;
	m_userParameterValue[ke2] = unext;
	ylabel = m_userParameterValue[ke2] + bwidy*.5;
//*-*-                move FCNB to FCNA and fill FCNB with next row
	for (ix = 1; ix <= nx + 1; ++ix) {
	    fcna[ix-1] = fcnb[ix-1];
//           fU[ke1-1] = xlo + Double_urt(ix-1)*bwidx;
	    m_userParameterValue[ke1] = xlo + Double_urt(ix-1)*bwidx;
//            Eval(nparx, fGin, ff, fU, 4);
	    Eval(nparx, fGin, ff, m_userParameterValue, 4);
	    fcnb[ix-1] = ff;
	}
//*-*-                look for contours crossing the FCNxy squares
	for (ix = 1; ix <= nx; ++ix) {
	    d__1 = URMath::Max(fcna[ix-1],fcnb[ix-1]),
	    d__2 = URMath::Max(fcna[ix],fcnb[ix]);
	    fmx  = URMath::Max(d__1,d__2);
	    d__1 = URMath::Min(fcna[ix-1],fcnb[ix-1]),
            d__2 = URMath::Min(fcna[ix],fcnb[ix]);
	    fmn  = URMath::Min(d__1,d__2);
	    for (ics = 1; ics <= 20; ++ics) {
		if (contur[ics-1] > fmn)  goto L240;
	    }
	    continue;
L240:
	    if (contur[ics-1] < fmx) chln[ix-1] = clabel[ics-1];
	}
//*-*-                print a row of the contour plot
        *m_logStream << setw(12) << setprecision(4) << ylabel << chln << '\n';
        // std::printf(" %12.4g %s\n",ylabel,chln.c_str() );
    }
//*-*-                contours printed, label x-axis
    chln            = " ";
    chln.resize(nx+1, ' ');
    chln.replace(0,1,"I");
    chln.replace(ixmid-1,1,"I");
    chln.replace(nx-1,1,"I");
    // std::printf("              %s\n", chln.c_str());
    *m_logStream << "              " << chln << '\n';

//*-*-               the hardest of all: print x-axis scale!
    chln =  " ";
    if (nx <= 26) {
        //std::printf("        %12.4g%s%12.4g\n",xlo,chln.c_str(),xup);
       *m_logStream << "        " << setw(12) << setprecision(4) << xlo << chln << setw(12) << xup << '\n';
        //std::printf("              %s%12.4g\n",chln.c_str(),xsav);
       *m_logStream << "              " << chln << setw(12) << setprecision(4) << xsav << '\n';
    } else {
        //std::printf("        %12.4g%s%12.4g%s%12.4g\n",xlo,chln.c_str(),xsav,chln.c_str(),xup);
       *m_logStream << "        " << setw(12) << setprecision(4) << xlo << chln 
            << setw(12) << setprecision(4) << xsav << chln
            << setw(12) << setprecision(4) << xup << '\n';
    }
//    std::printf("       X-AXIS: PARAMETER%3d: %s  ONE COLUMN=%12.4g\n"
//            ,ke1,(fCpnam[ke1-1]).c_str(),bwidx);
    //std::printf("       X-AXIS: PARAMETER%3d: %s  ONE COLUMN=%12.4g\n"
    //        ,ke1,(m_userParameterName[ke1]).c_str(),bwidx);
    *m_logStream << "       X-AXIS: PARAMETER" << setw(3) << ke1 << ": "<< m_userParameterName[ke1]
         << "  ONE COLUMN=" << setw(12) << setprecision(4) << bwidx << '\n';
    //std::printf(" FUNCTION VALUES: F(I)=%12.4g +%12.4g *I**2\n",fAmin,fUp);
    *m_logStream << " FUNCTION VALUES: F(I)=" << setw(12) << setprecision(4) << fAmin
         << " +" << setw(12) << setprecision(4) << fUp << " *I**2\n";
//*-*-                finished.  reset input values
//    fU[ke1-1] = xsav;
//    fU[ke2-1] = ysav;
    m_userParameterValue[ke1] = xsav;
    m_userParameterValue[ke2] = ysav;
    ierrf     = 0;
    *m_logStream << endl;
    return;
L1350:
       *m_logStream << " INVALID PARAMETER NUMBER(S) REQUESTED.  IGNORED." << endl;
    //std::printf(" INVALID PARAMETER NUMBER(S) REQUESTED.  IGNORED.\n");
    ierrf = 1;
} /* mncntr_ */

//______________________________________________________________________________
void URMinuit::mncomd(const string& crdbin, Int_urt &icondn)
{
//*-*-*-*-*-*-*-*-*-*-*Reads a command string and executes*-*-*-*-*-*-*-*-*-*
//*-*                  ===================================
//*-*        Called by user.  'Reads' a command string and executes.
//*-*     Equivalent to MNEXCM except that the command is given as a
//*-*          character string.
//*-*
//*-*     ICONDN = 0: command executed normally
//*-*              1: command is blank, ignored
//*-*              2: command line unreadable, ignored
//*-*              3: unknown command, ignored
//*-*              4: abnormal termination (e.g., MIGRAD not converged)
//*-*              5: command is a request to read PARAMETER definitions
//*-*              6: 'SET INPUT' command
//*-*              7: 'SET TITLE' command
//*-*              8: 'SET COVAR' command
//*-*              9: reserved
//*-*             10: END command
//*-*             11: EXIT or STOP command
//*-*             12: RETURN command
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Int_urt ierr, ipos, i, llist, lenbuf, lnc;
    Bool_urt leader;
    string comand, crdbuf, ctemp;

    crdbuf = crdbin;
    int (*pf)(int)=toupper; 
    transform(crdbuf.begin(), crdbuf.end(), crdbuf.begin(), pf); 
    lenbuf = crdbuf.length();
    icondn = 0;
//*-*-    record not case-sensitive, get upper case, strip leading blanks
    leader = kurTRUE;
    ipos = 1;
    for (i = 1; i <= URMath::Min(20,lenbuf); ++i) {
	if (crdbuf[i-1] == '\'') break;
	if (crdbuf[i-1] == ' ') {
	    if (leader) ++ipos;
	    continue;
	}
	leader = kurFALSE;
    }

//*-*-                    blank or null command
    if (ipos > lenbuf) {
        *m_logStream << " BLANK COMMAND IGNORED." << endl;
	//std::printf(" BLANK COMMAND IGNORED.\n");
	icondn = 1;
	return;
    }
//*-*-                                          . .   preemptive commands
//*-*-              if command is 'PARAMETER'
    if (crdbuf.substr(ipos-1,3) == "PAR") {
	icondn = 5;
	fLphead = kurTRUE;
	return;
    }
//*-*-              if command is 'SET INPUT'
    if (crdbuf.substr(ipos-1,7) == "SET INP") {
	icondn = 6;
	fLphead = kurTRUE;
	return;
    }
//*-*-              if command is 'SET TITLE'
    if (crdbuf.substr(ipos-1,7) == "SET TIT") {
	icondn = 7;
	fLphead = kurTRUE;
	return;
    }
//*-*-              if command is 'SET COVARIANCE'
    if (crdbuf.substr(ipos-1,7) == "SET COV") {
	icondn = 8;
	fLphead = kurTRUE;
	return;
    }
//*-*-              crack the command . . . . . . . . . . . . . . . .
    ctemp = crdbuf.substr(ipos-1,lenbuf-ipos+1);
    mncrck(ctemp, 20, comand, lnc, fMaxpar, fCOMDplist, llist, ierr, fIsyswr);
    if (ierr > 0) {
	//std::printf(" COMMAND CANNOT BE INTERPRETED\n");
        *m_logStream << " COMMAND CANNOT BE INTERPRETED" << endl;
	icondn = 2;
	return;
    }

    mnexcm(comand, fCOMDplist, llist, ierr);
    icondn = ierr;
} /* mncomd_ */

//______________________________________________________________________________
void URMinuit::mncont(Int_urt ike1, Int_urt ike2, Int_urt nptu, Double_urt *xptu, Double_urt *yptu, Int_urt &ierrf)
{
//*-*-*-*-*-*-*Find points along a contour where FCN is minimum*-*-*-*-*-*-*
//*-*          ================================================
//*-*       Find NPTU points along a contour where the function
//*-*             FMIN (X(KE1),X(KE2)) =  AMIN+UP
//*-*       where FMIN is the minimum of FCN with respect to all
//*-*       the other NPAR-2 variable parameters (if any).
//*-*   IERRF on return will be equal to the number of points found:
//*-*     NPTU if normal termination with NPTU points found
//*-*     -1   if errors in the calling sequence (KE1, KE2 not variable)
//*-*      0   if less than four points can be found (using MNMNOT)
//*-*     n>3  if only n points can be found (n < NPTU)
//*-*
//*-*                 input arguments: parx, pary, devs, ngrid
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    /* System generated locals */
    Int_urt i__1;

    /* Local variables */
    Double_urt d__1, d__2;
    Double_urt dist, xdir, ydir, aopt,  u1min, u2min;
    Double_urt abest, scalx, scaly;
    Double_urt a1, a2, val2mi, val2pl, dc, sclfac, bigdis, sigsav;
    Int_urt nall, iold, line, mpar, ierr, inew, move, next, i, j, nfcol, iercr;
    Int_urt idist=0, npcol, kints, i2, i1, lr, nfcnco=0, ki1, ki2, ki3, ke3;
    Int_urt nowpts, istrav, nfmxin, isw2, isw4;
    Bool_urt ldebug;

    /* Function Body */
//    Int_urt ke1 = ike1+1;
//    Int_urt ke2 = ike2+1;
    Int_urt ke1 = ike1;
    Int_urt ke2 = ike2;
    ldebug = fIdbg[6] >= 1;
//    if (ke1 <= 0 || ke2 <= 0) goto L1350;
//    if (ke1 > fNu || ke2 > fNu) goto L1350;
    if (ke1 <= 0 || ke2 <= 0) goto L1350;
    if (ke1 > fNu || ke2 > fNu) goto L1350;
//    ki1 = fNiofex[ke1-1];
//    ki2 = fNiofex[ke2-1];
    ki1 = m_userParameterIdToInternalId[ke1];
    ki2 = m_userParameterIdToInternalId[ke2];
    if (ki1 <= 0 || ki2 <= 0) goto L1350;
    if (ki1 == ki2) goto L1350;
    if (nptu < 4)  goto L1400;

    nfcnco  = fNfcn;
    fNfcnmx = (nptu + 5)*100*(fNpar + 1);
//*-*-          The minimum
    mncuve();
//    u1min  = fU[ke1-1];
//    u2min  = fU[ke2-1];
    u1min  = m_userParameterValue[ke1];
    u2min  = m_userParameterValue[ke2];
    ierrf  = 0;
    fCfrom = "MNContour ";
    fNfcnfr = nfcnco;
    if (fISW[4] >= 0) {
        //std::printf(" START MNCONTOUR CALCULATION OF%4d POINTS ON CONTOUR.\n",nptu);
        *m_logStream << " START MNCONTOUR CALCULATION OF" << setw(4) << nptu << " POINTS ON CONTOUR." << endl;
	if (fNpar > 2) {
	    if (fNpar == 3) {
		ki3 = 6 - ki1 - ki2;
		ke3 = fNexofi[ki3-1];
//                std::printf(" EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER %3d  %s\n",ke3,fCpnam[ke3-1].c_str());
		//std::printf(" EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER %3d  %s\n",ke3,m_userParameterName[ke3].c_str());
                *m_logStream << " EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER " << setw(3) << ke3
                     << "  " << m_userParameterName[ke3] << endl;
	    } else {
		//std::printf(" EACH POINT IS A MINIMUM WITH RESPECT TO THE OTHER%3d VARIABLE PARAMETERS.\n",fNpar - 2);
                *m_logStream << " EACH POINT IS A MINIMUM WITH RESPECT TO THE OTHER" << setw(3) << (fNpar - 2)
                   << " VARIABLE PARAMETERS." << endl;
	    }
	}
    }

//*-*-          Find the first four points using MNMNOT
//*-*-             ........................ first two points
    mnmnot(ke1, ke2, val2pl, val2mi);
    if (fErn[ki1-1] == fUndefi) {
	xptu[0] = fAlim[ke1-1];
	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");
    } else {
	if (fErn[ki1-1] >= 0) goto L1500;
	xptu[0] = u1min + fErn[ki1-1];
    }
    yptu[0] = val2mi;

    if (fErp[ki1-1] == fUndefi) {
	xptu[2] = fBlim[ke1-1];
	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");
    } else {
	if (fErp[ki1-1] <= 0) goto L1500;
	xptu[2] = u1min + fErp[ki1-1];
    }
    yptu[2] = val2pl;
    scalx = 1 / (xptu[2] - xptu[0]);
//*-*-             ........................... next two points
    mnmnot(ke2, ke1, val2pl, val2mi);
    if (fErn[ki2-1] == fUndefi) {
	yptu[1] = fAlim[ke2-1];
	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");
    } else {
	if (fErn[ki2-1] >= 0) goto L1500;
	yptu[1] = u2min + fErn[ki2-1];
    }
    xptu[1] = val2mi;
    if (fErp[ki2-1] == fUndefi) {
	yptu[3] = fBlim[ke2-1];
	mnwarn("W", "MNContour ", "Contour squeezed by parameter limits.");
    } else {
	if (fErp[ki2-1] <= 0) goto L1500;
	yptu[3] = u2min + fErp[ki2-1];
    }
    xptu[3] = val2pl;
    scaly   = 1 / (yptu[3] - yptu[1]);
    nowpts  = 4;
    next    = 5;
    if (ldebug) {
	//std::printf(" Plot of four points found by MINOS\n");
        *m_logStream << " Plot of four points found by MINOS" << endl;
	fXpt[0]  = u1min;
	fYpt[0]  = u2min;
	fChpt[0] = ' ';
//*-*  Computing MIN
	nall = URMath::Min(nowpts + 1,101);
	for (i = 2; i <= nall; ++i) {
	    fXpt[i-1] = xptu[i-2];
	    fYpt[i-1] = yptu[i-2];
	}
	//std::printf(fChpt,"%s"," ABCD");
        *m_logStream << fChpt << " ABCD" << endl;
	mnplot(fXpt, fYpt, fChpt, nall, fNpagwd, fNpagln);
    }

//*-*-              ..................... save some values before fixing
    isw2   = fISW[1];
    isw4   = fISW[3];
    sigsav = fEDM;
    istrav = fIstrat;
    dc     = fDcovar;
    fApsi  = fEpsi*.5;
    abest  = fAmin;
    mpar   = fNpar;
    nfmxin = fNfcnmx;
    for (i = 1; i <= mpar; ++i) { fXt[i-1] = fX[i-1]; }
    i__1 = mpar*(mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) { fVthmat[j-1] = fVhmat[j-1]; }
    for (i = 1; i <= mpar; ++i) {
	fCONTgcc[i-1] = fGlobcc[i-1];
	fCONTw[i-1]   = fWerr[i-1];
    }
//*-*-                          fix the two parameters in question
//    kints = fNiofex[ke1-1];
    kints = m_userParameterIdToInternalId[ke1];
    mnfixp(kints-1, ierr);
//    kints = fNiofex[ke2-1];
    kints = m_userParameterIdToInternalId[ke2];
    mnfixp(kints-1, ierr);
//*-*-              ......................Fill in the rest of the points
    for (inew = next; inew <= nptu; ++inew) {
//*-*            find the two neighbouring points with largest separation
	bigdis = 0;
	for (iold = 1; iold <= inew - 1; ++iold) {
	    i2 = iold + 1;
	    if (i2 == inew) i2 = 1;
	    d__1 = scalx*(xptu[iold-1] - xptu[i2-1]);
	    d__2 = scaly*(yptu[iold-1] - yptu[i2-1]);
	    dist = d__1*d__1 + d__2*d__2;
	    if (dist > bigdis) {
		bigdis = dist;
		idist  = iold;
	    }
	}
	i1 = idist;
	i2 = i1 + 1;
	if (i2 == inew) i2 = 1;
//*-*-                  next point goes between I1 and I2
	a1 = .5;
	a2 = .5;
L300:
	fXmidcr = a1*xptu[i1-1] + a2*xptu[i2-1];
	fYmidcr = a1*yptu[i1-1] + a2*yptu[i2-1];
	xdir    = yptu[i2-1] - yptu[i1-1];
	ydir    = xptu[i1-1] - xptu[i2-1];
	sclfac  = URMath::Max(URMath::Abs(xdir*scalx),URMath::Abs(ydir*scaly));
	fXdircr = xdir / sclfac;
	fYdircr = ydir / sclfac;
	fKe1cr  = ke1;
	fKe2cr  = ke2;
//*-*-               Find the contour crossing point along DIR
	fAmin = abest;
	mncros(aopt, iercr);
	if (iercr > 1) {
//*-*-             If cannot find mid-point, try closer to point 1
	    if (a1 > .5) {
		if (fISW[4] >= 0) {
		    //std::printf(" MNCONT CANNOT FIND NEXT POINT ON CONTOUR.  ONLY%3d POINTS FOUND.\n",nowpts);
                    *m_logStream << " MNCONT CANNOT FIND NEXT POINT ON CONTOUR.  ONLY" << setw(3) << nowpts
                                 << " POINTS FOUND." << endl;
		}
		goto L950;
	    }
	    mnwarn("W", "MNContour ", "Cannot find midpoint, try closer.");
	    a1 = .75;
	    a2 = .25;
	    goto L300;
	}
//*-*-               Contour has been located, insert new point in list
	for (move = nowpts; move >= i1 + 1; --move) {
	    xptu[move] = xptu[move-1];
	    yptu[move] = yptu[move-1];
	}
	++nowpts;
	xptu[i1] = fXmidcr + fXdircr*aopt;
	yptu[i1] = fYmidcr + fYdircr*aopt;
    }
L950:

    ierrf = nowpts;
    fCstatu = "SUCCESSFUL";
    if (nowpts < nptu) 	fCstatu = "INCOMPLETE";

//*-*-               make a lineprinter plot of the contour
    if (fISW[4] >= 0) {
	fXpt[0]  = u1min;
	fYpt[0]  = u2min;
	fChpt[0] = ' ';
	nall = URMath::Min(nowpts + 1,101);
	for (i = 2; i <= nall; ++i) {
	    fXpt[i-1]  = xptu[i-2];
	    fYpt[i-1]  = yptu[i-2];
	    fChpt[i-1] = 'X';
	}
        fChpt[nall] = 0;
//      std::printf(" Y-AXIS: PARAMETER %3d  %s\n",ke2,fCpnam[ke2-1].c_str());
        //std::printf(" Y-AXIS: PARAMETER %3d  %s\n",ke2,m_userParameterName[ke2].c_str());
	*m_logStream << " Y-AXIS: PARAMETER " << setw(3) << ke2 << "  " << m_userParameterName[ke2] << '\n';

	mnplot(fXpt, fYpt, fChpt, nall, fNpagwd, fNpagln);

//      std::printf("                         X-AXIS: PARAMETER %3d  %s\n",ke1,fCpnam[ke1-1].c_str());
	//std::printf("                         X-AXIS: PARAMETER %3d  %s\n",ke1,m_userParameterName[ke1].c_str());
        *m_logStream << "                         X-AXIS: PARAMETER " << setw(3) << ke1 << "  " << m_userParameterName[ke1] << '\n';
    }
//*-*-                print out the coordinates around the contour
    if (fISW[4] >= 1) {
	npcol = (nowpts + 1) / 2;
	nfcol = nowpts / 2;
	//std::printf("%5d POINTS ON CONTOUR.   FMIN=%13.5e   ERRDEF=%11.3g\n",nowpts,abest,fUp);
	*m_logStream << setw(5) << nowpts << " POINTS ON CONTOUR.   FMIN=" << setw(13) << setprecision(5) << abest
	     << "   ERRDEF=" << setw(11) << setprecision(3) << fUp << '\n';
//      std::printf("         %s%s%s%s\n",fCpnam[ke1-1].c_str(),
//            fCpnam[ke2-1].c_str(),
//            fCpnam[ke1-1].c_str(),
//            fCpnam[ke2-1].c_str());
        //std::printf("         %s%s%s%s\n",m_userParameterName[ke1].c_str(),
        //                           m_userParameterName[ke2].c_str(),
        //                           m_userParameterName[ke1].c_str(),
        //                           m_userParameterName[ke2].c_str());
	*m_logStream << "         " << m_userParameterName[ke1] << m_userParameterName[ke2]
	                    << m_userParameterName[ke1] << m_userParameterName[ke1] << '\n';
	for (line = 1; line <= nfcol; ++line) {
	    lr = line + npcol;
	    //std::printf(" %5d%13.5e%13.5e          %5d%13.5e%13.5e\n",line,xptu[line-1],yptu[line-1],lr,xptu[lr-1],yptu[lr-1]);
	    *m_logStream << setw(6) << line << setw(13) << setprecision(5) << xptu[line-1] << setw(13) << setprecision(5) << yptu[line-1]
		 << setw(15) << lr  << setw(13) << setprecision(5) << xptu[lr-1]   << setw(13) << setprecision(5) << yptu[lr-1] << '\n';
	}
	if (nfcol < npcol) {
	  //std::printf(" %5d%13.5e%13.5e\n",npcol,xptu[npcol-1],yptu[npcol-1]);
	    *m_logStream << setw(6) << npcol << setw(13) << setprecision(5) << xptu[npcol-1]   << setw(13) << setprecision(5) << yptu[npcol-1] << '\n';
	}
    }
//*-*-                                   . . contour finished. reset v
    fItaur = 1;
    mnfree(1);
    mnfree(1);
    i__1 = mpar*(mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) { fVhmat[j-1] = fVthmat[j-1]; }
    for (i = 1; i <= mpar; ++i) {
	fGlobcc[i-1] = fCONTgcc[i-1];
	fWerr[i-1]   = fCONTw[i-1];
	fX[i-1]      = fXt[i-1];
    }
    mninex(fX);
    fEDM    = sigsav;
    fAmin   = abest;
    fISW[1] = isw2;
    fISW[3] = isw4;
    fDcovar = dc;
    fItaur  = 0;
    fNfcnmx = nfmxin;
    fIstrat = istrav;
//    fU[ke1-1] = u1min;
//    fU[ke2-1] = u2min;
    m_userParameterValue[ke1] = u1min;
    m_userParameterValue[ke2] = u2min;
    goto L2000;
//*-*-                                    Error returns
L1350:
    //std::printf(" INVALID PARAMETER NUMBERS.\n");
    *m_logStream << " INVALID PARAMETER NUMBERS.\n";
    goto L1450;
L1400:
    //std::printf(" LESS THAN FOUR POINTS REQUESTED.\n");
    *m_logStream << " LESS THAN FOUR POINTS REQUESTED.\n";
L1450:
    ierrf   = -1;
    fCstatu = "USER ERROR";
    goto L2000;
L1500:
    //std::printf(" MNCONT UNABLE TO FIND FOUR POINTS.\n");
    *m_logStream << " MNCONT UNABLE TO FIND FOUR POINTS.\n";
//    fU[ke1-1] = u1min;
//    fU[ke2-1] = u2min;
    m_userParameterValue[ke1] = u1min;
    m_userParameterValue[ke2] = u2min;
    ierrf     = 0;
    fCstatu   = "FAILED";
L2000:
    fCfrom  = "MNContour ";
    fNfcnfr = nfcnco;
    *m_logStream << endl;
} /* mncont_ */

//______________________________________________________________________________
//void URMinuit::mncrck(const string& cardbuf, Int_urt maxcwd, const string &comand, Int_urt &lnc,
// header change (mrs43): 
//   seems as if the point of this is to modify comand - should not be a const reference
//   cardbuf also appears to be an input only -- pass by value instead of const reference
void URMinuit::mncrck(string cardbuf, Int_urt maxcwd, string& comand, Int_urt &lnc,
        Int_urt mxp, Double_urt *plist, Int_urt &llist, Int_urt &ierr, Int_urt)
{
//*-*-*-*-*-*-*-*-*-*-*-*Cracks the free-format input*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    ============================
//*-*       Cracks the free-format input, expecting zero or more
//*-*         alphanumeric fields (which it joins into COMAND(1:LNC))
//*-*         followed by one or more numeric fields separated by
//*-*         blanks and/or one comma.  The numeric fields are put into
//*-*         the LLIST (but at most MXP) elements of PLIST.
//*-*      IERR = 0 if no errors,
//*-*           = 1 if error(s).
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    /* Initialized data */

    char *cnull  = 0;
    const char *cnumer = "123456789-.0+";

    /* Local variables */
    Int_urt ifld, iend, lend, left, nreq, ipos, kcmnd, nextb, ic, ibegin, ltoadd;
    Int_urt ielmnt, lelmnt[25], nelmnt;
    string ctemp;
    char *celmnt[25];
    char command[25];

    /* Function Body */
    char *crdbuf = const_cast<char*>( cardbuf.c_str() );
    lend   = cardbuf.length();
    ielmnt = 0;
    nextb  = 1;
    ierr   = 0;
//*-*-                                  . . . .  loop over words CELMNT
L10:
    for (ipos = nextb; ipos <= lend; ++ipos) {
	ibegin = ipos;
	if (crdbuf[ipos-1] == ' ') continue;
	if (crdbuf[ipos-1] == ',') goto L250;
	goto L150;
    }
    goto L300;
L150:
//*-*-              found beginning of word, look for end
    for (ipos = ibegin + 1; ipos <= lend; ++ipos) {
	if (crdbuf[ipos-1] == ' ') goto L250;
	if (crdbuf[ipos-1] == ',') goto L250;
    }
    ipos = lend + 1;
L250:
    iend = ipos - 1;
    ++ielmnt;
	if (iend >= ibegin) celmnt[ielmnt-1] = &crdbuf[ibegin-1];
    else                celmnt[ielmnt-1] = cnull;
    lelmnt[ielmnt-1] = iend - ibegin + 1;
    if (lelmnt[ielmnt-1] > 19) {
      //std::printf(" MINUIT WARNING: INPUT DATA WORD TOO LONG.\n");
	*m_logStream << " MINUIT WARNING: INPUT DATA WORD TOO LONG." << endl;
        ctemp = cardbuf.substr(ibegin-1,iend-ibegin+1);
        //std::printf("     ORIGINAL:%s\n",ctemp.c_str());
	*m_logStream << "     ORIGINAL:" << ctemp << '\n';
	//std::printf(" TRUNCATED TO:%s\n",celmnt[ielmnt-1]);
	*m_logStream << " TRUNCATED TO:" << celmnt[ielmnt-1] << endl;
	lelmnt[ielmnt-1] = 19;
    }
    if (ipos >= lend) goto L300;
    if (ielmnt >= 25) goto L300;
//*-*-                    look for comma or beginning of next word
    for (ipos = iend + 1; ipos <= lend; ++ipos) {
	if (crdbuf[ipos-1] == ' ') continue;
	nextb = ipos;
	if (crdbuf[ipos-1] == ',') nextb = ipos + 1;
	goto L10;
    }
//*-*-                All elements found, join the alphabetic ones to
//*-*-                               form a command
L300:
    nelmnt      = ielmnt;
    command[0]  = ' '; command[1] = 0;
    lnc         = 1;
    plist[0]    = 0;
    llist       = 0;
    if (ielmnt == 0) goto L900;
    kcmnd = 0;
    for (ielmnt = 1; ielmnt <= nelmnt; ++ielmnt) {
	if ( celmnt[ielmnt-1] == cnull) goto L450;
	for (ic = 1; ic <= 13; ++ic) {
	    if (*celmnt[ielmnt-1] == cnumer[ic-1]) goto L450;
	}
	if (kcmnd >= maxcwd) continue;
	left = maxcwd - kcmnd;
	ltoadd = lelmnt[ielmnt-1];
	if (ltoadd > left) ltoadd = left;
	strncpy(&command[kcmnd],celmnt[ielmnt-1],ltoadd);
	kcmnd += ltoadd;
	if (kcmnd == maxcwd) continue;
	command[kcmnd] = ' ';
	++kcmnd;
	command[kcmnd] = 0;
    }
    lnc = kcmnd;
    goto L900;
L450:
    lnc = kcmnd;
//*-*-                     . . . .  we have come to a numeric field
    llist = 0;
    for (ifld = ielmnt; ifld <= nelmnt; ++ifld) {
	++llist;
	if (llist > mxp) {
	    nreq = nelmnt - ielmnt + 1;
	    //std::printf(" MINUIT WARNING IN MNCRCK: \n");
	    *m_logStream << " MINUIT WARNING IN MNCRCK: \n";
            //std::printf(" COMMAND HAS INPUT%5d NUMERIC FIELDS, BUT MINUIT CAN ACCEPT ONLY%3d\n",nreq,mxp);
	    *m_logStream << " COMMAND HAS INPUT" << setw(5) << nreq << " NUMERIC FIELDS, BUT MINUIT CAN ACCEPT ONLY"
		 << setw(3) << mxp << endl;
	    goto L900;
	}
	if (celmnt[ifld-1] == cnull) plist[llist-1] = 0;
	else {
            sscanf(celmnt[ifld-1],"%lf",&plist[llist-1]);
	}
    }
//*-*-                                 end loop over numeric fields
L900:
    if (lnc <= 0) lnc = 1;
	comand = command;
} /* mncrck_ */

//______________________________________________________________________________
void URMinuit::mncros(Double_urt &aopt, Int_urt &iercr)
{
//*-*-*-*-*-*-*-*-*-*-*Find point where MNEVAL=AMIN+UP*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ===============================
//*-*       Find point where MNEVAL=AMIN+UP, along the line through
//*-*       XMIDCR,YMIDCR with direction XDIRCR,YDIRCR,   where X and Y
//*-*       are parameters KE1CR and KE2CR.  If KE2CR=0 (from MINOS),
//*-*       only KE1CR is varied.  From MNCONT, both are varied.
//*-*       Crossing point is at
//*-*        (U(KE1),U(KE2)) = (XMID,YMID) + AOPT*(XDIR,YDIR)
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt alsb[3], flsb[3], bmin, bmax, zmid, sdev, zdir, zlim;
    Double_urt coeff[3], aleft, aulim, fdist, adist, aminsv;
    Double_urt anext, fnext, slope, s1, s2, x1, x2, ecarmn, ecarmx;
    Double_urt determ, rt, smalla, aright, aim, tla, tlf, dfda,ecart;
    Int_urt iout=0, i, ileft, ierev, maxlk, ibest, ik, it;
    Int_urt noless, iworst=0, iright, itoohi, kex, ipt;
    Bool_urt ldebug;
    const char *chsign;
    x2 = 0;

    ldebug = fIdbg[6] >= 1;
    aminsv = fAmin;
//*-*-       convergence when F is within TLF of AIM and next prediction
//*-*-       of AOPT is within TLA of previous value of AOPT
    aim      = fAmin + fUp;
    tlf      = fUp*.01;
    tla      = .01;
    fXpt[0]  = 0;
    fYpt[0]  = aim;
    fChpt[0] = ' ';
    ipt = 1;
    if (fKe2cr == 0) {
	fXpt[1]  = -1;
	fYpt[1]  = fAmin;
	fChpt[1] = '.';
	ipt      = 2;
    }
//*-*-                   find the largest allowed A
    aulim = 100;
    for (ik = 1; ik <= 2; ++ik) {
	if (ik == 1) {
	    kex  = fKe1cr;
	    zmid = fXmidcr;
	    zdir = fXdircr;
	} else {
	    if (fKe2cr == 0) continue;
	    kex  = fKe2cr;
	    zmid = fYmidcr;
	    zdir = fYdircr;
	}
//       if (fNvarl[kex-1] <= 1) continue;
	if (m_userParameterFlag[kex] <= 1) continue;
	if (zdir == 0) continue;
	zlim = fAlim[kex-1];
	if (zdir > 0) zlim = fBlim[kex-1];
	aulim = URMath::Min(aulim,(zlim - zmid) / zdir);
    }
//*-*-                 LSB = Line Search Buffer
//*-*-         first point
    anext   = 0;
    aopt    = anext;
    fLimset = kurFALSE;
    if (aulim < aopt + tla) fLimset = kurTRUE;
    mneval(anext, fnext, ierev);
//*-* debug printout:
    if (ldebug) {
	std::printf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f\n",fNfcn,aim,fnext,aopt);
    }
    if (ierev > 0) goto L900;
    if (fLimset && fnext <= aim) goto L930;
    ++ipt;
    fXpt[ipt-1]  = anext;
    fYpt[ipt-1]  = fnext;
    fChpt[ipt-1] = charal[ipt-1];
    alsb[0] = anext;
    flsb[0] = fnext;
    fnext   = URMath::Max(fnext,aminsv + fUp*.1);
    aopt    = URMath::Sqrt(fUp / (fnext - aminsv)) - 1;
    if (URMath::Abs(fnext - aim) < tlf) goto L800;

    if (aopt < -.5)aopt = -.5;
    if (aopt > 1)  aopt = 1;
    fLimset = kurFALSE;
    if (aopt > aulim) {
	aopt    = aulim;
	fLimset = kurTRUE;
    }
    mneval(aopt, fnext, ierev);
//*-* debug printout:
    if (ldebug) {
	std::printf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f\n",fNfcn,aim,fnext,aopt);
    }
    if (ierev > 0) goto L900;
    if (fLimset && fnext <= aim) goto L930;
    alsb[1] = aopt;
    ++ipt;
    fXpt[ipt-1]  = alsb[1];
    fYpt[ipt-1]  = fnext;
    fChpt[ipt-1] = charal[ipt-1];
    flsb[1]      = fnext;
    dfda         = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
//*-*-                  DFDA must be positive on the contour
    if (dfda > 0) goto L460;
L300:
    mnwarn("D", "MNCROS    ", "Looking for slope of the right sign");
    maxlk = 15 - ipt;
    for (it = 1; it <= maxlk; ++it) {
	alsb[0] = alsb[1];
	flsb[0] = flsb[1];
	aopt    = alsb[0] + Double_urt(it)*.2;
	fLimset = kurFALSE;
	if (aopt > aulim) {
	    aopt    = aulim;
	    fLimset = kurTRUE;
	}
	mneval(aopt, fnext, ierev);
//*-* debug printout:
	if (ldebug) {
            std::printf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f\n",fNfcn,aim,fnext,aopt);
	}
	if (ierev > 0) goto L900;
	if (fLimset && fnext <= aim) goto L930;
	alsb[1] = aopt;
	++ipt;
	fXpt[ipt-1]  = alsb[1];
	fYpt[ipt-1]  = fnext;
	fChpt[ipt-1] = charal[ipt-1];
	flsb[1]      = fnext;
	dfda         = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
	if (dfda > 0) goto L450;
    }
    mnwarn("W", "MNCROS    ", "Cannot find slope of the right sign");
    goto L950;
L450:
//*-*-                   we have two points with the right slope
L460:
    aopt  = alsb[1] + (aim - flsb[1]) / dfda;
    fdist = URMath::Min(URMath::Abs(aim - flsb[0]),URMath::Abs(aim - flsb[1]));
    adist = URMath::Min(URMath::Abs(aopt - alsb[0]),URMath::Abs(aopt - alsb[1]));
    tla = .01;
    if (URMath::Abs(aopt) > 1) tla = URMath::Abs(aopt)*.01;
    if (adist < tla && fdist < tlf) goto L800;
    if (ipt >= 15) goto L950;
    bmin = URMath::Min(alsb[0],alsb[1]) - 1;
    if (aopt < bmin) aopt = bmin;
    bmax = URMath::Max(alsb[0],alsb[1]) + 1;
    if (aopt > bmax) aopt = bmax;
//*-*-                   Try a third point
    fLimset = kurFALSE;
    if (aopt > aulim) {
	aopt    = aulim;
	fLimset = kurTRUE;
    }
    mneval(aopt, fnext, ierev);
//*-* debug printout:
    if (ldebug) {
        std::printf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f\n",fNfcn,aim,fnext,aopt);
    }
    if (ierev > 0) goto L900;
    if (fLimset && fnext <= aim) goto L930;
    alsb[2] = aopt;
    ++ipt;
    fXpt[ipt-1]  = alsb[2];
    fYpt[ipt-1]  = fnext;
    fChpt[ipt-1] = charal[ipt-1];
    flsb[2]      = fnext;
//*-*-               now we have three points, ask how many <AIM
    ecarmn = URMath::Abs(fnext-aim);
    ibest  = 3;
    ecarmx = 0;
    noless = 0;
    for (i = 1; i <= 3; ++i) {
	ecart = URMath::Abs(flsb[i-1] - aim);
	if (ecart > ecarmx) { ecarmx = ecart; iworst = i; }
	if (ecart < ecarmn) { ecarmn = ecart; ibest = i; }
	if (flsb[i-1] < aim) ++noless;
    }
//*-*-          if at least one on each side of AIM, fit a parabola
    if (noless == 1 || noless == 2) goto L500;
//*-*-          if all three are above AIM, third must be closest to AIM
    if (noless == 0 && ibest != 3) goto L950;
//*-*-          if all three below, and third is not best, then slope
//*-*-            has again gone negative, look for positive slope.
    if (noless == 3 && ibest != 3) {
	alsb[1] = alsb[2];
	flsb[1] = flsb[2];
	goto L300;
    }
//*-*-          in other cases, new straight line thru last two points
    alsb[iworst-1] = alsb[2];
    flsb[iworst-1] = flsb[2];
    dfda = (flsb[1] - flsb[0]) / (alsb[1] - alsb[0]);
    goto L460;
//*-*-               parabola fit
L500:
    mnpfit(alsb, flsb, 3, coeff, sdev);
    if (coeff[2] <= 0) {
	mnwarn("D", "MNCROS    ", "Curvature is negative near contour line.");
    }
    determ = coeff[1]*coeff[1] - coeff[2]*4*(coeff[0] - aim);
    if (determ <= 0) {
	mnwarn("D", "MNCROS    ", "Problem 2, impossible determinant");
	goto L950;
    }
//*-*-               Find which root is the right one
    rt = URMath::Sqrt(determ);
    x1 = (-coeff[1] + rt) / (coeff[2]*2);
    x2 = (-coeff[1] - rt) / (coeff[2]*2);
    s1 = coeff[1] + x1*2*coeff[2];
    s2 = coeff[1] + x2*2*coeff[2];
    if (s1*s2 > 0) {
      //std::printf(" MNCONTour problem 1\n");
	*m_logStream << " MNCONTour problem 1" << endl;
    }
    aopt  = x1;
    slope = s1;
    if (s2 > 0) {
	aopt  = x2;
	slope = s2;
    }
//*-*-        ask if converged
    tla = .01;
    if (URMath::Abs(aopt) > 1) tla = URMath::Abs(aopt)*.01;
    if (URMath::Abs(aopt - alsb[ibest-1]) < tla && URMath::Abs(flsb[ibest-1] - aim) < tlf) {
	goto L800;
    }
    if (ipt >= 15) goto L950;

//*-*-        see if proposed point is in acceptable zone between L and R
//*-*-        first find ILEFT, IRIGHT, IOUT and IBEST
    ileft  = 0;
    iright = 0;
    ibest  = 1;
    ecarmx = 0;
    ecarmn = URMath::Abs(aim - flsb[0]);
    for (i = 1; i <= 3; ++i) {
	ecart = URMath::Abs(flsb[i-1] - aim);
	if (ecart < ecarmn) { ecarmn = ecart; ibest = i; }
	if (ecart > ecarmx) { ecarmx = ecart; }
	if (flsb[i-1] > aim) {
	    if (iright == 0) iright = i;
	    else if (flsb[i-1] > flsb[iright-1]) iout = i;
	    else { iout = iright; iright = i; }
	}
        else if (ileft == 0) ileft = i;
	else if (flsb[i-1] < flsb[ileft-1]) iout = i;
	else { iout = ileft; ileft = i;	}
    }
//*-*-      avoid keeping a very bad point next time around
    if (ecarmx > URMath::Abs(flsb[iout-1] - aim)*10) {
	aopt = aopt*.5 + (alsb[iright-1] + alsb[ileft-1])*.25;
    }
//*-*-        knowing ILEFT and IRIGHT, get acceptable window
    smalla = tla*.1;
    if (slope*smalla > tlf) smalla = tlf / slope;
    aleft  = alsb[ileft-1] + smalla;
    aright = alsb[iright-1] - smalla;
//*-*-        move proposed point AOPT into window if necessary
    if (aopt < aleft)   aopt = aleft;
    if (aopt > aright)  aopt = aright;
    if (aleft > aright) aopt = (aleft + aright)*.5;

//*-*-        see if proposed point outside limits (should be impossible!)
    fLimset = kurFALSE;
    if (aopt > aulim) {
	aopt    = aulim;
	fLimset = kurTRUE;
    }
//*-*-                 Evaluate function at new point AOPT
    mneval(aopt, fnext, ierev);
//*-* debug printout:
    if (ldebug) {
        std::printf(" MNCROS: calls=%8d   AIM=%10.5f  F,A=%10.5f%10.5f\n",fNfcn,aim,fnext,aopt);
    }
    if (ierev > 0) goto L900;
    if (fLimset && fnext <= aim) goto L930;
    ++ipt;
    fXpt[ipt-1]  = aopt;
    fYpt[ipt-1]  = fnext;
    fChpt[ipt-1] = charal[ipt-1];
//*-*-               Replace odd point by new one
    alsb[iout-1] = aopt;
    flsb[iout-1] = fnext;
//*-*-         the new point may not be the best, but it is the only one
//*-*-         which could be good enough to pass convergence criteria
    ibest = iout;
    goto L500;

//*-*-      Contour has been located, return point to MNCONT OR MINOS
L800:
    iercr = 0;
    goto L1000;
//*-*-               error in the minimization
L900:
    if (ierev == 1) goto L940;
    goto L950;
//*-*-               parameter up against limit
L930:
    iercr = 1;
    goto L1000;
//*-*-               too many calls to FCN
L940:
    iercr = 2;
    goto L1000;
//*-*-               cannot find next point
L950:
    iercr = 3;
//*-*-               in any case
L1000:
    if (ldebug) {
	itoohi = 0;
	for (i = 1; i <= ipt; ++i) {
	    if (fYpt[i-1] > aim + fUp) {
		fYpt[i-1]  = aim + fUp;
		fChpt[i-1] = '+';
		itoohi     = 1;
	    }
	}
	fChpt[ipt] = 0;
        chsign = "POSI";
	if (fXdircr < 0) chsign = "NEGA";
	if (fKe2cr == 0) {
	  //std::printf("  %sTIVE MINOS ERROR, PARAMETER %3d\n",chsign,fKe1cr);
	    *m_logStream << chsign << "TIVE MINOS ERROR, PARAMETER " << setw(3) << fKe1cr << endl;
	}
	if (itoohi == 1) {
	  //std::printf("POINTS LABELLED '+' WERE TOO HIGH TO PLOT.\n");
	    *m_logStream << "POINTS LABELLED '+' WERE TOO HIGH TO PLOT." << endl;
	}
	if (iercr == 1) {
	  //std::printf("RIGHTMOST POINT IS UP AGAINST LIMIT.\n");
	    *m_logStream << "RIGHTMOST POINT IS UP AGAINST LIMIT." << endl;
	}
	mnplot(fXpt, fYpt, fChpt, ipt, fNpagwd, fNpagln);
    }
} /* mncros_ */

//______________________________________________________________________________
void URMinuit::mncuve()
{
//*-*-*-*-*-*-*-*Makes sure that the current point is a local minimum*-*-*-*-*
//*-*            ====================================================
//*-*        Makes sure that the current point is a local
//*-*        minimum and that the error matrix exists,
//*-*        or at least something good enough for MINOS and MNCONT
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dxdi, wint;
    Int_urt ndex, iext, i, j;

    if (fISW[3] < 1) {
      //std::printf(" FUNCTION MUST BE MINIMIZED BEFORE CALLING %s\n",fCfrom.c_str());
	*m_logStream << " FUNCTION MUST BE MINIMIZED BEFORE CALLING " << fCfrom << endl;
	fApsi = fEpsi;
	mnmigr();
    }
    if (fISW[1] < 3) {
	mnhess();
	if (fISW[1] < 1) {
	    mnwarn("W", fCfrom.c_str(), "NO ERROR MATRIX.  WILL IMPROVISE.\n");
	    for (i = 1; i <= fNpar; ++i) {
		ndex = i*(i-1) / 2;
		for (j = 1; j <= i-1; ++j) {
		    ++ndex;
		    fVhmat[ndex-1] = 0;
		}
		++ndex;
		if (fG2[i-1] <= 0) {
		    wint = fWerr[i-1];
		    iext = fNexofi[i-1];
//                    if (fNvarl[iext-1] > 1) {
		    if (m_userParameterFlag[iext] > 1) {
			mndxdi(fX[i-1], i-1, dxdi);
			if (URMath::Abs(dxdi) < .001) wint = .01;
			else                   wint /= URMath::Abs(dxdi);
		    }
		    fG2[i-1] = fUp / (wint*wint);
		}
		fVhmat[ndex-1] = 2 / fG2[i-1];
	    }
	    fISW[1] = 1;
	    fDcovar = 1;
	} else  mnwerr();
    }
} /* mncuve_ */

//______________________________________________________________________________
void URMinuit::mnderi()
{
//*-*-*-*-*-*-*-*Calculates the first derivatives of FCN (GRD)*-*-*-*-*-*-*-*
//*-*            =============================================
//*-*        Calculates the first derivatives of FCN (GRD),
//*-*        either by finite differences or by transforming the user-
//*-*        supplied derivatives to internal coordinates,
//*-*        according to whether ISW(3) is zero or one.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt step, dfmin, stepb4, dd, df, fs1;
    Double_urt tlrstp, tlrgrd, epspri, optstp, stpmax, stpmin, fs2, grbfor=0, d1d2, xtf;
    Int_urt icyc, ncyc, iint, iext, i, nparx;
    Bool_urt ldebug;

    nparx = fNpar;
    ldebug = fIdbg[2] >= 1;
    if (fAmin == fUndefi) mnamin();
    if (fISW[2] == 1) goto L100;

    if (ldebug) {
//*-*-                      make sure starting at the right place
	mninex(fX);
	nparx = fNpar;
//	Eval(nparx, fGin, fs1, fU, 4);	++fNfcn;
        Eval(nparx, fGin, fs1, m_userParameterValue, 4);	++fNfcn;
	if (fs1 != fAmin) {
	    df    = fAmin - fs1;
           ostringstream warning;
           warning << "function value differs from AMIN by " << df;
	    mnwarn("D", "MNDERI", warning.str().c_str());
	    fAmin = fs1;
	}
	//std::printf("  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI\n");
        //std::printf(" PAR    DERIV     STEP      MINSTEP   OPTSTEP  D1-D2    2ND DRV\n");
	*m_logStream << "  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI\n";
	*m_logStream << " PAR    DERIV     STEP      MINSTEP   OPTSTEP  D1-D2    2ND DRV" << endl;
    }
    dfmin = fEpsma2*8*(URMath::Abs(fAmin) + fUp);
    if (fIstrat <= 0) {
	ncyc   = 2;
	tlrstp = .5;
	tlrgrd = .1;
    } else if (fIstrat == 1) {
	ncyc   = 3;
	tlrstp = .3;
	tlrgrd = .05;
    } else {
	ncyc   = 5;
	tlrstp = .1;
	tlrgrd = .02;
    }
//*-*-                               loop over variable parameters
    for (i = 1; i <= fNpar; ++i) {
	epspri = fEpsma2 + URMath::Abs(fGrd[i-1]*fEpsma2);
//*-*-        two-point derivatives always assumed necessary
//*-*-        maximum number of cycles over step size depends on strategy
	xtf = fX[i-1];
	stepb4 = 0;
//*-*-                              loop as little as possible here!/
	for (icyc = 1; icyc <= ncyc; ++icyc) {
//*-*-                ........ theoretically best step
	    optstp = URMath::Sqrt(dfmin / (URMath::Abs(fG2[i-1]) + epspri));
//*-*-                    step cannot decrease by more than a factor of ten
	    step = URMath::Max(optstp,URMath::Abs(fGstep[i-1]*.1));
//*-*-                but if parameter has limits, max step size = 0.5
	    if (fGstep[i-1] < 0 && step > .5) step = .5;
//*-*-                and not more than ten times the previous step
	    stpmax = URMath::Abs(fGstep[i-1])*10;
	    if (step > stpmax) step = stpmax;
//*-*-                minimum step size allowed by machine precision
	    stpmin = URMath::Abs(fEpsma2*fX[i-1])*8;
	    if (step < stpmin) step = stpmin;
//*-*-                end of iterations if step change less than factor 2
	    if (URMath::Abs((step - stepb4) / step) < tlrstp) goto L50;
//*-*-        take step positive
	    stepb4 = step;
	    if (fGstep[i-1] > 0) fGstep[i-1] =  URMath::Abs(step);
	    else                 fGstep[i-1] = -URMath::Abs(step);
	    stepb4  = step;
	    fX[i-1] = xtf + step;
	    mninex(fX);
//            Eval(nparx, fGin, fs1, fU, 4);	    ++fNfcn;
	    Eval(nparx, fGin, fs1, m_userParameterValue, 4);	    ++fNfcn;
//*-*-        take step negative
	    fX[i-1] = xtf - step;
	    mninex(fX);
//            Eval(nparx, fGin, fs2, fU, 4);	    ++fNfcn;
	    Eval(nparx, fGin, fs2, m_userParameterValue, 4);	    ++fNfcn;
	    grbfor = fGrd[i-1];
	    fGrd[i-1] = (fs1 - fs2) / (step*2);
	    fG2[i-1]  = (fs1 + fs2 - fAmin*2) / (step*step);
	    fX[i-1]   = xtf;
	    if (ldebug) {
		d1d2 = (fs1 + fs2 - fAmin*2) / step;
		std::printf("%4d%11.3g%11.3g%10.2g%10.2g%10.2g%10.2g\n",i,fGrd[i-1],step,stpmin,optstp,d1d2,fG2[i-1]);
	    }
//*-*-        see if another iteration is necessary
	    if (URMath::Abs(grbfor - fGrd[i-1]) / (URMath::Abs(fGrd[i-1]) + dfmin/step) < tlrgrd)
		goto L50;
	}
//*-*-                          end of ICYC loop. too many iterations
	if (ncyc == 1) goto L50;
        { ostringstream warning2;
        warning2 << "First derivative not converged. " << fGrd[i-1] << grbfor;
        mnwarn("D", "MNDERI", warning2.str().c_str()); }
L50:
	;
    }
    mninex(fX);
    return;
//*-*-                                       .  derivatives calc by fcn
L100:
    for (iint = 1; iint <= fNpar; ++iint) {
	iext = fNexofi[iint-1];
//       if (fNvarl[iext-1] <= 1) {
	if (m_userParameterFlag[iext] <= 1) {
           fGrd[iint-1] = fGin[iext-1];
	} else {
           dd = (fBlim[iext-1] - fAlim[iext-1])*.5*URMath::Cos(fX[iint-1]);
           fGrd[iint-1] = fGin[iext-1]*dd;
        }
    }
} /* mnderi_ */

//______________________________________________________________________________
void URMinuit::mndxdi(Double_urt pint, Int_urt ipar, Double_urt &dxdi)
{
//*-*-*-*Calculates the transformation factor between ext/internal values*-*
//*-*    =====================================================================
//*-*        calculates the transformation factor between external and
//*-*        internal parameter values.     this factor is one for
//*-*        parameters which are not limited.     called from MNEMAT.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt i = fNexofi[ipar];
    dxdi = 1;
//    if (fNvarl[i-1] > 1) {
    if (m_userParameterFlag[i] > 1) {
	dxdi = URMath::Abs((fBlim[i-1] - fAlim[i-1])*URMath::Cos(pint))*.5;
    }
} /* mndxdi_ */

//______________________________________________________________________________
void URMinuit::mneig(Double_urt *a, Int_urt ndima, Int_urt n, Int_urt mits, Double_urt *work, Double_urt precis, Int_urt &ifault)
{
//*-*-*-*-*-*-*-*-*-*-*-*Compute matrix eigen values*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    ===========================
    /* System generated locals */
    Int_urt a_offset;
    Double_urt d__1;

    /* Local variables */
    Double_urt b, c, f, h, r, s, hh, gl, pr, pt;
    Int_urt i, j, k, l, m=0, i0, i1, j1, m1, n1;

//*-*-         PRECIS is the machine precision EPSMAC
    /* Parameter adjustments */
    a_offset = ndima + 1;
    a -= a_offset;
    --work;

    /* Function Body */
    ifault = 1;

    i = n;
    for (i1 = 2; i1 <= n; ++i1) {
	l  = i-2;
	f  = a[i + (i-1)*ndima];
	gl = 0;

	if (l < 1) goto L25;

	for (k = 1; k <= l; ++k) {
	    d__1 = a[i + k*ndima];
	    gl  += d__1*d__1;
	}
L25:
	h = gl + f*f;

	if (gl > 1e-35) goto L30;

	work[i]     = 0;
	work[n + i] = f;
	goto L65;
L30:
	++l;
	gl = URMath::Sqrt(h);
	if (f >= 0) gl = -gl;
	work[n + i] = gl;
	h -= f*gl;
	a[i + (i-1)*ndima] = f - gl;
	f = 0;
	for (j = 1; j <= l; ++j) {
	    a[j + i*ndima] = a[i + j*ndima] / h;
	    gl = 0;
	    for (k = 1; k <= j; ++k) { gl += a[j + k*ndima]*a[i + k*ndima]; }
	    if (j >= l) goto L47;
	    j1 = j + 1;
	    for (k = j1; k <= l; ++k) {	gl += a[k + j*ndima]*a[i + k*ndima]; }
L47:
	    work[n + j] = gl / h;
	    f += gl*a[j + i*ndima];
	}
	hh = f / (h + h);
	for (j = 1; j <= l; ++j) {
	    f  = a[i + j*ndima];
	    gl = work[n + j] - hh*f;
	    work[n + j] = gl;
	    for (k = 1; k <= j; ++k) {
		a[j + k*ndima] = a[j + k*ndima] - f*work[n + k] - gl*a[i + k*ndima];
	    }
	}
	work[i] = h;
L65:
	--i;
    }
    work[1] = 0;
    work[n + 1] = 0;
    for (i = 1; i <= n; ++i) {
	l = i-1;
	if (work[i] == 0 || l == 0) goto L100;

	for (j = 1; j <= l; ++j) {
	    gl = 0;
	    for (k = 1; k <= l; ++k) { gl += a[i + k*ndima]*a[k + j*ndima]; }
	    for (k = 1; k <= l; ++k) { a[k + j*ndima] -= gl*a[k + i*ndima]; }
	}
L100:
	work[i] = a[i + i*ndima];
	a[i + i*ndima] = 1;
	if (l == 0) continue;

	for (j = 1; j <= l; ++j) {
	    a[i + j*ndima] = 0;
	    a[j + i*ndima] = 0;
	}
    }

    n1 = n - 1;
    for (i = 2; i <= n; ++i) {
	i0 = n + i-1;
	work[i0] = work[i0 + 1];
    }
    work[n + n] = 0;
    b = 0;
    f = 0;
    for (l = 1; l <= n; ++l) {
	j = 0;
	h = precis*(URMath::Abs(work[l]) + URMath::Abs(work[n + l]));
	if (b < h) b = h;
	for (m1 = l; m1 <= n; ++m1) {
	    m = m1;
	    if (URMath::Abs(work[n + m]) <= b)	goto L150;
	}

L150:
	if (m == l) goto L205;

L160:
	if (j == mits) return;
	++j;
	pt = (work[l + 1] - work[l]) / (work[n + l]*2);
	r  = URMath::Sqrt(pt*pt + 1);
	pr = pt + r;
	if (pt < 0) pr = pt - r;

	h = work[l] - work[n + l] / pr;
	for (i = l; i <= n; ++i) { work[i] -= h; }
	f += h;
	pt = work[m];
	c  = 1;
	s  = 0;
	m1 = m - 1;
	i  = m;
	for (i1 = l; i1 <= m1; ++i1) {
	    j = i;
	    --i;
	    gl = c*work[n + i];
	    h  = c*pt;
	    if (URMath::Abs(pt) >= URMath::Abs(work[n + i])) goto L180;

	    c = pt / work[n + i];
	    r = URMath::Sqrt(c*c + 1);
	    work[n + j] = s*work[n + i]*r;
	    s  = 1 / r;
	    c /= r;
	    goto L190;
L180:
	    c = work[n + i] / pt;
	    r = URMath::Sqrt(c*c + 1);
	    work[n + j] = s*pt*r;
	    s = c / r;
	    c = 1 / r;
L190:
	    pt = c*work[i] - s*gl;
	    work[j] = h + s*(c*gl + s*work[i]);
	    for (k = 1; k <= n; ++k) {
		h = a[k + j*ndima];
		a[k + j*ndima] = s*a[k + i*ndima] + c*h;
		a[k + i*ndima] = c*a[k + i*ndima] - s*h;
	    }
	}
	work[n + l] = s*pt;
	work[l]     = c*pt;

	if (URMath::Abs(work[n + l]) > b) goto L160;

L205:
	work[l] += f;
    }
    for (i = 1; i <= n1; ++i) {
	k  = i;
	pt = work[i];
	i1 = i + 1;
	for (j = i1; j <= n; ++j) {
	    if (work[j] >= pt) continue;
	    k  = j;
	    pt = work[j];
	}

	if (k == i) continue;

	work[k] = work[i];
	work[i] = pt;
	for (j = 1; j <= n; ++j) {
	    pt = a[j + i*ndima];
	    a[j + i*ndima] = a[j + k*ndima];
	    a[j + k*ndima] = pt;
	}
    }
    ifault = 0;
} /* mneig_ */

//______________________________________________________________________________
void URMinuit::mnemat(Double_urt *emat, Int_urt ndim)
{
// Calculates the external error matrix from the internal matrix
//
// Note that if the matrix is declared like Double_urt matrix[5][5]
// in the calling program, one has to call mnemat with, eg
//     gMinuit->mnemat(&matrix[0][0],5);

    /* System generated locals */
    Int_urt emat_dim1, emat_offset;

    /* Local variables */
    Double_urt dxdi, dxdj;
    Int_urt i, j, k, npard, k2, kk, iz, nperln, kga, kgb;
    string ctemp;

    /* Parameter adjustments */
    emat_dim1 = ndim;
    emat_offset = emat_dim1 + 1;
    emat -= emat_offset;

    /* Function Body */
    if (fISW[1] < 1) return;
    if (fISW[4] >= 2) {
      //std::printf(" EXTERNAL ERROR MATRIX.    NDIM=%4d    NPAR=%3d    ERR DEF=%g\n",ndim,fNpar,fUp);
	*m_logStream << " EXTERNAL ERROR MATRIX.    NDIM=" << setw(4) << ndim << "    NPAR=" << setw(3)
	     << fNpar << "    ERR DEF=" << fUp << endl;
    }
//*-*-                   size of matrix to be printed
    npard = fNpar;
    if (ndim < fNpar) {
	npard = ndim;
	if (fISW[4] >= 0) {
	    std::printf(" USER-DIMENSIONED  ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.\n");
	}
    }
//*-*-                NPERLN is the number of elements that fit on one line

    nperln = (fNpagwd - 5) / 10;
    nperln = URMath::Min(nperln,13);
    if (fISW[4] >= 1 && npard > nperln) {
      //  disable annoying printout when obtaining error matrix for fits with
      //  large numbers of parameters
      //	std::printf(" ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.\n");
    }
//*-*-                I counts the rows of the matrix
    for (i = 1; i <= npard; ++i) {
	mndxdi(fX[i-1], i-1, dxdi);
	kga = i*(i-1) / 2;
	for (j = 1; j <= i; ++j) {
	    mndxdi(fX[j-1], j-1, dxdj);
	    kgb = kga + j;
	    emat[i + j*emat_dim1] = dxdi*fVhmat[kgb-1]*dxdj*fUp;
	    emat[j + i*emat_dim1] = emat[i + j*emat_dim1];
	}
    }
//*-*-                   IZ is number of columns to be printed in row I
    if (fISW[4] >= 2) {
	for (i = 1; i <= npard; ++i) {
	    iz = npard;
	    if (npard >= nperln) iz = i;
	    for (k = 1; nperln < 0 ? k >= iz : k <= iz; k += nperln) {
		k2 = k + nperln - 1;
		if (k2 > iz) k2 = iz;
                ctemp = " ";
		for (kk = k; kk <= k2; ++kk) {
                   ostringstream localTemp;
                   localTemp << " " << setw(10) << setprecision(3) << emat[i + kk*emat_dim1];
		    ctemp += localTemp.str();
	        }
                std::printf("%s\n",ctemp.c_str());
	    }
	}
    }
} /* mnemat_ */

//______________________________________________________________________________
void URMinuit::mnerrs(Int_urt number, Double_urt &eplus, Double_urt &eminus, Double_urt &eparab, Double_urt &gcc)
{
//*-*-*-*-*-*-*-*-*-*Utility routine to get MINOS errors*-*-*-*-*-*-*-*-*-*-*
//*-*                ===================================
//*-*    Called by user.
//*-*    NUMBER is the parameter number
//*-*    values returned by MNERRS:
//*-*       EPLUS, EMINUS are MINOS errors of parameter NUMBER,
//*-*       EPARAB is 'parabolic' error (from error matrix).
//*-*                 (Errors not calculated are set = 0)
//*-*       GCC is global correlation coefficient from error matrix
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Double_urt dxdi;
    Int_urt ndiag, iin, iex;

//    iex = number+1;
    iex = number;

    if (iex > fNu || iex <= 0) goto L900;
//    iin = fNiofex[iex-1];
//    cout << "mnerr: request to get error for param = " << iex << endl;
    iin = m_userParameterIdToInternalId[iex];
//    cout << "   found internal parameter iin " << iin << endl;
    if (iin <= 0) goto L900;

//*-*-            IEX is external number, IIN is internal number
    eplus  = fErp[iin-1];
    if (eplus == fUndefi)  eplus = 0;
    eminus = fErn[iin-1];
    if (eminus == fUndefi) eminus = 0;
    mndxdi(fX[iin-1], iin-1, dxdi);
    ndiag  = iin*(iin + 1) / 2;
    eparab = URMath::Abs(dxdi*URMath::Sqrt(URMath::Abs(fUp*fVhmat[ndiag- 1])));
//*-*-             global correlation coefficient
    gcc = 0;
    if (fISW[1] < 2) return;
    gcc = fGlobcc[iin-1];
    return;
//*-*-                 ERROR.  parameter number not valid
L900:
    eplus  = 0;
    eminus = 0;
    eparab = 0;
    gcc    = 0;
} /* mnerrs_ */

//______________________________________________________________________________
void URMinuit::mneval(Double_urt anext, Double_urt &fnext, Int_urt &ierev)
{
//*-*-*-*-*-*-*Evaluates the function being analyzed by MNCROS*-*-*-*-*-*-*-*
//*-*          ===============================================
//*-*      Evaluates the function being analyzed by MNCROS, which is
//*-*      generally the minimum of FCN with respect to all remaining
//*-*      variable parameters.  The class data members contains the
//*-*      data necessary to know the values of U(KE1CR) and U(KE2CR)
//*-*      to be used, namely     U(KE1CR) = XMIDCR + ANEXT*XDIRCR
//*-*      and (if KE2CR .NE. 0)  U(KE2CR) = YMIDCR + ANEXT*YDIRCR
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt nparx;

//   fU[fKe1cr-1] = fXmidcr + anext*fXdircr;
    m_userParameterValue[fKe1cr] = fXmidcr + anext*fXdircr;
//    if (fKe2cr != 0) fU[fKe2cr-1] = fYmidcr + anext*fYdircr;
    if (fKe2cr != 0) m_userParameterValue[fKe2cr] = fYmidcr + anext*fYdircr;
    mninex(fX);
    nparx = fNpar;
//    Eval(nparx, fGin, fnext, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, fnext, m_userParameterValue, 4);    ++fNfcn;
    ierev = 0;
    if (fNpar > 0) {
	fItaur = 1;
	fAmin = fnext;
	fISW[0] = 0;
	mnmigr();
	fItaur = 0;
	fnext = fAmin;
	if (fISW[0] >= 1) ierev = 1;
	if (fISW[3] < 1)  ierev = 2;
    }
} /* mneval_ */

//______________________________________________________________________________
void URMinuit::mnexcm(const string& command, Double_urt *plist, Int_urt llist, Int_urt &ierflg)
{
//*-*-*-*-*-*Interprets a command and takes appropriate action*-*-*-*-*-*-*-*
//*-*        =================================================
//*-*        either directly by skipping to the corresponding code in
//*-*        MNEXCM, or by setting up a call to a function
//*-*
//*-*  recognized MINUIT commands:
//*-*  obsolete commands:
//*-*      IERFLG is now (94.5) defined the same as ICONDN in MNCOMD
//*-*            = 0: command executed normally
//*-*              1: command is blank, ignored
//*-*              2: command line unreadable, ignored
//*-*              3: unknown command, ignored
//*-*              4: abnormal termination (e.g., MIGRAD not converged)
//*-*              9: reserved
//*-*             10: END command
//*-*             11: EXIT or STOP command
//*-*             12: RETURN command
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    string comand = command;
    static string clower = "abcdefghijklmnopqrstuvwxyz";
    static string cupper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const string cname[40] = {
      "MINImize  ",
      "SEEk      ",
      "SIMplex   ",
      "MIGrad    ",
      "MINOs     ",
      "SET xxx   ",
      "SHOw xxx  ",
      "TOP of pag",
      "FIX       ",
      "REStore   ",
      "RELease   ",
      "SCAn      ",
      "CONtour   ",
      "HESse     ",
      "SAVe      ",
      "IMProve   ",
      "CALl fcn  ",
      "STAndard  ",
      "END       ",
      "EXIt      ",
      "RETurn    ",
      "CLEar     ",
      "HELP      ",
      "MNContour ",
      "STOp      ",
      "JUMp      ",
      "          ",
      "          ",
      "          ",
      "          ",
      "          ",
      "          ",
      "          ",
      "COVARIANCE",
      "PRINTOUT  ",
      "GRADIENT  ",
      "MATOUT    ",
      "ERROR DEF ",
      "LIMITS    ",
      "PUNCH     "};

    Int_urt nntot = 40;

    /* Local variables */
    Double_urt step, xptu[101], yptu[101], f, rno;
    Int_urt icol, kcol, ierr, iint, iext, lnow, nptu, i, iflag, ierrf;
    Int_urt ilist, nparx, izero, nf, lk, it, iw, inonde, nsuper;
    Int_urt it2, ke1, ke2, nowprt, kll, krl;
    string chwhy, c26, cvblnk, cneway, comd;
    string ctemp;
    Bool_urt lfreed, ltofix, lfixed;

//*-*  alphabetical order of command names!

    /* Function Body */

    lk = comand.length();
    if (lk > 20) lk = 20;
    fCword =  comand;
    int (*pf)(int)=toupper; 
    transform(fCword.begin(), fCword.end(), fCword.begin(), pf); 
//*-*-          Copy the first MAXP arguments into WORD7, making
//*-*-          sure that WORD7(1)=0 if LLIST=0
    for (iw = 1; iw <= fMaxpar; ++iw) {
	fWord7[iw-1] = 0;
	if (iw <= llist) fWord7[iw-1] = plist[iw-1];
    }
    ++fIcomnd;
    fNfcnlc = fNfcn;
    if (fCword.substr(0,7) != "SET PRI" || fWord7[0] >= 0) {
	if (fISW[4] >= 0) {
	    lnow = llist;
	    if (lnow > 4) lnow = 4;
	    std::printf(" **********\n");
            ostringstream octemp;
            octemp << " **" << setw(5) << fIcomnd << " **" << fCword;
            ctemp = octemp.str();
//            ctemp = Form(" **%5d **%s",.c_str());
	    for (i = 1; i <= lnow; ++i) {
               ostringstream localTemp;
               localTemp << " " << plist[i-1];
                ctemp += localTemp.str();
	    }
            std::printf("%s\n",ctemp.c_str());
	    inonde = 0;
	    if (llist > lnow) {
		kll = llist;
		if (llist > fMaxpar) {
		    inonde = 1;
		    kll = fMaxpar;
		}
                std::printf(" ***********\n");
		for (i = lnow + 1; i <= kll; ++i) {
                   std::printf("%12.4g\n",plist[i-1]);
		}
	    }
	    std::printf(" **********\n");
	    if (inonde > 0) {
		std::printf("  ERROR: ABOVE CALL TO MNEXCM TRIED TO PASS MORE THAN %d PARAMETERS.\n", fMaxpar);
	    }
	}
    }
    fNfcnmx = Int_urt(fWord7[0]);
    if (fNfcnmx <= 0) {
	fNfcnmx = fNpar*100 + 200 + fNpar*fNpar*5;
    }
    fEpsi = fWord7[1];
    if (fEpsi <= 0) {
	fEpsi = fUp*.1;
    }
    fLnewmn = kurFALSE;
    fLphead = kurTRUE;
    fISW[0] = 0;
    ierflg = 0;
//*-*-               look for command in list CNAME . . . . . . . . . .
    ctemp = fCword.substr(0,3);
    for (i = 1; i <= nntot; ++i) {
	if (strncmp(ctemp.c_str(),cname[i-1].c_str(),3) == 0) goto L90;
    }
    std::printf("UNKNOWN COMMAND IGNORED:%s\n", comand.c_str());
    ierflg = 3;
    return;
//*-*-               normal case: recognized MINUIT command . . . . . . .
L90:
    if (fCword.substr(0,4) == "MINO") i = 5;
    if (i != 6 && i != 7 && i != 8 && i != 23) {
	fCfrom  = cname[i-1];
	fNfcnfr = fNfcn;
    }
//*-*-             1    2    3    4    5    6    7    8    9   10
    switch (i) {
	case 1:  goto L400;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L700;
	case 7:  goto L700;
	case 8:  goto L800;
	case 9:  goto L900;
	case 10:  goto L1000;
	case 11:  goto L1100;
	case 12:  goto L1200;
	case 13:  goto L1300;
	case 14:  goto L1400;
	case 15:  goto L1500;
	case 16:  goto L1600;
	case 17:  goto L1700;
	case 18:  goto L1800;
	case 19:  goto L1900;
	case 20:  goto L1900;
	case 21:  goto L1900;
	case 22:  goto L2200;
	case 23:  goto L2300;
	case 24:  goto L2400;
	case 25:  goto L1900;
	case 26:  goto L2600;
	case 27:  goto L3300;
	case 28:  goto L3300;
	case 29:  goto L3300;
	case 30:  goto L3300;
	case 31:  goto L3300;
	case 32:  goto L3300;
	case 33:  goto L3300;
	case 34:  goto L3400;
	case 35:  goto L3500;
	case 36:  goto L3600;
	case 37:  goto L3700;
	case 38:  goto L3800;
	case 39:  goto L3900;
	case 40:  goto L4000;
    }
//*-*-                                       . . . . . . . . . . seek
L200:
    mnseek();
    return;
//*-*-                                       . . . . . . . . . . simplex
L300:
    mnsimp();
    if (fISW[3] < 1) ierflg = 4;
    return;
//*-*-                                       . . . . . . migrad, minimize
L400:
    nf = fNfcn;
    fApsi = fEpsi;
    mnmigr();
    mnwerr();
    if (fISW[3] >= 1) return;
    ierflg = 4;
    if (fISW[0] == 1) return;
    if (fCword.substr(0,3) == "MIG") return;

    fNfcnmx = fNfcnmx + nf - fNfcn;
    nf = fNfcn;
    mnsimp();
    if (fISW[0] == 1) return;
    fNfcnmx = fNfcnmx + nf - fNfcn;
    mnmigr();
    if (fISW[3] >= 1) ierflg = 0;
    mnwerr();
    return;
//*-*-                                       . . . . . . . . . . minos
L500:
    nsuper = fNfcn + ( ( fNpar + 1 ) << 1)*fNfcnmx;
//*-*-         possible loop over new minima
    fEpsi = fUp*.1;
L510:
    mncuve();
    mnmnos();
    if (! fLnewmn) return;
    mnrset(0);
    mnmigr();
    mnwerr();
    if (fNfcn < nsuper) goto L510;
    std::printf(" TOO MANY FUNCTION CALLS. MINOS GIVES UP\n");
    ierflg = 4;
    return;
//*-*-                                       . . . . . . . . . .set, show
L700:
    mnset();
    return;
//*-*-                                       . . . . . . . . . . top of page

L800:
    std::printf("1\n");
    return;
//*-*-                                       . . . . . . . . . . fix
L900:
    ltofix = kurTRUE;
//*-*-                                       . . (also release) ....
L901:
    lfreed = kurFALSE;
    lfixed = kurFALSE;
    if (llist == 0) {
	std::printf("%s:  NO PARAMETERS REQUESTED \n",fCword.c_str());
	return;
    }
    for (ilist = 1; ilist <= llist; ++ilist) {
	iext = Int_urt(plist[ilist-1]);
	chwhy = " IS UNDEFINED.";
        if (iext <= 0) goto L930;
        if (iext > fNu) goto L930;
//	if (fNvarl[iext-1] < 0) goto L930;
        if (m_userParameterFlag[iext] < 0) goto L930;
	chwhy = " IS CONSTANT.  ";
//        if (fNvarl[iext-1] == 0) goto L930;
	if (m_userParameterFlag[iext] == 0) goto L930;
//        iint = fNiofex[iext-1];
	iint = m_userParameterIdToInternalId[iext];
	if (ltofix) {
	    chwhy = " ALREADY FIXED.";
	    if (iint == 0) goto L930;
	    mnfixp(iint-1, ierr);
	    if (ierr == 0) lfixed = kurTRUE;
	    else           ierflg = 4;
	} else {
	    chwhy = " ALREADY VARIABLE.";
	    if (iint > 0) goto L930;
	    krl = -abs(iext);
	    mnfree(krl);
	    lfreed = kurTRUE;
	}
	continue;
L930:
	std::printf(" PARAMETER%4d %s IGNORED.\n",iext,chwhy.c_str());
    }
    if (lfreed || lfixed) mnrset(0);
    if (lfreed) {
	fISW[1] = 0;
	fDcovar = 1;
	fEDM = fBigedm;
	fISW[3] = 0;
    }
    mnwerr();
    if (fISW[4] > 1) mnprin(5, fAmin);
    return;
//*-*-                                       . . . . . . . . . . restore
L1000:
    it = Int_urt(fWord7[0]);
    if (it > 1 || it < 0) goto L1005;
    lfreed = fNpfix > 0;
    mnfree(it);
    if (lfreed) {
	mnrset(0);
	fISW[1] = 0;
	fDcovar = 1;
	fEDM    = fBigedm;
    }
    return;
L1005:
    std::printf(" IGNORED.  UNKNOWN ARGUMENT:%4d\n",it);
    ierflg = 3;
    return;
//*-*-                                       . . . . . . . . . . release
L1100:
    ltofix = kurFALSE;
    goto L901;
//*-*-                                      . . . . . . . . . . scan . . .
L1200:
    iext = Int_urt(fWord7[0]);
    if (iext <= 0) goto L1210;
    it2 = 0;
//    if (iext <= fNu) it2 = fNiofex[iext-1];
    if (iext <= fNu) it2 = m_userParameterIdToInternalId[iext];
    if (it2 <= 0) goto L1250;

L1210:
    mnscan();
    return;
L1250:
    std::printf(" PARAMETER%4d NOT VARIABLE.\n",iext);
    ierflg = 3;
    return;
//*-*-                                       . . . . . . . . . . contour
L1300:
    ke1 = Int_urt(fWord7[0]);
    ke2 = Int_urt(fWord7[1]);
    if (ke1 == 0) {
	if (fNpar == 2) {
	    ke1 = fNexofi[0];
	    ke2 = fNexofi[1];
	} else {
	    std::printf("%s:  NO PARAMETERS REQUESTED \n",fCword.c_str());
	    ierflg = 3;
	    return;
	}
    }
    fNfcnmx = 1000;
//    mncntr(ke1-1, ke2-1, ierrf);
    mncntr(ke1, ke2, ierrf);
    if (ierrf > 0) ierflg = 3;
    return;
//*-*-                                       . . . . . . . . . . hesse
L1400:
    mnhess();
    mnwerr();
    if (fISW[4] >= 0) mnprin(2, fAmin);
    if (fISW[4] >= 1) mnmatu(1);
    return;
//*-*-                                       . . . . . . . . . . save
L1500:
    mnsave();
    return;
//*-*-                                       . . . . . . . . . . improve
L1600:
    mncuve();
    mnimpr();
    if (fLnewmn) goto L400;
    ierflg = 4;
    return;
//*-*-                                       . . . . . . . . . . call fcn
L1700:
    iflag = Int_urt(fWord7[0]);
    nparx = fNpar;
    f = fUndefi;
//    Eval(nparx, fGin, f, fU, iflag);    ++fNfcn;
    Eval(nparx, fGin, f, m_userParameterValue, iflag);    ++fNfcn;
    nowprt = 0;
    if (f != fUndefi) {
	if (fAmin == fUndefi) {
	    fAmin  = f;
	    nowprt = 1;
	} else if (f < fAmin) {
	    fAmin  = f;
	    nowprt = 1;
	}
	if (fISW[4] >= 0 && iflag <= 5 && nowprt == 1) {
	    mnprin(5, fAmin);
	}
	if (iflag == 3)  fFval3 = f;
    }
    if (iflag > 5) mnrset(1);
    return;
//*-*-                                       . . . . . . . . . . standard
L1800:
//    stand();
    return;
//*-*-                                      . . . return, stop, end, exit
L1900:
    it = Int_urt(fWord7[0]);
    if (fFval3 != fAmin && it == 0) {
	iflag = 3;
	std::printf(" CALL TO USER FUNCTION WITH IFLAG = 3\n");
	nparx = fNpar;
//        Eval(nparx, fGin, f, fU, iflag);	++fNfcn;
	Eval(nparx, fGin, f, m_userParameterValue, iflag);	++fNfcn;
    }
    ierflg = 11;
    if (fCword.substr(0,3) == "END") ierflg = 10;
    if (fCword.substr(0,3) == "RET") ierflg = 12;
    return;
//*-*-                                       . . . . . . . . . . clear
L2200:
    mncler();
    if (fISW[4] >= 1) {
	std::printf(" MINUIT MEMORY CLEARED. NO PARAMETERS NOW DEFINED.\n");
    }
    return;
//*-*-                                       . . . . . . . . . . help
L2300:
    kcol = 0;
    for (icol = 5; icol <= lk; ++icol) {
	if (fCword[icol-1] == ' ') continue;
	kcol = icol;
	goto L2320;
    }
L2320:
    if (kcol == 0) comd = "*   ";
    else           comd = fCword.substr(kcol-1,lk-kcol+1);
    mnhelp(comd);
    return;
//*-*-                                      . . . . . . . . . . MNContour
L2400:
    fEpsi = fUp*.05;
    ke1 = Int_urt(fWord7[0]);
    ke2 = Int_urt(fWord7[1]);
    if (ke1 == 0 && fNpar == 2) {
	ke1 = fNexofi[0];
	ke2 = fNexofi[1];
    }
    nptu = Int_urt(fWord7[2]);
    if (nptu <= 0)  nptu = 20;
    if (nptu > 101) nptu = 101;
    fNfcnmx = (nptu + 5)*100*(fNpar + 1);
//    mncont(ke1-1, ke2-1, nptu, xptu, yptu, ierrf);
    mncont(ke1, ke2, nptu, xptu, yptu, ierrf);
    if (ierrf < nptu) ierflg = 4;
    if (ierrf == -1)  ierflg = 3;
    return;
//*-*-                                     . . . . . . . . . . jump
L2600:
    step = fWord7[0];
    if (step <= 0) step = 2;
    rno = 0;
    izero = 0;
    for (i = 1; i <= fNpar; ++i) {
	mnrn15(rno, izero);
	rno      = rno*2 - 1;
	fX[i-1] += rno*step*fWerr[i-1];
    }
    mninex(fX);
    mnamin();
    mnrset(0);
    return;
//*-*-                                     . . . . . . . . . . blank line
L3300:
    std::printf(" BLANK COMMAND IGNORED.\n");
    ierflg = 1;
    return;
//*-*  . . . . . . . . obsolete commands     . . . . . . . . . . . . . .
//*-*-                                     . . . . . . . . . . covariance
L3400:
    std::printf(" THE *COVARIANCE* COMMAND IS OSBSOLETE. THE COVARIANCE MATRIX IS NOW SAVED IN A DIFFERENT FORMAT WITH THE *SAVE* COMMAND AND READ IN WITH:*SET COVARIANCE*\n");
    ierflg = 3;
    return;
//*-*-                                       . . . . . . . . . . printout
L3500:
    cneway = "SET PRInt ";
    goto L3100;
//*-*-                                       . . . . . . . . . . gradient
L3600:
    cneway = "SET GRAd  ";
    goto L3100;
//*-*-                                       . . . . . . . . . . matout
L3700:
    cneway = "SHOW COVar";
    goto L3100;
//*-*-                                       . . . . . . . . . error def
L3800:
    cneway = "SET ERRdef";
    goto L3100;
//*-*-                                       . . . . . . . . . . limits
L3900:
    cneway = "SET LIMits";
    goto L3100;
//*-*-                                       . . . . . . . . . . punch
L4000:
    cneway = "SAVE      ";
//*-*-                               ....... come from obsolete commands
L3100:
    std::printf(" OBSOLETE COMMAND:%s   PLEASE USE: %s\n",fCword.c_str()
                                                  ,cneway.c_str());
    fCword = cneway;
    if (fCword == "SAVE      ") goto L1500;
    goto L700;
//*-*                                 . . . . . . . . . . . . . . . . . .
} /* mnexcm_ */

//______________________________________________________________________________
void URMinuit::mnexin(Double_urt *pint)
{
//*-*-*-*-*Transforms the external parameter values U to internal values*-*-*
//*-*      =============================================================
//*-*        Transforms the external parameter values U to internal
//*-*        values in the dense array PINT.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Double_urt pinti;
    Int_urt iint, iext;

    fLimset = kurFALSE;
    for (iint = 1; iint <= fNpar; ++iint) {
	iext = fNexofi[iint-1];
//	mnpint(fU[iext-1], iext-1, pinti);
        mnpint(m_userParameterValue[iext], iext, pinti);
	pint[iint-1] = pinti;
    }
} /* mnexin_ */

//______________________________________________________________________________
void URMinuit::mnfixp(Int_urt iint1, Int_urt &ierr)
{
//*-*-*-*-*-*-*Removes parameter IINT from the internal parameter list*-*-*
//*-*          =======================================================
//*-*        and arranges the rest of the list to fill the hole.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt yyover;
    Int_urt kold, nold, ndex, knew, iext, i, j, m, n, lc, ik;

//*-*-                          first see if it can be done
    ierr = 0;
    Int_urt iint = iint1+1;
    if (iint > fNpar || iint <= 0) {
	ierr = 1;
        std::printf(" MINUIT ERROR.  ARGUMENT TO MNFIXP=%4d\n",iint);
	return;
    }
    iext = fNexofi[iint-1];
    if (fNpfix >= fMaxpar) {
	ierr = 1;
	std::printf(" MINUIT CANNOT FIX PARAMETER%4d MAXIMUM NUMBER THAT CAN BE FIXED IS %d\n",iext,fMaxpar);
	return;
    }
//*-*-                          reduce number of variable parameters by one

//    fNiofex[iext-1] = 0;
    m_userParameterIdToInternalId[iext] = 0;
    nold = fNpar;
    --fNpar;
//*-*-                      save values in case parameter is later restored

    ++fNpfix;
    fIpfix[fNpfix-1]  = iext;
    lc                = iint;
    fXs[fNpfix-1]     = fX[lc-1];
    fXts[fNpfix-1]    = fXt[lc-1];
    fDirins[fNpfix-1] = fWerr[lc-1];
    fGrds[fNpfix-1]   = fGrd[lc-1];
    fG2s[fNpfix-1]    = fG2[lc-1];
    fGsteps[fNpfix-1] = fGstep[lc-1];
//*-*-                       shift values for other parameters to fill hole
    for (ik = iext + 1; ik <= fNu; ++ik) {
//       if (fNiofex[ik-1] > 0) {
       if (m_userParameterIdToInternalId[ik] > 0) {
//          lc = fNiofex[ik-1] - 1;
          lc = m_userParameterIdToInternalId[ik] - 1;
//          fNiofex[ik-1] = lc;
          m_userParameterIdToInternalId[ik] = lc;
          fNexofi[lc-1] = ik;
          fX[lc-1]      = fX[lc];
          fXt[lc-1]     = fXt[lc];
          fDirin[lc-1]  = fDirin[lc];
          fWerr[lc-1]   = fWerr[lc];
          fGrd[lc-1]    = fGrd[lc];
          fG2[lc-1]     = fG2[lc];
          fGstep[lc-1]  = fGstep[lc];
       }
    }
    if (fISW[1] <= 0) return;
//*-*-                   remove one row and one column from variance matrix
    if (fNpar <= 0)   return;
    for (i = 1; i <= nold; ++i) {
	m       = URMath::Max(i,iint);
	n       = URMath::Min(i,iint);
	ndex    = m*(m-1) / 2 + n;
	fFIXPyy[i-1] = fVhmat[ndex-1];
    }
    yyover = 1 / fFIXPyy[iint-1];
    knew   = 0;
    kold   = 0;
    for (i = 1; i <= nold; ++i) {
	for (j = 1; j <= i; ++j) {
	    ++kold;
	    if (j == iint || i == iint) continue;
	    ++knew;
	    fVhmat[knew-1] = fVhmat[kold-1] - fFIXPyy[j-1]*fFIXPyy[i-1]*yyover;
	}
    }
} /* mnfixp_ */

//______________________________________________________________________________
void URMinuit::mnfree(Int_urt k)
{
//*-*-*-*Restores one or more fixed parameter(s) to variable status*-*-*-*-*-*
//*-*    ==========================================================
//*-*        Restores one or more fixed parameter(s) to variable status
//*-*        by inserting it into the internal parameter list at the
//*-*        appropriate place.
//*-*
//*-*        K = 0 means restore all parameters
//*-*        K = 1 means restore the last parameter fixed
//*-*        K = -I means restore external parameter I (if possible)
//*-*        IQ = fix-location where internal parameters were stored
//*-*        IR = external number of parameter being restored
//*-*        IS = internal number of parameter being restored
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt grdv, xv, dirinv, g2v, gstepv, xtv;
    Int_urt i, ipsav, ka, lc, ik, iq, ir, is;

    if (k > 1) {
	std::printf(" CALL TO MNFREE IGNORED.  ARGUMENT GREATER THAN ONE\n");
    }
    if (fNpfix < 1) {
	std::printf(" CALL TO MNFREE IGNORED.  THERE ARE NO FIXED PARAMETERS\n");
    }
    if (k == 1 || k == 0) goto L40;

//*-*-                  release parameter with specified external number
    ka = abs(k);
//    if (fNiofex[ka-1] == 0) goto L15;
    if (m_userParameterIdToInternalId[ka] == 0) goto L15;
    std::printf(" IGNORED.  PARAMETER SPECIFIED IS ALREADY VARIABLE.\n");
    return;
L15:
    if (fNpfix < 1) goto L21;
    for (ik = 1; ik <= fNpfix; ++ik) { if (fIpfix[ik-1] == ka) goto L24; }
L21:
    std::printf(" PARAMETER%4d NOT FIXED.  CANNOT BE RELEASED.\n",ka);
    return;
L24:
    if (ik == fNpfix) goto L40;

//*-*-                  move specified parameter to end of list
    ipsav  = ka;
    xv     = fXs[ik-1];
    xtv    = fXts[ik-1];
    dirinv = fDirins[ik-1];
    grdv   = fGrds[ik-1];
    g2v    = fG2s[ik-1];
    gstepv = fGsteps[ik-1];
    for (i = ik + 1; i <= fNpfix; ++i) {
	fIpfix[i-2]  = fIpfix[i-1];
	fXs[i-2]     = fXs[i-1];
	fXts[i-2]    = fXts[i-1];
	fDirins[i-2] = fDirins[i-1];
	fGrds[i-2]   = fGrds[i-1];
	fG2s[i-2]    = fG2s[i-1];
	fGsteps[i-2] = fGsteps[i-1];
    }
    fIpfix[fNpfix-1]  = ipsav;
    fXs[fNpfix-1]     = xv;
    fXts[fNpfix-1]    = xtv;
    fDirins[fNpfix-1] = dirinv;
    fGrds[fNpfix-1]   = grdv;
    fG2s[fNpfix-1]    = g2v;
    fGsteps[fNpfix-1] = gstepv;
//*-*-               restore last parameter in fixed list  -- IPFIX(NPFIX)
L40:
    if (fNpfix < 1) goto L300;
    ir = fIpfix[fNpfix-1];
    is = 0;
    // shift internal one higher to create hole for parameter that has been freed
    for (ik = fNu; ik >= ir; --ik) {
//       if (fNiofex[ik-1] > 0) {
       if (m_userParameterIdToInternalId[ik] > 0) {
//	    lc = fNiofex[ik-1] + 1;
          lc = m_userParameterIdToInternalId[ik] + 1;
          is = lc - 1;
//	    fNiofex[ik-1] = lc;
          m_userParameterIdToInternalId[ik] = lc;
          fNexofi[lc-1] = ik;
          fX[lc-1]      = fX[lc-2];
          fXt[lc-1]     = fXt[lc-2];
          fDirin[lc-1]  = fDirin[lc-2];
          fWerr[lc-1]   = fWerr[lc-2];
          fGrd[lc-1]    = fGrd[lc-2];
          fG2[lc-1]     = fG2[lc-2];
          fGstep[lc-1]  = fGstep[lc-2];
       }
    }
    ++fNpar;
    if (is == 0) is = fNpar;
//    fNiofex[ir-1] = is;
    m_userParameterIdToInternalId[ir] = is;
    fNexofi[is-1] = ir;
    iq           = fNpfix;
    fX[is-1]     = fXs[iq-1];
    fXt[is-1]    = fXts[iq-1];
    fDirin[is-1] = fDirins[iq-1];
    fWerr[is-1]  = fDirins[iq-1];
    fGrd[is-1]   = fGrds[iq-1];
    fG2[is-1]    = fG2s[iq-1];
    fGstep[is-1] = fGsteps[iq-1];
    --fNpfix;
    fISW[1] = 0;
    fDcovar = 1;
    if (fISW[4] - fItaur >= 1) {
//      std::printf("                    PARAMETER%4d  %s RESTORED TO VARIABLE.\n",ir,
//            fCpnam[ir-1].c_str());
	std::printf("                    PARAMETER%4d  %s RESTORED TO VARIABLE.\n",ir,
                      m_userParameterName[ir].c_str());
    }
    if (k == 0) goto L40;
L300:
//*-*-        if different from internal, external values are taken
    mnexin(fX);
} /* mnfree_ */

//______________________________________________________________________________
void URMinuit::mngrad()
{
//*-*-*-*-*-*-*-*-*-*Interprets the SET GRAD command*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                ===============================
//*-*       Called from MNSET
//*-*       Interprets the SET GRAD command, which informs MINUIT whether
//*-*       the first derivatives of FCN will be calculated by the user
//*-*       inside FCN.  It can check the user derivative calculation
//*-*       by comparing it with a finite difference approximation.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt fzero, err;
    Int_urt i, nparx, lc, istsav;
    Bool_urt lnone;
    static string cwd = "    ";

    fISW[2] = 1;
    nparx   = fNpar;
    if (fWord7[0] > 0) goto L2000;

//*-*-                 get user-calculated first derivatives from FCN
    for (i = 1; i <= fNu; ++i) { fGin[i-1] = fUndefi; }
    mninex(fX);
//    Eval(nparx, fGin, fzero, fU, 2);    ++fNfcn;
    Eval(nparx, fGin, fzero, m_userParameterValue, 2);    ++fNfcn;
    mnderi();
    for (i = 1; i <= fNpar; ++i) { fGRADgf[i-1] = fGrd[i-1]; }
//*-*-                   get MINUIT-calculated first derivatives
    fISW[2] = 0;
    istsav  = fIstrat;
    fIstrat = 2;
    mnhes1();
    fIstrat = istsav;
    std::printf(" CHECK OF GRADIENT CALCULATION IN FCN\n");
    std::printf("            PARAMETER      G(IN FCN)   G(MINUIT)  DG(MINUIT)   AGREEMENT\n");
    fISW[2] = 1;
    lnone = kurFALSE;
    for (lc = 1; lc <= fNpar; ++lc) {
	i   = fNexofi[lc-1];
	cwd = "GOOD";
	err = fDgrd[lc-1];
	if (URMath::Abs(fGRADgf[lc-1] - fGrd[lc-1]) > err)  cwd = " BAD";
	if (fGin[i-1] == fUndefi) {
	    cwd      = "NONE";
	    lnone    = kurTRUE;
	    fGRADgf[lc-1] = 0;
	}
	if (cwd != "GOOD") fISW[2] = 0;
//         std::printf("       %5d  %10s%12.4e%12.4e%12.4e    %s\n",i
//            ,(fCpnam[i-1]).c_str()
//            ,fGRADgf[lc-1],fGrd[lc-1],err,cwd.c_str());
           std::printf("       %5d  %10s%12.4e%12.4e%12.4e    %s\n",i
                    ,(m_userParameterName[i]).c_str()
                    ,fGRADgf[lc-1],fGrd[lc-1],err,cwd.c_str());
    }
    if (lnone) {
	std::printf("  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DERIVATIVE\n");
    }
    if (fISW[2] == 0) {
	std::printf(" MINUIT DOES NOT ACCEPT DERIVATIVE CALCULATIONS BY FCN\n");
        std::printf(" TO FORCE ACCEPTANCE, ENTER *SET GRAD    1*\n");
    }

L2000:
    return;
} /* mngrad_ */

//______________________________________________________________________________
void URMinuit::mnhelp()
{
   string comd = "";
   mnhelp(comd);
}

//______________________________________________________________________________
void URMinuit::mnhelp(const string& comd_in)
{
//*-*-*-*-*-*-*-*HELP routine for MINUIT interactive commands*-*-*-*-*-*-*-*-*
//*-*            ============================================
//*-*
//*-*      COMD ='*' or "" prints a global help for all commands
//*-*      COMD =Command_name: print detailed help for one command.
//*-*         Note that at least 3 characters must be given for the command
//*-*         name.
//*-*
//*-*     Author: Rene Brun
//*-*             comments extracted from the MINUIT documentation file.
//*-*
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//*-*.......................................................................
//*-*
//*-*  Global HELP: Summary of all commands
//*-*  ====================================
//*-*
   string comd = comd_in;
   int (*pf)(int)=toupper; 
    transform(comd.begin(), comd.end(), comd.begin(), pf); 
//    if( comd.length() == 0 || comd[0] == "*" || comd[0] == "?" || comd[0] == 0 || comd=="HELP" ) {
    if( comd.length() == 0 || comd.substr(0,1) == "*" || comd.substr(0,1) == "?" || comd=="HELP" ) {
       std::printf("   ==>List of MINUIT Interactive commands:\n");
       std::printf(" CLEar     Reset all parameter names and values undefined\n");
       std::printf(" CONtour   Make contour map of the user function\n");
       std::printf(" EXIT      Exit from Interactive Minuit\n");
       std::printf(" FIX       Cause parameter(s) to remain constant\n");
       std::printf(" HESse     Calculate the Hessian or error matrix.\n");
       std::printf(" IMPROVE   Search for a new minimum around current minimum\n");
       std::printf(" MIGrad    Minimize by the method of Migrad\n");
       std::printf(" MINImize  MIGRAD + SIMPLEX method if Migrad fails\n");
       std::printf(" MINOs     Exact (non-linear) parameter error analysis\n");
       std::printf(" MNContour Calculate one MINOS function contour\n");
       std::printf(" PARameter Define or redefine new parameters and values\n");
       std::printf(" RELease   Make previously FIXed parameters variable again\n");
       std::printf(" REStore   Release last parameter fixed\n");
       std::printf(" SAVe      Save current parameter values on a file\n");
       std::printf(" SCAn      Scan the user function by varying parameters\n");
       std::printf(" SEEk      Minimize by the method of Monte Carlo\n");
       std::printf(" SET       Set various MINUIT constants or conditions\n");
       std::printf(" SHOw      Show values of current constants or conditions\n");
       std::printf(" SIMplex   Minimize by the method of Simplex\n");
	goto L99;
    }

//*-* __________________________________________________________________
//*-*
//*-* --  Command CLEAR
//*-* --  =============
//*-*
    if( !strncmp(comd.c_str(),"CLE",3) ) {
       std::printf(" ***>CLEAR");
       std::printf(" Resets all parameter names and values to undefined.\n");
       std::printf(" Must normally be followed by a PARameters command or \n");
       std::printf(" equivalent, in order to define parameter values.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command CONTOUR
//*-* --  ===============
//*-* .
    if( !strncmp(comd.c_str(),"CON",3) ) {
       std::printf(" ***>CONTOUR <par1>  <par2>  [devs]  [ngrid]\n");
       std::printf(" Instructs Minuit to trace contour lines of the user function\n");
       std::printf(" with respect to the two parameters whose external numbers\n");
       std::printf(" are <par1> and <par2>.\n");
       std::printf(" Other variable parameters of the function, if any, will have\n");
       std::printf(" their values fixed at the current values during the contour\n");
       std::printf(" tracing. The optional parameter [devs] (default value 2.)\n");
       std::printf(" gives the number of standard deviations in each parameter\n");
       std::printf(" which should lie entirely within the plotting area.\n");
       std::printf(" Optional parameter [ngrid] (default value 25 unless page\n");
       std::printf(" size is too small) determines the resolution of the plot,\n");
       std::printf(" i.e. the number of rows and columns of the grid at which the\n");
       std::printf(" function will be evaluated. [See also MNContour.]\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command END
//*-* --  ===========
//*-* .
    if( !strncmp(comd.c_str(),"END",3) ) {
       std::printf(" ***>END\n");
       std::printf(" Signals the end of a data block (i.e., the end of a fit),\n");
       std::printf(" and implies that execution should continue, because another\n");
       std::printf(" Data Block follows. A Data Block is a set of Minuit data\n");
       std::printf(" consisting of\n");
       std::printf("     (1) A Title,\n");
       std::printf("     (2) One or more Parameter Definitions,\n");
       std::printf("     (3) A blank line, and\n");
       std::printf("     (4) A set of Minuit Commands.\n");
       std::printf(" The END command is used when more than one Data Block is to\n");
       std::printf(" be used with the same FCN function. It first causes Minuit\n");
       std::printf(" to issue a CALL FCN with IFLAG=3, in order to allow FCN to\n");
       std::printf(" perform any calculations associated with the final fitted\n");
       std::printf(" parameter values, unless a CALL FCN 3 command has already\n");
       std::printf(" been executed at the current FCN value.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* .
//*-* --
//*-* --  Command EXIT
//*-* --  ============
    if( !strncmp(comd.c_str(),"EXI",3) ) {
       std::printf(" ***>EXIT\n");
       std::printf(" Signals the end of execution.\n");
       std::printf(" The EXIT command first causes Minuit to issue a CALL FCN\n");
       std::printf(" with IFLAG=3, to allow FCN to perform any calculations\n");
       std::printf(" associated with the final fitted parameter values, unless a\n");
       std::printf(" CALL FCN 3 command has already been executed.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command FIX
//*-* --  ===========
//*-* .
    if( !strncmp(comd.c_str(),"FIX",3) ) {
       std::printf(" ***>FIX} <parno> [parno] ... [parno]\n");
       std::printf(" Causes parameter(s) <parno> to be removed from the list of\n");
       std::printf(" variable parameters, and their value(s) will remain constant\n");
       std::printf(" during subsequent minimizations, etc., until another command\n");
       std::printf(" changes their value(s) or status.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command HESSE
//*-* --  =============
//*-* .
    if( !strncmp(comd.c_str(),"HES",3) ) {
       std::printf(" ***>HESse  [maxcalls]");
       std::printf(" Calculate, by finite differences, the Hessian or error matrix.\n");
       std::printf("  That is, it calculates the full matrix of second derivatives\n");
       std::printf(" of the function with respect to the currently variable\n");
       std::printf(" parameters, and inverts it, printing out the resulting error\n");
       std::printf(" matrix. The optional argument [maxcalls] specifies the\n");
       std::printf(" (approximate) maximum number of function calls after which\n");
       std::printf(" the calculation will be stopped.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command IMPROVE
//*-* --  ===============
//*-* .
    if( !strncmp(comd.c_str(),"IMP",3) ) {
       std::printf(" ***>IMPROVE  [maxcalls]");
       std::printf(" If a previous minimization has converged, and the current\n");
       std::printf(" values of the parameters therefore correspond to a local\n");
       std::printf(" minimum of the function, this command requests a search for\n");
       std::printf(" additional distinct local minima.\n");
       std::printf(" The optional argument [maxcalls] specifies the (approximate\n");
       std::printf(" maximum number of function calls after which the calculation\n");
       std::printf(" will be stopped.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command MIGRAD
//*-* --  ==============
//*-* .
    if( !strncmp(comd.c_str(),"MIG",3) ) {
       std::printf(" ***>MIGrad  [maxcalls]  [tolerance]\n");
       std::printf(" Causes minimization of the function by the method of Migrad,\n");
       std::printf(" the most efficient and complete single method, recommended\n");
       std::printf(" for general functions (see also MINImize).\n");
       std::printf(" The minimization produces as a by-product the error matrix\n");
       std::printf(" of the parameters, which is usually reliable unless warning\n");
       std::printf(" messages are produced.\n");
       std::printf(" The optional argument [maxcalls] specifies the (approximate)\n");
       std::printf(" maximum number of function calls after which the calculation\n");
       std::printf(" will be stopped even if it has not yet converged.\n");
       std::printf(" The optional argument [tolerance] specifies required tolerance\n");
       std::printf(" on the function value at the minimum.\n");
       std::printf(" The default tolerance is 0.1, and the minimization will stop\n");
       std::printf(" when the estimated vertical distance to the minimum (EDM) is\n");
       std::printf(" less than 0.001*[tolerance]*UP (see [SET ERRordef]).\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command MINIMIZE
//*-* --  ================
//*-* .
    if( !strncmp(comd.c_str(),"MINI",4) ) {
       std::printf(" ***>MINImize  [maxcalls] [tolerance]\n");
       std::printf(" Causes minimization of the function by the method of Migrad,\n");
       std::printf(" as does the MIGrad command, but switches to the SIMplex method\n");
       std::printf(" if Migrad fails to converge. Arguments are as for MIGrad.\n");
       std::printf(" Note that command requires four characters to be unambiguous.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command MINOS
//*-* --  =============
//*-* .
    if( !strncmp(comd.c_str(),"MIN0",4) ) {
       std::printf(" ***>MINOs  [maxcalls]  [parno] [parno] ...\n");
       std::printf(" Causes a Minos error analysis to be performed on the parameters\n");
       std::printf(" whose numbers [parno] are specified. If none are specified,\n");
       std::printf(" Minos errors are calculated for all variable parameters.\n");
       std::printf(" Minos errors may be expensive to calculate, but are very\n");
       std::printf(" reliable since they take account of non-linearities in the\n");
       std::printf(" problem as well as parameter correlations, and are in general\n");
       std::printf(" asymmetric.\n");
       std::printf(" The optional argument [maxcalls] specifies the (approximate)\n");
       std::printf(" maximum number of function calls per parameter requested,\n");
       std::printf(" after which the calculation will stop for that parameter.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command MNCONTOUR
//*-* --  =================
//*-* .
    if( !strncmp(comd.c_str(),"MNC",3) ) {
       std::printf(" ***>MNContour  <par1> <par2> [npts]\n");
       std::printf(" Calculates one function contour of FCN with respect to\n");
       std::printf(" parameters par1 and par2, with FCN minimized always with\n");
       std::printf(" respect to all other NPAR-2 variable parameters (if any).\n");
       std::printf(" Minuit will try to find npts points on the contour (default 20)\n");
       std::printf(" If only two parameters are variable at the time, it is not\n");
       std::printf(" necessary to specify their numbers. To calculate more than\n");
       std::printf(" one contour, it is necessary to SET ERRordef to the appropriate\n");
       std::printf(" value and issue the MNContour command for each contour.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command PARAMETER
//*-* --  =================
//*-* .
    if( !strncmp(comd.c_str(),"PAR",3) ) {
       std::printf(" ***>PARameters\n");
       std::printf(" followed by one or more parameter definitions.\n");
       std::printf(" Parameter definitions are of the form:\n");
       std::printf("   <number>  ''name''  <value>  <step>  [lolim] [uplim] \n");
       std::printf(" for example:\n");
       std::printf("  3  ''K width''  1.2   0.1\n");
       std::printf(" the last definition is followed by a blank line or a zero.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command RELEASE
//*-* --  ===============
//*-* .
    if( !strncmp(comd.c_str(),"REL",3) ) {
       std::printf(" ***>RELease  <parno> [parno] ... [parno]\n");
       std::printf(" If <parno> is the number of a previously variable parameter\n");
       std::printf(" which has been fixed by a command: FIX <parno>, then that\n");
       std::printf(" parameter will return to variable status.  Otherwise a warning\n");
       std::printf(" message is printed and the command is ignored.\n");
       std::printf(" Note that this command operates only on parameters which were\n");
       std::printf(" at one time variable and have been FIXed. It cannot make\n");
       std::printf(" constant parameters variable; that must be done by redefining\n");
       std::printf(" the parameter with a PARameters command.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command RESTORE
//*-* --  ===============
//*-* .
    if( !strncmp(comd.c_str(),"RES",3) ) {
       std::printf(" ***>REStore  [code]\n");
       std::printf(" If no [code] is specified, this command restores all previously\n");
       std::printf(" FIXed parameters to variable status. If [code]=1, then only\n");
       std::printf(" the last parameter FIXed is restored to variable status.\n");
       std::printf(" If code is neither zero nor one, the command is ignored.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command RETURN
//*-* --  ==============
//*-* .
    if( !strncmp(comd.c_str(),"RET",3) ) {
       std::printf(" ***>RETURN\n");
       std::printf(" Signals the end of a data block, and instructs Minuit to return\n");
       std::printf(" to the program which called it. The RETurn command first\n");
       std::printf(" causes Minuit to CALL FCN with IFLAG=3, in order to allow FCN\n");
       std::printf(" to perform any calculations associated with the final fitted\n");
       std::printf(" parameter values, unless a CALL FCN 3 command has already been\n");
       std::printf(" executed at the current FCN value.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SAVE
//*-* --  ============
//*-* .
    if( !strncmp(comd.c_str(),"SAV",3) ) {
       std::printf(" ***>SAVe\n");
       std::printf(" Causes the current parameter values to be saved on a file in\n");
       std::printf(" such a format that they can be read in again as Minuit\n");
       std::printf(" parameter definitions. If the covariance matrix exists, it is\n");
       std::printf(" also output in such a format. The unit number is by default 7,\n");
       std::printf(" or that specified by the user in his call to MINTIO or\n");
       std::printf(" MNINIT. The user is responsible for opening the file previous\n");
       std::printf(" to issuing the [SAVe] command (except where this can be done\n");
       std::printf(" interactively).\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SCAN
//*-* --  ============
//*-* .
    if( !strncmp(comd.c_str(),"SCA",3) ) {
       std::printf(" ***>SCAn  [parno]  [numpts] [from]  [to]\n");
       std::printf(" Scans the value of the user function by varying parameter\n");
       std::printf(" number [parno], leaving all other parameters fixed at the\n");
       std::printf(" current value. If [parno] is not specified, all variable\n");
       std::printf(" parameters are scanned in sequence.\n");
       std::printf(" The number of points [numpts] in the scan is 40 by default,\n");
       std::printf(" and cannot exceed 100. The range of the scan is by default\n");
       std::printf(" 2 standard deviations on each side of the current best value,\n");
       std::printf(" but can be specified as from [from] to [to].\n");
       std::printf(" After each scan, if a new minimum is found, the best parameter\n");
       std::printf(" values are retained as start values for future scans or\n");
       std::printf(" minimizations. The curve resulting from each scan is plotted\n");
       std::printf(" on the output unit in order to show the approximate behaviour\n");
       std::printf(" of the function.\n");
       std::printf(" This command is not intended for minimization, but is sometimes\n");
       std::printf(" useful for debugging the user function or finding a\n");
       std::printf(" reasonable starting point.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SEEK
//*-* --  ============
//*-* .
    if( !strncmp(comd.c_str(),"SEE",3) ) {
       std::printf(" ***>SEEk  [maxcalls]  [devs]\n");
       std::printf(" Causes a Monte Carlo minimization of the function, by choosing\n");
       std::printf(" random values of the variable parameters, chosen uniformly\n");
       std::printf(" over a hypercube centered at the current best value.\n");
       std::printf(" The region size is by default 3 standard deviations on each\n");
       std::printf(" side, but can be changed by specifying the value of [devs].\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SET
//*-* --  ===========
//*-* .
    if( !strncmp(comd.c_str(),"SET",3) ) {
       std::printf(" ***>SET <option_name>\n");
       std::printf("  SET BATch\n");
       std::printf("    Informs Minuit that it is running in batch mode.\n");

       std::printf(" \n");
       std::printf("  SET EPSmachine  <accuracy>\n");
       std::printf("    Informs Minuit that the relative floating point arithmetic\n");
       std::printf("    precision is <accuracy>. Minuit determines the nominal\n");
       std::printf("    precision itself, but the SET EPSmachine command can be\n");
       std::printf("    used to override Minuit own determination, when the user\n");
       std::printf("    knows that the FCN function value is not calculated to\n");
       std::printf("    the nominal machine accuracy. Typical values of <accuracy>\n");
       std::printf("    are between 10**-5 and 10**-14.\n");

       std::printf(" \n");
       std::printf("  SET ERRordef  <up>\n");
       std::printf("    Sets the value of UP (default value= 1.), defining parameter\n");
       std::printf("    errors. Minuit defines parameter errors as the change\n");
       std::printf("    in parameter value required to change the function value\n");
       std::printf("    by UP. Normally, for chisquared fits UP=1, and for negative\n");
       std::printf("    log likelihood, UP=0.5.\n");

       std::printf(" ");
       std::printf("   SET GRAdient  [force]\n");
       std::printf("    Informs Minuit that the user function is prepared to\n");
       std::printf("    calculate its own first derivatives and return their values\n");
       std::printf("    in the array GRAD when IFLAG=2 (see specs of FCN).\n");
       std::printf("    If [force] is not specified, Minuit will calculate\n");
       std::printf("    the FCN derivatives by finite differences at the current\n");
       std::printf("    point and compare with the user calculation at that point,\n");
       std::printf("    accepting the user values only if they agree.\n");
       std::printf("    If [force]=1, Minuit does not do its own derivative\n");
       std::printf("    calculation, and uses the derivatives calculated in FCN.\n");

       std::printf(" \n");
       std::printf("   SET INPut  [unitno]  [filename]\n");
       std::printf("    Causes Minuit, in data-driven mode only, to read subsequent\n");
       std::printf("    commands (or parameter definitions) from a different input\n");
       std::printf("    file. If no [unitno] is specified, reading reverts to the\n");
       std::printf("    previous input file, assuming that there was one.\n");
       std::printf("    If [unitno] is specified, and that unit has not been opened,\n");
       std::printf("    then Minuit attempts to open the file [filename]} if a\n");
       std::printf("    name is specified. If running in interactive mode and\n");
       std::printf("    [filename] is not specified and [unitno] is not opened,\n");
       std::printf("    Minuit prompts the user to enter a file name.\n");
       std::printf("    If the word REWIND is added to the command (note:no blanks\n");
       std::printf("    between INPUT and REWIND), the file is rewound before\n");
       std::printf("    reading. Note that this command is implemented in standard\n");
       std::printf("    Fortran 77 and the results may depend on the  system;\n");
       std::printf("    for example, if a filename is given under VM/CMS, it must\n");
       std::printf("    be preceeded by a slash.\n");

       std::printf(" \n");
       std::printf("   SET INTeractive\n");
       std::printf("    Informs Minuit that it is running interactively.\n");

       std::printf(" \n");
       std::printf("   SET LIMits  [parno]  [lolim]  [uplim]\n");
       std::printf("    Allows the user to change the limits on one or all\n");
       std::printf("    parameters. If no arguments are specified, all limits are\n");
       std::printf("    removed from all parameters. If [parno] alone is specified,\n");
       std::printf("    limits are removed from parameter [parno].\n");
       std::printf("    If all arguments are specified, then parameter [parno] will\n");
       std::printf("    be bounded between [lolim] and [uplim].\n");
       std::printf("    Limits can be specified in either order, Minuit will take\n");
       std::printf("    the smaller as [lolim] and the larger as [uplim].\n");
       std::printf("    However, if [lolim] is equal to [uplim], an error condition\n");
       std::printf("    results.\n");

       std::printf(" \n");
       std::printf("   SET LINesperpage\n");
       std::printf("     Sets the number of lines for one page of output.\n");
       std::printf("     Default value is 24 for interactive mode\n");

       std::printf(" \n");
       std::printf("   SET NOGradient");
       std::printf("    The inverse of SET GRAdient, instructs Minuit not to\n");
       std::printf("    use the first derivatives calculated by the user in FCN.\n");

       std::printf(" \n");
       std::printf("   SET NOWarnings\n");
       std::printf("    Supresses Minuit warning messages.\n");

       std::printf(" \n");
       std::printf("   SET OUTputfile  <unitno>\n");
       std::printf("    Instructs Minuit to write further output to unit <unitno>.\n");

       std::printf(" \n");
       std::printf("   SET PAGethrow  <integer>\n");
       std::printf("    Sets the carriage control character for ``new page'' to\n");
       std::printf("    <integer>. Thus the value 1 produces a new page, and 0\n");
       std::printf("    produces a blank line, on some devices (see TOPofpage)\n");


       std::printf(" \n");
       std::printf("   SET PARameter  <parno>  <value>\n");
       std::printf("    Sets the value of parameter <parno> to <value>.\n");
       std::printf("    The parameter in question may be variable, fixed, or\n");
       std::printf("    constant, but must be defined.\n");

       std::printf(" \n");
       std::printf("   SET PRIntout  <level>\n");
       std::printf("    Sets the print level, determining how much output will be\n");
       std::printf("    produced. Allowed values and their meanings are displayed\n");
       std::printf("    after a SHOw PRInt command, and are currently <level>=:\n");
       std::printf("      [-1]  no output except from SHOW commands\n");
       std::printf("       [0]  minimum output\n");
       std::printf("       [1]  default value, normal output\n");
       std::printf("       [2]  additional output giving intermediate results.\n");
       std::printf("       [3]  maximum output, showing progress of minimizations.\n");
       std::printf("    Note: See also the SET WARnings command.\n");

       std::printf(" \n");
       std::printf("   SET RANdomgenerator  <seed>\n");
       std::printf("    Sets the seed of the random number generator used in SEEk.\n");
       std::printf("    This can be any integer between 10000 and 900000000, for\n");
       std::printf("    example one which was output from a SHOw RANdom command of\n");
       std::printf("    a previous run.\n");

       std::printf(" \n");
       std::printf("   SET STRategy  <level>\n");
       std::printf("    Sets the strategy to be used in calculating first and second\n");
       std::printf("    derivatives and in certain minimization methods.\n");
       std::printf("    In general, low values of <level> mean fewer function calls\n");
       std::printf("    and high values mean more reliable minimization.\n");
       std::printf("    Currently allowed values are 0, 1 (default), and 2.\n");

       std::printf(" \n");
       std::printf("   SET TITle\n");
       std::printf("    Informs Minuit that the next input line is to be considered\n");
       std::printf("    the (new) title for this task or sub-task.  This is for\n");
       std::printf("    the convenience of the user in reading his output.\n");

       std::printf(" \n");
       std::printf("   SET WARnings\n");
       std::printf("    Instructs Minuit to output warning messages when suspicious\n");
       std::printf("    conditions arise which may indicate unreliable results.\n");
       std::printf("    This is the default.\n");

       std::printf(" \n");
       std::printf("    SET WIDthpage\n");
       std::printf("    Informs Minuit of the output page width.\n");
       std::printf("    Default values are 80 for interactive jobs\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SHOW
//*-* --  ============
//*-* .
    if( !strncmp(comd.c_str(),"SHO",3) ) {
       std::printf(" ***>SHOw  <option_name>\n");
       std::printf("  All SET XXXX commands have a corresponding SHOw XXXX command.\n");
       std::printf("  In addition, the SHOw commands listed starting here have no\n");
       std::printf("  corresponding SET command for obvious reasons.\n");

       std::printf(" \n");
       std::printf("   SHOw CORrelations\n");
       std::printf("    Calculates and prints the parameter correlations from the\n");
       std::printf("    error matrix.\n");

       std::printf(" \n");
       std::printf("   SHOw COVariance\n");
       std::printf("    Prints the (external) covariance (error) matrix.\n");

       std::printf(" \n");
       std::printf("   SHOw EIGenvalues\n");
       std::printf("    Calculates and prints the eigenvalues of the covariance\n");
       std::printf("    matrix.");

       std::printf(" \n");
       std::printf("   SHOw FCNvalue\n");
       std::printf("    Prints the current value of FCN.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command SIMPLEX
//*-* --  ===============
//*-* .
    if( !strncmp(comd.c_str(),"SIM",3) ) {
       std::printf(" ***>SIMplex  [maxcalls]  [tolerance]\n");
       std::printf(" Performs a function minimization using the simplex method of\n");
       std::printf(" Nelder and Mead. Minimization terminates either when the\n");
       std::printf(" function has been called (approximately) [maxcalls] times,\n");
       std::printf(" or when the estimated vertical distance to minimum (EDM) is\n");
       std::printf(" less than [tolerance].\n");
       std::printf(" The default value of [tolerance] is 0.1*UP(see SET ERRordef).\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command STANDARD
//*-* --  ================
//*-* .
    if( !strncmp(comd.c_str(),"STA",3) ) {
       std::printf(" ***>STAndard\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command STOP
//*-* --  ============
//*-* .
    if( !strncmp(comd.c_str(),"STO",3) ) {
       std::printf(" ***>STOP\n");
       std::printf(" Same as EXIT.\n");
	goto L99;
    }
//*-* __________________________________________________________________
//*-* --
//*-* --  Command TOPOFPAGE
//*-* --  =================
//*-* .
    if( !strncmp(comd.c_str(),"TOP",3) ) {
       std::printf(" ***>TOPofpage\n");
       std::printf(" Causes Minuit to write the character specified in a\n");
       std::printf(" SET PAGethrow command (default = 1) to column 1 of the output\n");
       std::printf(" file, which may or may not position your output medium to\n");
       std::printf(" the top of a page depending on the device and system.\n");
	goto L99;
    }
//*-* __________________________________________________________________
       std::printf(" Unknown MINUIT command. Type HELP for list of commands.\n");

L99:
    return;
} /* mnhelp_ */

//______________________________________________________________________________
void URMinuit::mnhess()
{
//*-*-*-*-*-*Calculates the full second-derivative matrix of FCN*-*-*-*-*-*-*-*
//*-*        ===================================================
//*-*        by taking finite differences. When calculating diagonal
//*-*        elements, it may iterate so that step size is nearly that
//*-*        which gives function change= UP/10. The first derivatives
//*-*        of course come as a free side effect, but with a smaller
//*-*        step size in order to obtain a known accuracy.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dmin_, dxdi, elem, wint, tlrg2, d, dlast, ztemp, g2bfor;
    Double_urt df, aimsag, fs1, tlrstp, fs2, stpinm, g2i, sag=0, xtf, xti, xtj;
    Int_urt icyc, ncyc, ndex, idrv, iext, npar2, i, j, ifail, npard, nparx, id, multpy;
    Bool_urt ldebug;

    ldebug = fIdbg[3] >= 1;
    if (fAmin == fUndefi) {
	mnamin();
    }
    if (fIstrat <= 0) {
	ncyc   = 3;
	tlrstp = .5;
	tlrg2  = .1;
    } else if (fIstrat == 1) {
	ncyc   = 5;
	tlrstp = .3;
	tlrg2  = .05;
    } else {
	ncyc   = 7;
	tlrstp = .1;
	tlrg2  = .02;
    }
    if (fISW[4] >= 2 || ldebug) {
	std::printf("   START COVARIANCE MATRIX CALCULATION.\n");
    }
    fCfrom  = "HESSE   ";
    fNfcnfr = fNfcn;
    fCstatu = "OK        ";
    npard   = fNpar;
//*-*-                make sure starting at the right place
    mninex(fX);
    nparx = fNpar;
//    Eval(nparx, fGin, fs1, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, fs1, m_userParameterValue, 4);    ++fNfcn;
    if (fs1 != fAmin) {
        df    = fAmin - fs1;
       ostringstream warning;
       warning << "function value differs from AMIN by " << df;
	mnwarn("D", "MNHESS", warning.str().c_str());
    }
    fAmin = fs1;
    if (ldebug) {
	std::printf(" PAR D   GSTEP           D          G2         GRD         SAG    \n");
    }
//*-*-                                       . . . . . . diagonal elements .

//*-*-        ISW(2) = 1 if approx, 2 if not posdef, 3 if ok
//*-*-        AIMSAG is the sagitta we are aiming for in second deriv calc.

    aimsag = URMath::Sqrt(fEpsma2)*(URMath::Abs(fAmin) + fUp);
//*-*-        Zero the second derivative matrix
    npar2 = fNpar*(fNpar + 1) / 2;
    for (i = 1; i <= npar2; ++i) { fVhmat[i-1] = 0; }

//*-*-        Loop over variable parameters for second derivatives
    idrv = 2;
    for (id = 1; id <= npard; ++id) {
	i = id + fNpar - npard;
	iext = fNexofi[i-1];
	if (fG2[i-1] == 0) {
           ostringstream warning;
           warning << "Second derivative enters zero, param " << iext;
	    mnwarn("W", "HESSE", warning.str().c_str());
	    wint = fWerr[i-1];
//            if (fNvarl[iext-1] > 1) {
	    if (m_userParameterFlag[iext] > 1) {
		mndxdi(fX[i-1], i-1, dxdi);
		if (URMath::Abs(dxdi) < .001) wint = .01;
		else                          wint /= URMath::Abs(dxdi);
	    }
	    fG2[i-1] = fUp / (wint*wint);
	}
	xtf   = fX[i-1];
	dmin_ = fEpsma2*8*URMath::Abs(xtf);

//*-*-                              find step which gives sagitta = AIMSAG
	d = URMath::Abs(fGstep[i-1]);
        int skip50 = 0;
	for (icyc = 1; icyc <= ncyc; ++icyc) {
//*-*-                              loop here only if SAG=0
	    for (multpy = 1; multpy <= 5; ++multpy) {
//*-*-          take two steps
		fX[i-1] = xtf + d;
		mninex(fX);
		nparx = fNpar;
//                Eval(nparx, fGin, fs1, fU, 4);    ++fNfcn;
		Eval(nparx, fGin, fs1, m_userParameterValue, 4);    ++fNfcn;
		fX[i-1] = xtf - d;
		mninex(fX);
//		Eval(nparx, fGin, fs2, fU, 4);    ++fNfcn;
                Eval(nparx, fGin, fs2, m_userParameterValue, 4);    ++fNfcn;
		fX[i-1] = xtf;
		sag = (fs1 + fs2 - fAmin*2)*.5;
		if (sag != 0) goto L30;
		if (fGstep[i-1] < 0) {
		    if (d >= .5) goto L26;
		    d *= 10;
		    if (d > .5) 	d = .51;
		    continue;
		}
		d *= 10;
	    }
L26:
           { ostringstream warning3;
            warning3 << "Second derivative zero for parameter " << iext;
            mnwarn("W", "HESSE", warning3.str().c_str()); }
	    goto L390;
//*-*-                            SAG is not zero
L30:
	    g2bfor    = fG2[i-1];
	    fG2[i-1]  = sag*2 / (d*d);
	    fGrd[i-1] = (fs1 - fs2) / (d*2);
	    if (ldebug) {
		std::printf("%4d%2d%12.5g%12.5g%12.5g%12.5g\n",i,idrv,fGstep[i-1],fG2[i-1],fGrd[i-1],sag);
	    }
	    if (fGstep[i-1] > 0) fGstep[i-1] =  URMath::Abs(d);
	    else                 fGstep[i-1] = -URMath::Abs(d);
	    fDirin[i-1] = d;
	    fHESSyy[i-1]= fs1;
	    dlast       = d;
	    d           = URMath::Sqrt(aimsag*2 / URMath::Abs(fG2[i-1]));
//*-*-        if parameter has limits, max int step size = 0.5
	    stpinm = .5;
	    if (fGstep[i-1] < 0) d = URMath::Min(d,stpinm);
	    if (d < dmin_) d = dmin_;
//*-*-          see if converged
	    if (URMath::Abs((d - dlast) / d) < tlrstp ||
	        URMath::Abs((fG2[i-1] - g2bfor) / fG2[i-1]) < tlrg2) {
               skip50 = 1;
               break;
            }
	    d = URMath::Min(d,dlast*102);
	    d = URMath::Max(d,dlast*.1);
	}
//*-*-                      end of step size loop
        if (!skip50) {
           ostringstream warning;
           warning << "Second Deriv. SAG,AIM= " << iext << " " << sag << " " << aimsag;
	   mnwarn("D", "MNHESS", warning.str().c_str());
        }
	ndex = i*(i + 1) / 2;
	fVhmat[ndex-1] = fG2[i-1];
    }
//*-*-                             end of diagonal second derivative loop
    mninex(fX);
//*-*-                                    refine the first derivatives
    if (fIstrat > 0) mnhes1();
    fISW[1] = 3;
    fDcovar = 0;
//*-*-                                       . . . .  off-diagonal elements

    if (fNpar == 1) goto L214;
    for (i = 1; i <= fNpar; ++i) {
	for (j = 1; j <= i-1; ++j) {
	    xti     = fX[i-1];
	    xtj     = fX[j-1];
	    fX[i-1] = xti + fDirin[i-1];
	    fX[j-1] = xtj + fDirin[j-1];
	    mninex(fX);
//	    Eval(nparx, fGin, fs1, fU, 4);	    ++fNfcn;
            Eval(nparx, fGin, fs1, m_userParameterValue, 4);	    ++fNfcn;
	    fX[i-1] = xti;
	    fX[j-1] = xtj;
	    elem = (fs1 + fAmin - fHESSyy[i-1] - fHESSyy[j-1]) / (
		    fDirin[i-1]*fDirin[j-1]);
	    ndex = i*(i-1) / 2 + j;
	    fVhmat[ndex-1] = elem;
	}
    }
L214:
    mninex(fX);
//*-*-                 verify matrix positive-definite
    mnpsdf();
    for (i = 1; i <= fNpar; ++i) {
	for (j = 1; j <= i; ++j) {
	    ndex = i*(i-1) / 2 + j;
      // factor of 0.5 here is based on observation not derivation
      // and noting that factor of 2 appears in setting fVhmat below
      // after inversion of fP (MRS)
      fSecDer[ndex-1] = fVhmat[ndex-1]*.5;
	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1];
	    fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];
	}
    }
    mnvert(fP, fMaxint, fMaxint, fNpar, ifail);
    if (ifail > 0) {
	mnwarn("W", "HESSE", "Matrix inversion fails.");
	goto L390;
    }
//*-*-                                       . . . . . . .  calculate  e d m
    fEDM = 0;

    for (i = 1; i <= fNpar; ++i) {
//*-*-                             off-diagonal elements
	ndex = i*(i-1) / 2;
	for (j = 1; j <= i-1; ++j) {
	    ++ndex;
	    ztemp = fP[i + j*fMaxpar - fMaxpar-1]*2;
	    fEDM += fGrd[i-1]*ztemp*fGrd[j-1];
	    fVhmat[ndex-1] = ztemp;
	}
//*-*-                             diagonal elements
	++ndex;
	fVhmat[ndex-1] = fP[i + i*fMaxpar - fMaxpar-1]*2;
	fEDM += fP[i + i*fMaxpar - fMaxpar-1]*(fGrd[i-1]*fGrd[i-1]);
    }
    if (fISW[4] >= 1 && fISW[1] == 3 && fItaur == 0) {
	std::printf(" COVARIANCE MATRIX CALCULATED SUCCESSFULLY\n");
    }
    goto L900;
//*-*-                             failure to invert 2nd deriv matrix
L390:
    fISW[1] = 1;
    fDcovar = 1;
    fCstatu = "FAILED    ";
    if (fISW[4] >= 0) {
	std::printf("  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. \n");
    }
    for (i = 1; i <= fNpar; ++i) {
	ndex = i*(i-1) / 2;
	for (j = 1; j <= i-1; ++j) {
	    ++ndex;
	    fVhmat[ndex-1] = 0;
	}
	++ndex;
	g2i = fG2[i-1];
	if (g2i <= 0) g2i = 1;
	fVhmat[ndex-1] = 2 / g2i;
    }
L900:
    return;
} /* mnhess_ */

//______________________________________________________________________________
void URMinuit::mnhes1()
{
//*-*-*-*Calculate first derivatives (GRD) and uncertainties (DGRD)*-*-*-*-*-*
//*-*    ==========================================================
//*-*         and appropriate step sizes GSTEP
//*-*      Called from MNHESS and MNGRAD
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dmin_, d, dfmin, dgmin=0, change, chgold, grdold=0, epspri;
    Double_urt fs1, optstp, fs2, grdnew=0, sag, xtf;
    Int_urt icyc, ncyc=0, idrv, i, nparx;
    Bool_urt ldebug;

    ldebug = fIdbg[5] >= 1;
    if (fIstrat <= 0) ncyc = 1;
    if (fIstrat == 1) ncyc = 2;
    if (fIstrat > 1)  ncyc = 6;
    idrv = 1;
    nparx = fNpar;
    dfmin = fEpsma2*4*(URMath::Abs(fAmin) + fUp);
//*-*-                                    main loop over parameters
    for (i = 1; i <= fNpar; ++i) {
	xtf    = fX[i-1];
	dmin_  = fEpsma2*4*URMath::Abs(xtf);
	epspri = fEpsma2 + URMath::Abs(fGrd[i-1]*fEpsma2);
	optstp = URMath::Sqrt(dfmin / (URMath::Abs(fG2[i-1]) + epspri));
	d = URMath::Abs(fGstep[i-1])*.2;
	if (d > optstp) d = optstp;
	if (d < dmin_)  d = dmin_;
	chgold = 1e4;
//*-*-                                      iterate reducing step size
	for (icyc = 1; icyc <= ncyc; ++icyc) {
	    fX[i-1] = xtf + d;
	    mninex(fX);
//	    Eval(nparx, fGin, fs1, fU, 4);	    ++fNfcn;
            Eval(nparx, fGin, fs1, m_userParameterValue, 4);	    ++fNfcn;
	    fX[i-1] = xtf - d;
	    mninex(fX);
//	    Eval(nparx, fGin, fs2, fU, 4);	    ++fNfcn;
            Eval(nparx, fGin, fs2, m_userParameterValue, 4);	    ++fNfcn;
	    fX[i-1] = xtf;
//*-*-                                      check if step sizes appropriate
	    sag    = (fs1 + fs2 - fAmin*2)*.5;
	    grdold = fGrd[i-1];
	    grdnew = (fs1 - fs2) / (d*2);
	    dgmin  = fEpsmac*(URMath::Abs(fs1) + URMath::Abs(fs2)) / d;
	    if (ldebug) {
		std::printf("%4d%2d%12.5g%12.5g%12.5g%12.5g%12.5g\n",i,idrv,fGstep[i-1],d,fG2[i-1],grdnew,sag);
	    }
	    if (grdnew == 0) goto L60;
	    change = URMath::Abs((grdold - grdnew) / grdnew);
	    if (change > chgold && icyc > 1) goto L60;
	    chgold    = change;
	    fGrd[i-1] = grdnew;
	    if (fGstep[i-1] > 0) fGstep[i-1] =  URMath::Abs(d);
	    else                 fGstep[i-1] = -URMath::Abs(d);
//*-*-                 decrease step until first derivative changes by <5%
	    if (change < .05) goto L60;
	    if (URMath::Abs(grdold - grdnew) < dgmin) goto L60;
	    if (d < dmin_) {
		mnwarn("D", "MNHES1", "Step size too small for 1st drv.");
		goto L60;
	    }
	    d *= .2;
	}
//*-*-                                      loop satisfied = too many iter
        { std::ostringstream warning;
        warning << "Too many iterations on D1. " << grdold << " " << grdnew;
        mnwarn("D", "MNHES1", warning.str().c_str()); }
L60:
	fDgrd[i-1] = URMath::Max(dgmin,URMath::Abs(grdold - grdnew));
    }
//*-*-                                       end of first deriv. loop
    mninex(fX);
} /* mnhes1_ */

//______________________________________________________________________________
void URMinuit::mnimpr()
{
//*-*-*-*-*-*-*Attempts to improve on a good local minimum*-*-*-*-*-*-*-*-*-*
//*-*          ===========================================
//*-*        Attempts to improve on a good local minimum by finding a
//*-*        better one.   The quadratic part of FCN is removed by MNCALF
//*-*        and this transformed function is minimized using the simplex
//*-*        method from several random starting points.
//*-*        ref. -- Goldstein and Price, Math.Comp. 25, 569 (1971)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    static Double_urt rnum = 0;

    /* Local variables */
    Double_urt amax, ycalf, ystar, ystst;
    Double_urt pb, ep, wg, xi, sigsav, reg, sig2;
    Int_urt npfn, ndex, loop=0, i, j, ifail, iseed;
    Int_urt jhold, nloop, nparx, nparp1, jh, jl, iswtr;

    if (fNpar <= 0) return;
    if (fAmin == fUndefi) mnamin();
    fCstatu = "UNCHANGED ";
    fItaur  = 1;
    fEpsi   = fUp*.1;
    npfn    = fNfcn;
    nloop   = Int_urt(fWord7[1]);
    if (nloop <= 0) nloop = fNpar + 4;
    nparx  = fNpar;
    nparp1 = fNpar + 1;
    wg = 1 / Double_urt(fNpar);
    sigsav = fEDM;
    fApsi  = fAmin;
    iswtr   = fISW[4] - 2*fItaur;
    for (i = 1; i <= fNpar; ++i) {
	fXt[i-1]       = fX[i-1];
	fIMPRdsav[i-1] = fWerr[i-1];
	for (j = 1; j <= i; ++j) {
	    ndex = i*(i-1) / 2 + j;
	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1];
	    fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];
	}
    }
    mnvert(fP, fMaxint, fMaxint, fNpar, ifail);
    if (ifail >= 1) goto L280;
//*-*-              Save inverted matrix in VT
    for (i = 1; i <= fNpar; ++i) {
	ndex = i*(i-1) / 2;
	for (j = 1; j <= i; ++j) {
	    ++ndex;
	    fVthmat[ndex-1] = fP[i + j*fMaxpar - fMaxpar-1];
	}
    }
    loop = 0;

L20:
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = fIMPRdsav[i-1]*2;
	mnrn15(rnum, iseed);
	fX[i-1] = fXt[i-1] + fDirin[i-1]*2*(rnum - .5);
    }
    ++loop;
    reg = 2;
    if (fISW[4] >= 0) {
	std::printf("START ATTEMPT NO.%2d TO FIND NEW MINIMUM\n",loop);
    }
L30:
    mncalf(fX, ycalf);
    fAmin = ycalf;
//*-*-                                       . . . . set up  random simplex
    jl = nparp1;
    jh = nparp1;
    fIMPRy[nparp1-1] = fAmin;
    amax = fAmin;
    for (i = 1; i <= fNpar; ++i) {
	xi = fX[i-1];
	mnrn15(rnum, iseed);
	fX[i-1] = xi - fDirin[i-1]*(rnum - .5);
	mncalf(fX, ycalf);
	fIMPRy[i-1] = ycalf;
	if (fIMPRy[i-1] < fAmin) {
	    fAmin = fIMPRy[i-1];
	    jl    = i;
	} else if (fIMPRy[i-1] > amax) {
	    amax = fIMPRy[i-1];
	    jh   = i;
	}
	for (j = 1; j <= fNpar; ++j) { fP[j + i*fMaxpar - fMaxpar-1] = fX[j-1]; }
	fP[i + nparp1*fMaxpar - fMaxpar-1] = xi;
	fX[i-1] = xi;
    }

    fEDM = fAmin;
    sig2 = fEDM;
//*-*-                                       . . . . . . .  start main loop
L50:
    if (fAmin < 0)   goto L95;
    if (fISW[1] <= 2) goto L280;
    ep = fAmin*.1;
    if (sig2 < ep && fEDM < ep) goto L100;
    sig2 = fEDM;
    if (fNfcn - npfn > fNfcnmx) goto L300;
//*-*-        calculate new point * by reflection
    for (i = 1; i <= fNpar; ++i) {
	pb = 0;
	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }
	fPbar[i-1]  = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];
	fPstar[i-1] = fPbar[i-1]*2 - fP[i + jh*fMaxpar - fMaxpar-1]*1;
    }
    mncalf(fPstar, ycalf);
    ystar = ycalf;
    if (ystar >= fAmin) goto L70;
//*-*-        point * better than jl, calculate new point **
    for (i = 1; i <= fNpar; ++i) {
	fPstst[i-1] = fPstar[i-1]*2 + fPbar[i- 1]*-1;
    }
    mncalf(fPstst, ycalf);
    ystst = ycalf;
    if (ystst < fIMPRy[jl-1]) goto L67;
    mnrazz(ystar, fPstar, fIMPRy, jh, jl);
    goto L50;
L67:
    mnrazz(ystst, fPstst, fIMPRy, jh, jl);
    goto L50;
//*-*-        point * is not as good as jl
L70:
    if (ystar >= fIMPRy[jh-1]) goto L73;
    jhold = jh;
    mnrazz(ystar, fPstar, fIMPRy, jh, jl);
    if (jhold != jh) goto L50;
//*-*-        calculate new point **
L73:
    for (i = 1; i <= fNpar; ++i) {
	fPstst[i-1] = fP[i + jh*fMaxpar - fMaxpar-1]*.5 + fPbar[i-1]*.5;
    }
    mncalf(fPstst, ycalf);
    ystst = ycalf;
    if (ystst > fIMPRy[jh-1]) goto L30;
//*-*-    point ** is better than jh
    if (ystst < fAmin) goto L67;
    mnrazz(ystst, fPstst, fIMPRy, jh, jl);
    goto L50;
//*-*-                                       . . . . . .  end main loop
L95:
    if (fISW[4] >= 0) {
	std::printf(" AN IMPROVEMENT ON THE PREVIOUS MINIMUM HAS BEEN FOUND\n");
    }
    reg = .1;
//*-*-                                       . . . . . ask if point is new
L100:
    mninex(fX);
//    Eval(nparx, fGin, fAmin, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, fAmin, m_userParameterValue, 4);    ++fNfcn;
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = reg*fIMPRdsav[i-1];
	if (URMath::Abs(fX[i-1] - fXt[i-1]) > fDirin[i-1])     goto L150;
    }
    goto L230;
L150:
    fNfcnmx = fNfcnmx + npfn - fNfcn;
    npfn    = fNfcn;
    mnsimp();
    if (fAmin >= fApsi) goto L325;
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = fIMPRdsav[i-1]*.1;
	if (URMath::Abs(fX[i-1] - fXt[i-1]) > fDirin[i-1])     goto L250;
    }
L230:
    if (fAmin < fApsi) 	goto L350;
    goto L325;
/*                                        . . . . . . truly new minimum */
L250:
    fLnewmn = kurTRUE;
    if (fISW[1] >= 1) {
	fISW[1] = 1;
	fDcovar = URMath::Max(fDcovar,.5);
    } else fDcovar = 1;
    fItaur  = 0;
    fNfcnmx = fNfcnmx + npfn - fNfcn;
    fCstatu = "NEW MINIMU";
    if (fISW[4] >= 0) {
	std::printf(" IMPROVE HAS FOUND A TRULY NEW MINIMUM\n");
        std::printf(" *************************************\n");
    }
    return;
//*-*-                                       . . . return to previous region
L280:
    if (fISW[4] > 0) {
	std::printf(" COVARIANCE MATRIX WAS NOT POSITIVE-DEFINITE\n");
    }
    goto L325;
L300:
    fISW[0] = 1;
L325:
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = fIMPRdsav[i-1]*.01;
	fX[i-1]     = fXt[i-1];
    }
    fAmin = fApsi;
    fEDM  = sigsav;
L350:
    mninex(fX);
    if (fISW[4] > 0) {
	std::printf(" IMPROVE HAS RETURNED TO REGION OF ORIGINAL MINIMUM\n");
    }
    fCstatu = "UNCHANGED ";
    mnrset(0);
    if (fISW[1] < 2) goto L380;
    if (loop < nloop && fISW[0] < 1) goto L20;
L380:
    if (iswtr >= 0) mnprin(5, fAmin);
    fItaur = 0;
} /* mnimpr_ */

//______________________________________________________________________________
void URMinuit::mninex(Double_urt *pint)
{
//*-*-*-*-*Transforms from internal coordinates (PINT) to external (U)*-*-*-*
//*-*      ===========================================================
//*-*        The minimizing routines which work in
//*-*        internal coordinates call this routine before calling FCN.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt i, j;

    for (j = 1; j <= fNpar; ++j) {
	i = fNexofi[j-1];
//       if (fNvarl[i-1] == 1) {
	if (m_userParameterFlag[i] == 1) {
//           fU[i-1] = pint[j-1];
	    m_userParameterValue[i] = pint[j-1];
	} else {
//           fU[i-1] = fAlim[i-1] + (URMath::Sin(pint[j-1]) + 1)*.5*(fBlim[i-1] - fAlim[i-1]);
	    m_userParameterValue[i] = fAlim[i-1] + (URMath::Sin(pint[j-1]) + 1)*.5*(fBlim[i-1] - fAlim[i-1]);
	}
    }
} /* mninex_ */

//______________________________________________________________________________
void URMinuit::mninit(Int_urt i1, Int_urt i2, Int_urt i3)
{
//*-*-*-*-*-*Main initialization member function for MINUIT*-*-*-*-*-*-*-*-*
//*-*        ==============================================
//*-*     It initializes some constants
//*-*                (including the logical I/O unit nos.),
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt piby2, epsp1, epsbak, epstry, distnn;
    Int_urt i, idb;

//*-*-           I/O unit numbers
    fIsysrd = i1;
    fIsyswr = i2;
    fIstkwr[0] = fIsyswr;
    fNstkwr = 1;
    fIsyssa = i3;
    fNstkrd = 0;
//*-*-              version identifier
    fCvrsn  = "95.03++ ";
//*-*-              some CONSTANT
    fMaxint = fMaxpar;
    fMaxext = 2*fMaxpar;
    fUndefi = -54321;
    fBigedm = 123456;
    fCundef = kUndefinedParameterName; // ")UNDEFINED";
    fCovmes[0] = "NO ERROR MATRIX       ";
    fCovmes[1] = "ERR MATRIX APPROXIMATE";
    fCovmes[2] = "ERR MATRIX NOT POS-DEF";
    fCovmes[3] = "ERROR MATRIX ACCURATE ";
//*-*-               some starting values
    fNblock = 0;
    fIcomnd = 0;
    fCtitl  = fCundef;
    fCfrom  = "INPUT   ";
    fNfcn   = 0;
    fNfcnfr = fNfcn;
    fCstatu = "INITIALIZE";
    fISW[2] = 0;
    fISW[3] = 0;
    fISW[4] = 1;
//*-*-        ISW(6)=0 for batch jobs,  =1 for interactive jobs
//*-*-                     =-1 for originally interactive temporarily batch

    fISW[5] = 0;
//    if (intrac(&dummy)) fISW[5] = 1;
//*-*-       DEBUG options set to default values
    for (idb = 0; idb <= 10; ++idb) { fIdbg[idb] = 0; }
    fLrepor = kurFALSE;
    fLwarn  = kurTRUE;
    fLimset = kurFALSE;
    fLnewmn = kurFALSE;
    fIstrat = 1;
    fItaur  = 0;
//*-*-       default page dimensions and 'new page' carriage control integer
    fNpagwd = 120;
    fNpagln = 56;
    fNewpag = 1;
    if (fISW[5] > 0) {
	fNpagwd = 80;
	fNpagln = 30;
	fNewpag = 0;
    }
    fUp = 1;
    fUpdflt = fUp;
//*-*-                  determine machine accuracy epsmac
    epstry = .5;
    for (i = 1; i <= 100; ++i) {
	epstry *= .5;
	epsp1 = epstry + 1;
	mntiny(epsp1, epsbak);
	if (epsbak < epstry) goto L35;
    }
    epstry = 1e-7;
    fEpsmac = epstry*4;
    std::printf(" MNINIT UNABLE TO DETERMINE ARITHMETIC PRECISION. WILL ASSUME:%g\n",fEpsmac);
L35:
    fEpsmac = epstry*8;
    fEpsma2 = URMath::Sqrt(fEpsmac)*2;
//*-*-                the vlims are a non-negligible distance from pi/2
//*-*-        used by MNPINT to set variables "near" the physical limits
    piby2   = URMath::ATan(1)*2;
    distnn  = URMath::Sqrt(fEpsma2)*8;
    fVlimhi =  piby2 - distnn;
    fVlimlo = -piby2 + distnn;
  //    *m_logStream << "size of m_userParameterName: " << m_userParameterName.size() << endl;
    mncler();
//    std::printf("  MINUIT RELEASE %s INITIALIZED.   DIMENSIONS 100/50  EPSMAC=%g",(const char*)fCvrsn,fEpsmac);
} /* mninit_ */

//______________________________________________________________________________
void URMinuit::mnlims()
{
//*-*-*-*Interprets the SET LIM command, to reset the parameter limits*-*-*-*
//*-*    =============================================================
//*-*       Called from MNSET
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dxdi, snew;
    Int_urt kint, i2, newcod, ifx=0, inu;

    fCfrom  = "SET LIM ";
    fNfcnfr = fNfcn;
    fCstatu = "NO CHANGE ";
    i2 = Int_urt(fWord7[0]);
    if (i2 > fMaxext || i2 < 0) goto L900;
    if (i2 > 0) goto L30;
//*-*-                                    set limits on all parameters
    newcod = 4;
    if (fWord7[1] == fWord7[2]) newcod = 1;
    for (inu = 1; inu <= fNu; ++inu) {
//       if (fNvarl[inu-1] <= 0) continue;
	if (m_userParameterFlag[inu] <= 0) continue;
//       if (fNvarl[inu-1] == 1 && newcod == 1) continue;
	if (m_userParameterFlag[inu] == 1 && newcod == 1) continue;
//        kint = fNiofex[inu-1];
	kint = m_userParameterIdToInternalId[inu];
//*-*-            see if parameter has been fixed
	if (kint <= 0) {
	    if (fISW[4] >= 0) {
		std::printf("           LIMITS NOT CHANGED FOR FIXED PARAMETER:%4d\n",inu);
	    }
	    continue;
	}
	if (newcod == 1) {
//*-*-           remove limits from parameter
	    if (fISW[4] > 0) {
                std::printf(" LIMITS REMOVED FROM PARAMETER  :%3d\n",inu);
	    }
	    fCstatu = "NEW LIMITS";
	    mndxdi(fX[kint-1], kint-1, dxdi);
	    snew           = fGstep[kint-1]*dxdi;
	    fGstep[kint-1] = URMath::Abs(snew);
//	    fNvarl[inu-1]  = 1;
            m_userParameterFlag[inu]  = 1;
	} else {
//*-*-            put limits on parameter
	    fAlim[inu-1] = URMath::Min(fWord7[1],fWord7[2]);
	    fBlim[inu-1] = URMath::Max(fWord7[1],fWord7[2]);
	    if (fISW[4] > 0) {
		std::printf(" PARAMETER %3d LIMITS SET TO  %15.5g%15.5g\n",inu,fAlim[inu-1],fBlim[inu-1]);
	    }
//            fNvarl[inu-1]  = 4;
	    m_userParameterFlag[inu]  = 4;
	    fCstatu        = "NEW LIMITS";
	    fGstep[kint-1] = -.1;
	}
    }
    goto L900;
//*-*-                                      set limits on one parameter
L30:
    //if (fNvarl[i2-1] <= 0) {
    if (m_userParameterFlag[i2] <= 0) {
        std::printf(" PARAMETER %3d IS NOT VARIABLE.\n", i2);
	goto L900;
    }
//    kint = fNiofex[i2-1];
    kint = m_userParameterIdToInternalId[i2];
//*-*-                                      see if parameter was fixed
    if (kint == 0) {
        std::printf(" REQUEST TO CHANGE LIMITS ON FIXED PARAMETER:%3d\n",i2);
	for (ifx = 1; ifx <= fNpfix; ++ifx) {
	    if (i2 == fIpfix[ifx-1]) goto L92;
	}
	std::printf(" MINUIT BUG IN MNLIMS. SEE F. JAMES\n");
L92:
	;
    }
    if (fWord7[1] != fWord7[2]) goto L235;
//*-*-                                      remove limits
//    if (fNvarl[i2-1] != 1) {
    if (m_userParameterFlag[i2] != 1) {
	if (fISW[4] > 0) {
            std::printf(" LIMITS REMOVED FROM PARAMETER  %2d\n",i2);
	}
	fCstatu = "NEW LIMITS";
	if (kint <= 0) {
	    fGsteps[ifx-1] = URMath::Abs(fGsteps[ifx-1]);
	} else {
	    mndxdi(fX[kint-1], kint-1, dxdi);
	    if (URMath::Abs(dxdi) < .01) dxdi = .01;
	    fGstep[kint-1] = URMath::Abs(fGstep[kint-1]*dxdi);
	    fGrd[kint-1]  *= dxdi;
	}
//        fNvarl[i2-1] = 1;
	m_userParameterFlag[i2] = 1;
    } else {
	std::printf(" NO LIMITS SPECIFIED.  PARAMETER %3d IS ALREADY UNLIMITED.  NO CHANGE.\n",i2);
    }
    goto L900;
//*-*-                                       put on limits
L235:
    fAlim[i2-1]  = URMath::Min(fWord7[1],fWord7[2]);
    fBlim[i2-1]  = URMath::Max(fWord7[1],fWord7[2]);
//    fNvarl[i2-1] = 4;
    m_userParameterFlag[i2] = 4;
    if (fISW[4] > 0) {
	std::printf(" PARAMETER %3d LIMITS SET TO  %15.5g%15.5g\n",i2,fAlim[i2-1],fBlim[i2-1]);
    }
    fCstatu = "NEW LIMITS";
    if (kint <= 0) fGsteps[ifx-1] = -.1;
    else           fGstep[kint-1] = -.1;

L900:
    if (fCstatu != "NO CHANGE ") {
	mnexin(fX);
	mnrset(1);
    }
} /* mnlims_ */

//______________________________________________________________________________
void URMinuit::mnline(Double_urt *start, Double_urt fstart, Double_urt *step, Double_urt slope, Double_urt toler)
{
//*-*-*-*-*-*-*-*-*-*Perform a line search from position START*-*-*-*-*-*-*-*
//*-*                =========================================
//*-*        along direction STEP, where the length of vector STEP
//*-*                   gives the expected position of minimum.
//*-*        FSTART is value of function at START
//*-*        SLOPE (if non-zero) is df/dx along STEP at START
//*-*        TOLER is initial tolerance of minimum in direction STEP
//*-*
//*-* SLAMBG and ALPHA control the maximum individual steps allowed.
//*-* The first step is always =1. The max length of second step is SLAMBG.
//*-* The max size of subsequent steps is the maximum previous successful
//*-*   step multiplied by ALPHA + the size of most recent successful step,
//*-*   but cannot be smaller than SLAMBG.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt xpq[12], ypq[12], slam, sdev, coeff[3], denom, flast;
    Double_urt fvals[3], xvals[3], f1, fvmin, xvmin, ratio, f2, f3, fvmax;
    Double_urt toler8, toler9, overal, undral, slamin, slamax, slopem;
    Int_urt i, nparx=0, nvmax=0, nxypt, kk, ipt;
    Bool_urt ldebug;
    string cmess;
    char chpq[13];
    int     l65, l70, l80;

    /* Function Body */
    l65 = 0; l70 = 0; l80 = 0;
    ldebug = fIdbg[1] >= 1;
//*-*-                 starting values for overall limits on total step SLAM
    overal = 1e3;
    undral = -100;
//*-*-                             debug check if start is ok
    if (ldebug) {
	mninex(&start[0]);
//	Eval(nparx, fGin, f1, fU, 4);	++fNfcn;
        Eval(nparx, fGin, f1, m_userParameterValue, 4);	++fNfcn;
	if (f1 != fstart) {
	    std::printf(" MNLINE start point not consistent, F values, parameters=\n");
	    for (kk = 1; kk <= fNpar; ++kk) {
                std::printf("  %14.5e\n",fX[kk-1]);
	    }
	}
    }
//*-*-                                     . set up linear search along STEP
    fvmin   = fstart;
    xvmin   = 0;
    nxypt   = 1;
    chpq[0] = charal[0];
    xpq[0]  = 0;
    ypq[0]  = fstart;
//*-*-              SLAMIN = smallest possible value of ABS(SLAM)
    slamin = 0;
    for (i = 1; i <= fNpar; ++i) {
	if (step[i-1] != 0) {
           ratio = URMath::Abs(start[i-1] / step[i-1]);
           if (slamin == 0)    slamin = ratio;
           if (ratio < slamin) slamin = ratio;
        }
	fX[i-1] = start[i-1] + step[i-1];
    }
    if (slamin == 0) slamin = fEpsmac;
    slamin *= fEpsma2;
    nparx = fNpar;

    mninex(fX);
//    Eval(nparx, fGin, f1, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, f1, m_userParameterValue, 4);    ++fNfcn;
    ++nxypt;
    chpq[nxypt-1] = charal[nxypt-1];
    xpq[nxypt-1] = 1;
    ypq[nxypt-1] = f1;
    if (f1 < fstart) {
	fvmin = f1;
	xvmin = 1;
    }
//*-*-                        . quadr interp using slope GDEL and two points
    slam   = 1;
    toler8 = toler;
    slamax = 5;
    flast  = f1;
//*-*-                        can iterate on two-points (cut) if no imprvmnt

    do {
       denom = (flast - fstart - slope*slam)*2 / (slam*slam);
       slam  = 1;
       if (denom != 0)    slam = -slope / denom;
       if (slam < 0)      slam = slamax;
       if (slam > slamax) slam = slamax;
       if (slam < toler8) slam = toler8;
       if (slam < slamin) {
          l80 = 1;
          break;
       }
       if (URMath::Abs(slam - 1) < toler8 && f1 < fstart) {
          l70 = 1;
          break;
       }
       if (URMath::Abs(slam - 1) < toler8) slam = toler8 + 1;
       if (nxypt >= 12) {
         l65 = 1;
         break;
       }
       for (i = 1; i <= fNpar; ++i) { fX[i-1] = start[i-1] + slam*step[i-1]; }
       mninex(fX);
       nparx = fNpar;
//       Eval(nparx, fGin, f2, fU, 4);    ++fNfcn;
       Eval(nparx, fGin, f2, m_userParameterValue, 4);    ++fNfcn;
       ++nxypt;
       chpq[nxypt-1] = charal[nxypt-1];
       xpq[nxypt-1]  = slam;
       ypq[nxypt-1]  = f2;
       if (f2 < fvmin) {
          fvmin = f2;
          xvmin = slam;
       }
       if (fstart == fvmin) {
          flast  = f2;
          toler8 = toler*slam;
          overal = slam - toler8;
          slamax = overal;
       }
    } while (fstart == fvmin);

    if (!l65 && !l70 && !l80) {
//*-*-                                       . quadr interp using 3 points
       xvals[0] = xpq[0];
       fvals[0] = ypq[0];
       xvals[1] = xpq[nxypt-2];
       fvals[1] = ypq[nxypt-2];
       xvals[2] = xpq[nxypt-1];
       fvals[2] = ypq[nxypt-1];
//*-*-                            begin iteration, calculate desired step
       do {
          slamax = URMath::Max(slamax,URMath::Abs(xvmin)*2);
          mnpfit(xvals, fvals, 3, coeff, sdev);
          if (coeff[2] <= 0) {
             slopem = coeff[2]*2*xvmin + coeff[1];
             if (slopem <= 0)  slam = xvmin + slamax;
             else              slam = xvmin - slamax;
          } else {
             slam = -coeff[1] / (coeff[2]*2);
             if (slam > xvmin + slamax) slam = xvmin + slamax;
             if (slam < xvmin - slamax) slam = xvmin - slamax;
          }
          if (slam > 0) { if (slam > overal) slam = overal; }
          else          { if (slam < undral) slam = undral; }

//*-*-              come here if step was cut below
          do {
             toler9 = URMath::Max(toler8,URMath::Abs(toler8*slam));
             for (ipt = 1; ipt <= 3; ++ipt) {
                if (URMath::Abs(slam - xvals[ipt-1]) < toler9) {
                   l70 = 1;
                   break;
                }
             }
             if (l70) break;
//*-*-               take the step
             if (nxypt >= 12) {
                l65 = 1;
                break;
             }
             for (i = 1; i <= fNpar; ++i) { fX[i-1] = start[i-1] + slam*step[i-1]; }
             mninex(fX);
//             Eval(nparx, fGin, f3, fU, 4);    ++fNfcn;
             Eval(nparx, fGin, f3, m_userParameterValue, 4);    ++fNfcn;
             ++nxypt;
             chpq[nxypt-1] = charal[nxypt-1];
             xpq[nxypt-1]  = slam;
             ypq[nxypt-1]  = f3;
//*-*-            find worst previous point out of three
             fvmax = fvals[0];
             nvmax = 1;
             if (fvals[1] > fvmax) {
                fvmax = fvals[1];
                nvmax = 2;
             }
             if (fvals[2] > fvmax) {
                fvmax = fvals[2];
                nvmax = 3;
             }
//*-*-             if latest point worse than all three previous, cut step
             if (f3 >= fvmax) {
                if (nxypt >= 12) {
                   l65 = 1;
                   break;
                }
                if (slam > xvmin) overal = URMath::Min(overal,slam - toler8);
                if (slam < xvmin) undral = URMath::Max(undral,slam + toler8);
                slam = (slam + xvmin)*.5;
             }
          } while (f3 >= fvmax);

//*-*-             prepare another iteration, replace worst previous point
          if (l65 || l70) break;

          xvals[nvmax-1] = slam;
          fvals[nvmax-1] = f3;
          if (f3 < fvmin) {
             fvmin = f3;
             xvmin = slam;
          } else {
             if (slam > xvmin) overal = URMath::Min(overal,slam - toler8);
             if (slam < xvmin) undral = URMath::Max(undral,slam + toler8);
          }
       } while (nxypt < 12);
    }

//*-*-                                           . . end of iteration . . .
//*-*-           stop because too many iterations
    if (!l70 && !l80) {
       cmess = " LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS ";
       if (ldebug) {
          std::printf(" MNLINE DEBUG: steps=\n");
          for (kk = 1; kk <= fNpar; ++kk) {
             std::printf("  %12.4g\n",step[kk-1]);
          }
       }
    }
//*-*-           stop because within tolerance
    if (l70) cmess = " LINE SEARCH HAS ATTAINED TOLERANCE ";
    if (l80) cmess = " STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM";

    fAmin = fvmin;
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = step[i-1]*xvmin;
	fX[i-1]     = start[i-1] + fDirin[i-1];
    }
    mninex(fX);
    if (xvmin < 0) {
	mnwarn("D", "MNLINE", " LINE MINIMUM IN BACKWARDS DIRECTION");
    }
    if (fvmin == fstart) {
	mnwarn("D", "MNLINE", " LINE SEARCH FINDS NO IMPROVEMENT ");
    }
    if (ldebug) {
        std::printf(" AFTER%3d POINTS,%s\n",nxypt,cmess.c_str());
	mnplot(xpq, ypq, chpq, nxypt, fNpagwd, fNpagln);
    }
} /* mnline_ */

//______________________________________________________________________________
void URMinuit::mnmatu(Int_urt kode)
{
//*-*-*-*-*-*-*-*Prints the covariance matrix v when KODE=1*-*-*-*-*-*-*-*-*
//*-*            ==========================================
//*-*        always prints the global correlations, and
//*-*        calculates and prints the individual correlation coefficients
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Int_urt ndex, i, j, m, n, ncoef, nparm, id, it, ix;
    Int_urt nsofar, ndi, ndj, iso, isw2, isw5;
    string ctemp;

    isw2 = fISW[1];
    if (isw2 < 1) {
	std::printf("%s\n",(fCovmes[isw2]).c_str());
	return;
    }
    if (fNpar == 0) {
	std::printf(" MNMATU: NPAR=0\n");
	return;
    }
//*-*-                                      . . . . .external error matrix
    if (kode == 1) {
	isw5    = fISW[4];
	fISW[4] = 2;
	mnemat(fP, fMaxint);
	if (isw2 < 3) {
	    std::printf("%s\n",(fCovmes[isw2]).c_str());
	}
	fISW[4] = isw5;
    }
//*-*-                                      . . . . . correlation coeffs. .
    if (fNpar <= 1) return;
    mnwerr();
//*-*-    NCOEF is number of coeff. that fit on one line, not to exceed 20
    ncoef = (fNpagwd - 19) / 6;
    ncoef = URMath::Min(ncoef,20);
    nparm = URMath::Min(fNpar,ncoef);
    std::printf(" PARAMETER  CORRELATION COEFFICIENTS  \n");
    ctemp = "       NO.  GLOBAL";
    for (id = 1; id <= nparm; ++id) {
       ostringstream localTemp;
       localTemp << setw(7) << right << fNexofi[id-1];
       ctemp += localTemp.str(); // Form(" %6d",);
    }
    std::printf("%s\n",ctemp.c_str());
    for (i = 1; i <= fNpar; ++i) {
	ix  = fNexofi[i-1];
	ndi = i*(i + 1) / 2;
	for (j = 1; j <= fNpar; ++j) {
	    m    = URMath::Max(i,j);
	    n    = URMath::Min(i,j);
	    ndex = m*(m-1) / 2 + n;
	    ndj  = j*(j + 1) / 2;
	    fMATUvline[j-1] = fVhmat[ndex-1] / URMath::Sqrt(URMath::Abs(fVhmat[ndi-1]*fVhmat[ndj-1]));
	}
	nparm = URMath::Min(fNpar,ncoef);
        ostringstream octemp;
        octemp << setw(7) << right << ix << ": " << fixed << setw(9) << setprecision(3) << fGlobcc[i-1];
        ctemp = octemp.str();
	for (it = 1; it <= nparm; ++it) {
           ostringstream localTemp;
           localTemp << fixed << setw(7) << setprecision(3) << fMATUvline[it-1];
           ctemp += localTemp.str(); //Form(" %6.3f",);
	}
        std::printf("%s\n",ctemp.c_str());
	if (i <= nparm) continue;
	for (iso = 1; iso <= 10; ++iso) {
            ctemp = "                  ";
	    nsofar = nparm;
	    nparm  = URMath::Min(fNpar,nsofar + ncoef);
	    for (it = nsofar + 1; it <= nparm; ++it) {
               ostringstream localTemp;
               localTemp << fixed << setw(7) << setprecision(3) << fMATUvline[it-1];
               ctemp = ctemp + localTemp.str(); // Form(" %6.3f",fMATUvline[it-1]);
	    }
            std::printf("%s\n",ctemp.c_str());
	    if (i <= nparm) break;
	}
    }
    if (isw2 < 3) {
        std::printf(" %s\n",(fCovmes[isw2]).c_str());
    }
} /* mnmatu_ */

//______________________________________________________________________________
void URMinuit::mnmigr()
{
//*-*-*-*-*-*-*-*-*Performs a local function minimization*-*-*-*-*-*-*-*-*-*
//*-*              ======================================
//*-*        Performs a local function minimization using basically the
//*-*        method of Davidon-Fletcher-Powell as modified by Fletcher
//*-*        ref. -- Fletcher, Comp.J. 13,317 (1970)   "switching method"
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt gdel, gami, vlen, dsum, gssq, vsum, d;
    Double_urt fzero, fs, ri, delgam, rhotol;
    Double_urt gdgssq, gvg, vgi;
    Int_urt npfn, ndex, iext, i, j, m, n, npsdf, nparx;
    Int_urt iswtr, lined2, kk, nfcnmg, nrstrt,iter;
    Bool_urt ldebug;
    Double_urt toler = 0.05;

    if (fNpar <= 0) return;
    if (fAmin == fUndefi) mnamin();
    ldebug  = kurFALSE; if ( fIdbg[4] >= 1) ldebug = kurTRUE;
    fCfrom  = "MIGRAD  ";
    fNfcnfr = fNfcn;
    nfcnmg  = fNfcn;
    fCstatu = "INITIATE  ";
    iswtr   = fISW[4] - 2*fItaur;
    npfn    = fNfcn;
    nparx   = fNpar;
    vlen    = (Double_urt) (fNpar*(fNpar + 1) / 2);
    nrstrt  = 0;
    npsdf   = 0;
    lined2  = 0;
    fISW[3] = -1;
    rhotol  = fApsi*.001;
    if (iswtr >= 1) {
	std::printf(" START MIGRAD MINIMIZATION.  STRATEGY%2d.  CONVERGENCE WHEN EDM .LT.%9.2e\n",fIstrat,rhotol);
    }
//*-*-                                          initialization strategy
    if (fIstrat < 2 || fISW[1] >= 3) goto L2;
//*-*-                               come (back) here to restart completely
L1:
    if (nrstrt > fIstrat) {
	fCstatu = "FAILED    ";
	fISW[3] = -1;
	goto L230;
    }
//*-*-                                     . get full covariance and gradient
    mnhess();
    mnwerr();
    npsdf = 0;
    if (fISW[1] >= 1) goto L10;
//*-*-                                       . get gradient at start point
L2:
    mninex(fX);
    if (fISW[2] == 1) {
//	Eval(nparx, fGin, fzero, fU, 2);	++fNfcn;
        Eval(nparx, fGin, fzero, m_userParameterValue, 2);	++fNfcn;
    }
    mnderi();
    if (fISW[1] >= 1) goto L10;
//*-*-                                  sometimes start with diagonal matrix
    for (i = 1; i <= fNpar; ++i) {
	fMIGRxxs[i-1]  = fX[i-1];
	fMIGRstep[i-1] = 0;
    }
//*-*-                          do line search if second derivative negative
    ++lined2;
    if (lined2 < (fIstrat + 1)*fNpar) {
	for (i = 1; i <= fNpar; ++i) {
	    if (fG2[i-1] > 0) continue;
	    if (fGrd[i-1] > 0) fMIGRstep[i-1] = -URMath::Abs(fGstep[i-1]);
	    else               fMIGRstep[i-1] =  URMath::Abs(fGstep[i-1]);
	    gdel = fMIGRstep[i-1]*fGrd[i-1];
	    fs   = fAmin;
	    mnline(fMIGRxxs, fs, fMIGRstep, gdel, toler);
	    mnwarn("D", "MNMIGR", "Negative G2 line search");
	    iext = fNexofi[i-1];
	    if (ldebug) {
		std::printf(" Negative G2 line search, param %3d %13.3g%13.3g\n",iext,fs,fAmin);
	    }
	    goto L2;
	}
    }
//*-*-                          make diagonal error matrix
    for (i = 1; i <= fNpar; ++i) {
	ndex = i*(i-1) / 2;
	for (j = 1; j <= i-1; ++j) {
	    ++ndex;
	    fVhmat[ndex-1] = 0;
	}
	++ndex;
	if (fG2[i-1] <= 0) fG2[i-1] = 1;
	fVhmat[ndex-1] = 2 / fG2[i-1];
    }
    fDcovar = 1;
    if (ldebug) {
	std::printf(" DEBUG MNMIGR, STARTING MATRIX DIAGONAL,  VHMAT=\n");
	for (kk = 1; kk <= Int_urt(vlen); ++kk) {
	    std::printf(" %10.2g\n",fVhmat[kk-1]);
	}
    }
//*-*-                                        ready to start first iteration
L10:
    ++nrstrt;
    if (nrstrt > fIstrat + 1) {
	fCstatu = "FAILED    ";
	goto L230;
    }
    fs = fAmin;
//*-*-                                       . . . get EDM and set up loop
    fEDM = 0;
    for (i = 1; i <= fNpar; ++i) {
	fMIGRgs[i-1]  = fGrd[i-1];
	fMIGRxxs[i-1] = fX[i-1];
	ndex     = i*(i-1) / 2;
	for (j = 1; j <= i-1; ++j) {
	    ++ndex;
	    fEDM += fMIGRgs[i-1]*fVhmat[ndex-1]*fMIGRgs[j-1];
	}
	++ndex;
	fEDM += fMIGRgs[i-1]*fMIGRgs[i-1]*.5*fVhmat[ndex-1];
    }
    fEDM = fEDM*.5*(fDcovar*3 + 1);
    if (fEDM < 0) {
	mnwarn("W", "MIGRAD", "STARTING MATRIX NOT POS-DEFINITE.");
	fISW[1] = 0;
	fDcovar = 1;
	goto L2;
    }
    if (fISW[1] == 0) fEDM = fBigedm;
    iter = 0;
    mninex(fX);
    mnwerr();
    if (iswtr >= 1) mnprin(3, fAmin);
    if (iswtr >= 2) mnmatu(0);
//*-*-                                       . . . . .  start main loop
L24:
    if (fNfcn - npfn >= fNfcnmx) goto L190;
    gdel = 0;
    gssq = 0;
    for (i = 1; i <= fNpar; ++i) {
	ri = 0;
	gssq += fMIGRgs[i-1]*fMIGRgs[i-1];
	for (j = 1; j <= fNpar; ++j) {
	    m    = URMath::Max(i,j);
	    n    = URMath::Min(i,j);
	    ndex = m*(m-1) / 2 + n;
	    ri  += fVhmat[ndex-1]*fMIGRgs[j-1];
	}
	fMIGRstep[i-1] = ri*-.5;
	gdel += fMIGRstep[i-1]*fMIGRgs[i-1];
    }
    if (gssq == 0) {
	mnwarn("D", "MIGRAD", " FIRST DERIVATIVES OF FCN ARE ALL ZERO");
	goto L300;
    }
//*-*-                if gdel positive, V not posdef
    if (gdel >= 0) {
	mnwarn("D", "MIGRAD", " NEWTON STEP NOT DESCENT.");
	if (npsdf == 1) goto L1;
	mnpsdf();
	npsdf = 1;
	goto L24;
    }
//*-*-                                       . . . . do line search
    mnline(fMIGRxxs, fs, fMIGRstep, gdel, toler);
    if (fAmin == fs) goto L200;
    fCfrom  = "MIGRAD  ";
    fNfcnfr = nfcnmg;
    fCstatu = "PROGRESS  ";
//*-*-                                       . get gradient at new point
    mninex(fX);
    if (fISW[2] == 1) {
//	Eval(nparx, fGin, fzero, fU, 2);	++fNfcn;
        Eval(nparx, fGin, fzero, m_userParameterValue, 2);	++fNfcn;
    }
    mnderi();
//*-*-                                        . calculate new EDM
    npsdf = 0;
L81:
    fEDM = 0;
    gvg = 0;
    delgam = 0;
    gdgssq = 0;
    for (i = 1; i <= fNpar; ++i) {
	ri  = 0;
	vgi = 0;
	for (j = 1; j <= fNpar; ++j) {
	    m    = URMath::Max(i,j);
	    n    = URMath::Min(i,j);
	    ndex = m*(m-1) / 2 + n;
	    vgi += fVhmat[ndex-1]*(fGrd[j-1] - fMIGRgs[j-1]);
	    ri  += fVhmat[ndex-1]*fGrd[j-1];
	}
	fMIGRvg[i-1] = vgi*.5;
	gami    = fGrd[i-1] - fMIGRgs[i-1];
	gdgssq += gami*gami;
	gvg    += gami*fMIGRvg[i-1];
	delgam += fDirin[i-1]*gami;
	fEDM   += fGrd[i-1]*ri*.5;
    }
    fEDM = fEDM*.5*(fDcovar*3 + 1);
//*-*-                         . if EDM negative,  not positive-definite
    if (fEDM < 0 || gvg <= 0) {
	mnwarn("D", "MIGRAD", "NOT POS-DEF. EDM OR GVG NEGATIVE.");
	fCstatu = "NOT POSDEF";
	if (npsdf == 1) goto L230;
	mnpsdf();
	npsdf = 1;
	goto L81;
    }
//*-*-                           print information about this iteration
    ++iter;
    if (iswtr >= 3 || ( iswtr == 2 && iter % 10 == 1 ) ) {
	mnwerr();
	mnprin(3, fAmin);
    }
    if (gdgssq == 0) {
	mnwarn("D", "MIGRAD", "NO CHANGE IN FIRST DERIVATIVES OVER LAST STEP");
    }
    if (delgam < 0) {
	mnwarn("D", "MIGRAD", "FIRST DERIVATIVES INCREASING ALONG SEARCH LINE");
    }
//*-*-                                       .  update covariance matrix
    fCstatu = "IMPROVEMNT";
    if (ldebug) {
	std::printf(" VHMAT 1 =\n");
	for (kk = 1; kk <= 10; ++kk) {
            std::printf(" %10.2g\n",fVhmat[kk-1]);
	}
    }
    dsum = 0;
    vsum = 0;
    for (i = 1; i <= fNpar; ++i) {
	for (j = 1; j <= i; ++j) {
	    if(delgam == 0 || gvg == 0) d = 0;
            else d = fDirin[i-1]*fDirin[j-1] / delgam - fMIGRvg[i-1]*fMIGRvg[j-1] / gvg;
	    dsum += URMath::Abs(d);
	    ndex  = i*(i-1) / 2 + j;
	    fVhmat[ndex-1] += d*2;
	    vsum += URMath::Abs(fVhmat[ndex-1]);
	}
    }
//*-*-               smooth local fluctuations by averaging DCOVAR
    fDcovar = (fDcovar + dsum / vsum)*.5;
    if (iswtr >= 3 || ldebug) {
	std::printf(" RELATIVE CHANGE IN COV. MATRIX=%5.1f per cent\n",fDcovar*100);
    }
    if (ldebug) {
	std::printf(" VHMAT 2 =\n");
	for (kk = 1; kk <= 10; ++kk) {
            std::printf(" %10.3g\n",fVhmat[kk-1]);
	}
    }
    if (delgam <= gvg) goto L135;
    for (i = 1; i <= fNpar; ++i) {
	fMIGRflnu[i-1] = fDirin[i-1] / delgam - fMIGRvg[i-1] / gvg;
    }
    for (i = 1; i <= fNpar; ++i) {
	for (j = 1; j <= i; ++j) {
	    ndex = i*(i-1) / 2 + j;
	    fVhmat[ndex-1] += gvg*2*fMIGRflnu[i-1]*fMIGRflnu[j-1];
	}
    }
L135:
//*-*-                                             and see if converged
    if (fEDM < rhotol*.1) goto L300;
//*-*-                                   if not, prepare next iteration
    for (i = 1; i <= fNpar; ++i) {
	fMIGRxxs[i-1] = fX[i-1];
	fMIGRgs[i-1]  = fGrd[i-1];
    }
    fs = fAmin;
    if (fISW[1] == 0 && fDcovar < .5)  fISW[1] = 1;
    if (fISW[1] == 3 && fDcovar > .1)  fISW[1] = 1;
    if (fISW[1] == 1 && fDcovar < .05) fISW[1] = 3;
    goto L24;
//*-*-                                       . . . . .  end main loop
//*-*-                                        . . call limit in MNMIGR
L190:
    fISW[0] = 1;
    if (fISW[4] >= 0) {
	std::printf(" CALL LIMIT EXCEEDED IN MIGRAD.\n");
    }
    fCstatu = "CALL LIMIT";
    goto L230;
//*-*-                                        . . fails to improve . .
L200:
    if (iswtr >= 1) {
	std::printf(" MIGRAD FAILS TO FIND IMPROVEMENT\n");
    }
    for (i = 1; i <= fNpar; ++i) { fX[i-1] = fMIGRxxs[i-1]; }
    if (fEDM < rhotol) goto L300;
    if (fEDM < URMath::Abs(fEpsma2*fAmin)) {
	if (iswtr >= 0) {
	    std::printf(" MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT.\n");
	}
	goto L300;
    }
    if (fIstrat < 1) {
	if (fISW[4] >= 0) {
	    std::printf(" MIGRAD FAILS WITH STRATEGY=0.   WILL TRY WITH STRATEGY=1.\n");
	}
	fIstrat = 1;
    }
    goto L1;
//*-*-                                        . . fails to converge
L230:
    if (iswtr >= 0) {
	std::printf(" MIGRAD TERMINATED WITHOUT CONVERGENCE.\n");
    }
    if (fISW[1] == 3) fISW[1] = 1;
    fISW[3] = -1;
    goto L400;
//*-*-                                        . . apparent convergence
L300:
    if (iswtr >= 0) {
	std::printf(" MIGRAD MINIMIZATION HAS CONVERGED.\n");
    }
    if (fItaur == 0) {
	if (fIstrat >= 2 || (fIstrat == 1 && fISW[1] < 3)) {
	    if (fISW[4] >= 0) {
		std::printf(" MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.\n");
	    }
	    mnhess();
	    mnwerr();
	    npsdf = 0;
	    if (fEDM > rhotol) goto L10;
	}
    }
    fCstatu = "CONVERGED ";
    fISW[3] = 1;
//*-*-                                          come here in any case
L400:
    fCfrom  = "MIGRAD  ";
    fNfcnfr = nfcnmg;
    mninex(fX);
    mnwerr();
    if (iswtr >= 0) mnprin(3, fAmin);
    if (iswtr >= 1) mnmatu(1);
} /* mnmigr_ */

//______________________________________________________________________________
void URMinuit::mnmnos()
{
//*-*-*-*-*-*-*-*-*-*-*Performs a MINOS error analysis*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ===============================
//*-*        Performs a MINOS error analysis on those parameters for
//*-*        which it is requested on the MINOS command by calling
//*-*        MNMNOT for each parameter requested.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt val2mi, val2pl;
    Int_urt nbad, ilax, ilax2, ngood, nfcnmi, iin, knt;

    if (fNpar <= 0) goto L700;
    ngood = 0;
    nbad = 0;
    nfcnmi = fNfcn;
//*-*-                                     . loop over parameters requested
    for (knt = 1; knt <= fNpar; ++knt) {
	if (Int_urt(fWord7[1]) == 0) {
	    ilax = fNexofi[knt-1];
	} else {
	    if (knt >= 7) break;
	    ilax = Int_urt(fWord7[knt]);
	    if (ilax == 0) break;
	    if (ilax > 0 && ilax <= fNu) {
//               if (fNiofex[ilax-1] > 0) goto L565;
               if (m_userParameterIdToInternalId[ilax] > 0) goto L565;
	    }
            std::printf(" PARAMETER NUMBER %3d NOT VARIABLE. IGNORED.\n",ilax);
	    continue;
	}
L565:
//*-*-                                        calculate one pair of M E s
	ilax2 = 0;
	mnmnot(ilax, ilax2, val2pl, val2mi);
	if (fLnewmn) goto L650;
//*-*-                                         update NGOOD and NBAD
//        iin = fNiofex[ilax-1];
	iin = m_userParameterIdToInternalId[ilax];
	if (fErp[iin-1] > 0) ++ngood;
	else                   ++nbad;
	if (fErn[iin-1] < 0) ++ngood;
	else                   ++nbad;
    }
//*-*-                                          end of loop . . . . . . .
//*-*-                                       . . . . printout final values .
    fCfrom  = "MINOS   ";
    fNfcnfr = nfcnmi;
    fCstatu = "UNCHANGED ";
    if (ngood == 0 && nbad == 0) goto L700;
    if (ngood > 0 && nbad == 0)  fCstatu = "SUCCESSFUL";
    if (ngood == 0 && nbad > 0)  fCstatu = "FAILURE   ";
    if (ngood > 0 && nbad > 0)   fCstatu = "PROBLEMS  ";
    if (fISW[4] >= 0)    mnprin(4, fAmin);
    if (fISW[4] >= 2)    mnmatu(0);
    return;
//*-*-                                       . . . new minimum found . . . .
L650:
    fCfrom  = "MINOS   ";
    fNfcnfr = nfcnmi;
    fCstatu = "NEW MINIMU";
    if (fISW[4] >= 0) mnprin(4, fAmin);
    std::printf(" NEW MINIMUM FOUND.  GO BACK TO MINIMIZATION STEP.\n");
    std::printf(" =================================================\n");
    std::printf("                                                  V\n");
    std::printf("                                                  V\n");
    std::printf("                                                  V\n");
    std::printf("                                               VVVVVVV\n");
    std::printf("                                                VVVVV\n");
    std::printf("                                                 VVV\n");
    std::printf("                                                  V\n");
    std::printf("\n");
    return;
L700:
    std::printf(" THERE ARE NO MINOS ERRORS TO CALCULATE.\n");
} /* mnmnos_ */

//______________________________________________________________________________
void URMinuit::mnmnot(Int_urt ilax, Int_urt ilax2, Double_urt &val2pl, Double_urt &val2mi)
{
//*-*-*-*-*-*Performs a MINOS error analysis on one parameter*-*-*-*-*-*-*-*-*
//*-*        ================================================
//*-*        The parameter ILAX is varied, and the minimum of the
//*-*        function with respect to the other parameters is followed
//*-*        until it crosses the value FMIN+UP.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* System generated locals */
    Int_urt i__1;

    /* Local variables */
    Double_urt delu, aopt, eros;
    Double_urt abest, xunit, dc, ut, sigsav, du1;
    Double_urt fac, sig, sav;
    Int_urt marc, isig, mpar, ndex, imax, indx, ierr, i, j;
    Int_urt iercr, it, istrav, nfmxin, nlimit, isw2, isw4;
    string csig;

//*-*-                                       . . save and prepare start vals
    isw2    = fISW[1];
    isw4    = fISW[3];
    sigsav  = fEDM;
    istrav  = fIstrat;
    dc      = fDcovar;
    fLnewmn = kurFALSE;
    fApsi   = fEpsi*.5;
    abest   = fAmin;
    mpar    = fNpar;
    nfmxin  = fNfcnmx;
    for (i = 1; i <= mpar; ++i) { fXt[i-1] = fX[i-1]; }
    i__1 = mpar*(mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) { fVthmat[j-1] = fVhmat[j-1]; }
    for (i = 1; i <= mpar; ++i) {
	fMNOTgcc[i-1] = fGlobcc[i-1];
	fMNOTw[i-1]   = fWerr[i-1];
    }
//    it = fNiofex[ilax-1];
    it = m_userParameterIdToInternalId[ilax];
    fErp[it-1] = 0;
    fErn[it-1] = 0;
    mninex(fXt);
//    ut = fU[ilax-1];
//    if (fNvarl[ilax-1] == 1) {
    ut = m_userParameterValue[ilax];
    if (m_userParameterFlag[ilax] == 1) {
	fAlim[ilax-1] = ut - fMNOTw[it-1]*100;
	fBlim[ilax-1] = ut + fMNOTw[it-1]*100;
    }
    ndex  = it*(it + 1) / 2;
    xunit = URMath::Sqrt(fUp / fVthmat[ndex-1]);
    marc  = 0;
    for (i = 1; i <= mpar; ++i) {
	if (i == it) continue;
	++marc;
	imax = URMath::Max(it,i);
	indx = imax*(imax-1) / 2 + URMath::Min(it,i);
	fMNOTxdev[marc-1] = xunit*fVthmat[indx-1];
    }
//*-*-                          fix the parameter in question
    mnfixp(it-1, ierr);
    if (ierr > 0) {
	std::printf(" MINUIT ERROR. CANNOT FIX PARAMETER%4d    INTERNAL%3d\n",ilax,it);
	goto L700;
    }
//*-*-                      . . . . . Nota Bene: from here on, NPAR=MPAR-1
//*-*-     Remember: MNFIXP squeezes IT out of X, XT, WERR, and VHMAT,
//*-*-                                                   not W, VTHMAT
    for (isig = 1; isig <= 2; ++isig) {
	if (isig == 1) {
	    sig  = 1;
	    csig = "POSI";
	} else {
	    sig  = -1;
	    csig = "NEGA";
	}
//*-*-                                       . sig=sign of error being calcd
	if (fISW[4] > 1) {
//          std::printf(" DETERMINATION OF %sTIVE MINOS ERROR FOR PARAMETER%d\n"
//            ,csig.c_str(),ilax
//            ,(fCpnam[ilax-1]).c_str());
           
            std::printf(" DETERMINATION OF %sTIVE MINOS ERROR FOR PARAMETER%d (%s)\n"
                            ,csig.c_str(),ilax
                            ,(m_userParameterName[ilax]).c_str());
	}
	if (fISW[1] <= 0) {
	    mnwarn("D", "MINOS", "NO COVARIANCE MATRIX.");
	}
	nlimit     = fNfcn + nfmxin;
	fIstrat    = URMath::Max(istrav-1,0);
	du1        = fMNOTw[it-1];
//	fU[ilax-1] = ut + sig*du1;
//	fU[ilax-1] = URMath::Min(fU[ilax-1],fBlim[ilax-1]);
//	fU[ilax-1] = URMath::Max(fU[ilax-1],fAlim[ilax-1]);
        m_userParameterValue[ilax] = ut + sig*du1;
        m_userParameterValue[ilax] = URMath::Min(m_userParameterValue[ilax],fBlim[ilax-1]);
        m_userParameterValue[ilax] = URMath::Max(m_userParameterValue[ilax],fAlim[ilax-1]);
//        delu = fU[ilax-1] - ut;
        delu = m_userParameterValue[ilax] - ut;
//*-*-        stop if already at limit with negligible step size
//        if (URMath::Abs(delu) / (URMath::Abs(ut) + URMath::Abs(fU[ilax-1])) < fEpsmac) goto L440;
	if (URMath::Abs(delu) / (URMath::Abs(ut) + URMath::Abs(m_userParameterValue[ilax])) < fEpsmac) goto L440;
	fac = delu / fMNOTw[it-1];
	for (i = 1; i <= fNpar; ++i) {
	    fX[i-1] = fXt[i-1] + fac*fMNOTxdev[i-1];
	}
	if (fISW[4] > 1) {
//	    std::printf(" PARAMETER%4d SET TO%11.3e + %10.3e = %12.3e\n",ilax,ut,delu,fU[ilax-1]);
            std::printf(" PARAMETER%4d SET TO%11.3e + %10.3e = %12.3e\n",ilax,ut,delu,m_userParameterValue[ilax]);
	}
//*-*-                                       loop to hit AMIN+UP
	fKe1cr  = ilax;
	fKe2cr  = 0;
//	fXmidcr = fU[ilax-1];
        fXmidcr = m_userParameterValue[ilax];
	fXdircr = delu;

	fAmin = abest;
	fNfcnmx = nlimit - fNfcn;
	mncros(aopt, iercr);
	if (abest - fAmin > fUp*.01) goto L650;
	if (iercr == 1) goto L440;
	if (iercr == 2) goto L450;
	if (iercr == 3) goto L460;
//*-*-                                       . error successfully calculated
	eros = fXmidcr - ut + aopt*fXdircr;
	if (fISW[4] > 1) {
//          std::printf("         THE %4sTIVE MINOS ERROR OF PARAMETER%3d  %10s, IS %12.4e\n"
//            ,csig.c_str(),ilax
//            ,(fCpnam[ilax-1]).c_str(),eros);
            std::printf("         THE %4sTIVE MINOS ERROR OF PARAMETER%3d  %10s, IS %12.4e\n"
                           ,csig.c_str(),ilax
                           ,(m_userParameterName[ilax]).c_str(),eros);
	}
	goto L480;
//*-*-                                       . . . . . . . . failure returns
L440:
	if (fISW[4] >= 1) {
//          std::printf("     THE %4sTIVE MINOS ERROR OF PARAMETER%3d, %s EXCEEDS ITS LIMIT.\n"
//            ,csig.c_str(),ilax
//            ,(fCpnam[ilax-1]).c_str());
            std::printf("     THE %4sTIVE MINOS ERROR OF PARAMETER%3d, %s EXCEEDS ITS LIMIT.\n"
                              ,csig.c_str(),ilax
                              ,(m_userParameterName[ilax]).c_str());
	}
	eros = fUndefi;
	goto L480;
L450:
	if (fISW[4] >= 1) {
	    std::printf("         THE %4sTIVE MINOS ERROR%4d REQUIRES MORE THAN%5d FUNCTION CALLS.\n"
                         ,csig.c_str(),ilax,nfmxin);
	}
	eros = 0;
	goto L480;
L460:
	if (fISW[4] >= 1) {
	    std::printf("                         %4sTIVE MINOS ERROR NOT CALCULATED FOR PARAMETER%d\n"
                         ,csig.c_str(),ilax);
	}
	eros = 0;

L480:
	if (fISW[4] > 1) {
	    std::printf("     **************************************************************************\n");
	}
	if (sig < 0) {
	    fErn[it-1] = eros;
//	    if (ilax2 > 0 && ilax2 <= fNu) val2mi = fU[ilax2-1];
            if (ilax2 > 0 && ilax2 <= fNu) val2mi = m_userParameterValue[ilax2];
	} else {
	    fErp[it-1] = eros;
//	    if (ilax2 > 0 && ilax2 <= fNu) val2pl = fU[ilax2-1];
            if (ilax2 > 0 && ilax2 <= fNu) val2pl = m_userParameterValue[ilax2];
	}
    }
//*-*-                                       . . parameter finished. reset v
//*-*-                      normal termination */
    fItaur = 1;
    mnfree(1);
    i__1 = mpar*(mpar + 1) / 2;
    for (j = 1; j <= i__1; ++j) { fVhmat[j-1] = fVthmat[j-1]; }
    for (i = 1; i <= mpar; ++i) {
	fWerr[i-1]   = fMNOTw[i-1];
	fGlobcc[i-1] = fMNOTgcc[i-1];
	fX[i-1]      = fXt[i-1];
    }
    mninex(fX);
    fEDM    = sigsav;
    fAmin   = abest;
    fISW[1] = isw2;
    fISW[3] = isw4;
    fDcovar = dc;
    goto L700;
//*-*-                      new minimum
L650:
    fLnewmn = kurTRUE;
    fISW[1] = 0;
    fDcovar = 1;
    fISW[3] = 0;
//    sav     = fU[ilax-1];
    sav     = m_userParameterValue[ilax];
    fItaur  = 1;
    mnfree(1);
//    fU[ilax-1] = sav;
    m_userParameterValue[ilax] = sav;
    mnexin(fX);
    fEDM = fBigedm;
//*-*-                      in any case
L700:
    fItaur  = 0;
    fNfcnmx = nfmxin;
    fIstrat = istrav;
} /* mnmnot_ */

//______________________________________________________________________________
void URMinuit::mnparm(Int_urt k1, const string& cnamj, Double_urt uk, Double_urt wk, Double_urt a, Double_urt b, Int_urt &ierflg)
{
//*-*-*-*-*-*-*-*-*Implements one parameter definition*-*-*-*-*-*-*-*-*-*-*-*
//*-*              ===================================
//*-*        Called from MNPARS and user-callable
//*-*    Implements one parameter definition, that is:
//*-*          K     (external) parameter number
//*-*          CNAMK parameter name
//*-*          UK    starting value
//*-*          WK    starting step size or uncertainty
//*-*          A, B  lower and upper physical parameter limits
//*-*    and sets up (updates) the parameter lists.
//*-*    Output: IERFLG=0 if no problems
//*-*                  >0 if MNPARM unable to implement definition
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt vplu, a_small, gsmin, pinti, vminu, danger, sav, sav2;
    Int_urt ierr, kint, in, ix, ktofix, lastin, kinfix, nvl;
    string cnamk, chbufi;

//    Int_urt k = k1+1;
    Int_urt k = k1;
    cnamk   = cnamj;
    kint    = fNpar;
    if (k < 1 || k > fMaxext) {
//*-*-                    parameter number exceeds allowed maximum value
	//std::printf(" MINUIT USER ERROR.  PARAMETER NUMBER IS %3d  ALLOWED RANGE IS ONE TO %4d\n",k,fMaxext);
        *m_logStream << " MINUIT USER ERROR.  PARAMETER NUMBER IS " << setw(3) << k
                    << "  ALLOWED RANGE IS ONE TO " << setw(4) << fMaxext << '\n';
	goto L800;
    }
//*-*-                    normal parameter request
    ktofix = 0;
//    if (fNvarl[k-1] < 0) goto L50;
    if (m_userParameterFlag[k] < 0) goto L50;
//*-*-        previously defined parameter is being redefined
//*-*-                                    find if parameter was fixed
    for (ix = 1; ix <= fNpfix; ++ix) {
	if (fIpfix[ix-1] == k) ktofix = k;
    }
    if (ktofix > 0) {
	mnwarn("W", "PARAM DEF", "REDEFINING A FIXED PARAMETER.");
	if (kint >= fMaxint) {
	    //std::printf(" CANNOT RELEASE. MAX NPAR EXCEEDED.\n");
            *m_logStream << " CANNOT RELEASE. MAX NPAR EXCEEDED.\n";
	    goto L800;
	}
	mnfree(-k);
    }
//*-*-                      if redefining previously variable parameter
//    if (fNiofex[k-1] > 0) kint = fNpar - 1;
    if (m_userParameterIdToInternalId[k] > 0) kint = fNpar - 1;
L50:

//*-*-                                     . . .print heading
    if (fLphead && fISW[4] >= 0) {
	//std::printf(" PARAMETER DEFINITIONS:\n");
	//std::printf("    NO.   NAME         VALUE      STEP SIZE      LIMITS\n");
        *m_logStream << " PARAMETER DEFINITIONS:\n";
        *m_logStream << "    NO.   NAME                 VALUE      STEP SIZE      LIMITS\n";
	fLphead = kurFALSE;
    }
    if (wk > 0) goto L122;
//*-*-                                       . . .constant parameter . . . .
    if (fISW[4] >= 0) {
        //std::printf(" %5d %-10s %13.5e  constant\n",k,cnamk.c_str(),uk);
        *m_logStream << setw(6) << k << " " << setw(20) << cnamk << setw(13) << setprecision(5) << uk << "  constant\n";
    }
    nvl = 0;
    goto L200;
L122:
    if (a == 0 && b == 0) {
//*-*-                                     variable parameter without limits
	nvl = 1;
	if (fISW[4] >= 0) {
            //std::printf(" %5d %-10s %13.5e%13.5e     no limits\n",k,cnamk.c_str(),uk,wk);
           *m_logStream << setw(6) << k << " " << setw(20) << cnamk 
                        << setw(13) << setprecision(5) << uk
                        << setw(13) << setprecision(5) << wk << "     no limits\n";
	}
    } else {
//*-*-                                        variable parameter with limits
	nvl = 4;
	fLnolim = kurFALSE;
	if (fISW[4] >= 0) {
            //std::printf(" %5d %-10s %13.5e%13.5e  %13.5e%13.5e\n",k,cnamk.c_str(),uk,wk,a,b);
           *m_logStream << setw(6) << k << " " << setw(20) << cnamk 
              << setw(13) << setprecision(5) << uk
              << setw(13) << setprecision(5) << wk
              << setw(13) << setprecision(5) << a
              << setw(13) << setprecision(5) << b << '\n';
	}
    }
//*-*-                            . . request for another variable parameter
    ++kint;
    if (kint > fMaxint) {
	//std::printf(" MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.\n");
        *m_logStream << " MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.\n";
	goto L800;
    }
    if (nvl == 1) goto L200;
    if (a == b) {
	std::printf(" USER ERROR IN MINUIT PARAMETER\n");
	std::printf(" DEFINITION\n");
	std::printf(" UPPER AND LOWER LIMITS EQUAL.\n");
        *m_logStream << " USER ERROR IN MINUIT PARAMETER\n" << " DEFINITION\n" << " UPPER AND LOWER LIMITS EQUAL.\n";
	goto L800;
    }
    if (b < a) {
	sav = b;
	b = a;
	a = sav;
	mnwarn("W", "PARAM DEF", "PARAMETER LIMITS WERE REVERSED.");
	if (fLwarn) fLphead = kurTRUE;
    }
    if (b - a > 1e7) {
       ostringstream warning;
       warning << "LIMITS ON PARAM " << k << " TOO FAR APART.";
	mnwarn("W", "PARAM DEF", warning.str().c_str());
	if (fLwarn) fLphead = kurTRUE;
    }
    danger = (b - uk)*(uk - a);
    if (danger < 0) {
	mnwarn("W", "PARAM DEF", "STARTING VALUE OUTSIDE LIMITS.");
    }
    if (danger == 0) {
	mnwarn("W", "PARAM DEF", "STARTING VALUE IS AT LIMIT.");
    }
L200:
//*-*-                          . . . input OK, set values, arrange lists,
//*-*-                                   calculate step sizes GSTEP, DIRIN
    fCfrom      = "PARAMETR";
    fNfcnfr     = fNfcn;
    fCstatu     = "NEW VALUES";
    fNu         = URMath::Max(fNu,k);
//    fCpnam[k-1] = cnamk;
//    fU[k-1]     = uk;
    m_userParameterName[k] = cnamk;
    m_userParameterValue[k]     = uk;
    fAlim[k-1]  = a;
    fBlim[k-1]  = b;
//    fNvarl[k-1] = nvl;
    m_userParameterFlag[k] = nvl;
    mnrset(1);
//*-*-                            K is external number of new parameter
//*-*-          LASTIN is the number of var. params with ext. param. no.< K
    lastin = 0;
//    for (ix = 1; ix <= k-1; ++ix) { if (fNiofex[ix-1] > 0) ++lastin; }
    for (ix = 1; ix <= k-1; ++ix) { if (m_userParameterIdToInternalId[ix] > 0) ++lastin; }
//*-*-                KINT is new number of variable params, NPAR is old
    if (kint == fNpar) goto L280;
    if (kint > fNpar) {
//*-*-                         insert new variable parameter in list
	for (in = fNpar; in >= lastin + 1; --in) {
	    ix            = fNexofi[in-1];
//           fNiofex[ix-1] = in + 1;
	    m_userParameterIdToInternalId[ix] = in + 1;
	    fNexofi[in]   = ix;
	    fX[in]        = fX[in-1];
	    fXt[in]       = fXt[in-1];
	    fDirin[in]    = fDirin[in-1];
	    fG2[in]       = fG2[in-1];
	    fGstep[in]    = fGstep[in-1];
	}
    } else {
//*-*-                         remove variable parameter from list
	for (in = lastin + 1; in <= kint; ++in) {
	    ix            = fNexofi[in];
//           fNiofex[ix-1] = in;
	    m_userParameterIdToInternalId[ix] = in;
	    fNexofi[in-1] = ix;
	    fX[in-1]      = fX[in];
	    fXt[in-1]     = fXt[in];
	    fDirin[in-1]  = fDirin[in];
	    fG2[in-1]     = fG2[in];
	    fGstep[in-1]  = fGstep[in];
	}
    }
L280:
    ix = k;
//    fNiofex[ix-1] = 0;
    m_userParameterIdToInternalId[ix] = 0;
    fNpar = kint;
//*-*-                                      lists are now arranged . . . .
    if (nvl > 0) {
	in            = lastin + 1;
	fNexofi[in-1] = ix;
//        fNiofex[ix-1] = in;
	m_userParameterIdToInternalId[ix] = in;
//	sav           = fU[ix-1];
//        mnpint(sav, ix-1, pinti);
        sav           = m_userParameterValue[ix];
	mnpint(sav, ix, pinti);
	fX[in-1]    = pinti;
	fXt[in-1]   = fX[in-1];
	fWerr[in-1] = wk;
	sav2        = sav + wk;
//	mnpint(sav2, ix-1, pinti);
        mnpint(sav2, ix, pinti);
	vplu = pinti - fX[in-1];
	sav2 = sav - wk;
//	mnpint(sav2, ix-1, pinti);
        mnpint(sav2, ix, pinti);
	vminu = pinti - fX[in-1];
	fDirin[in-1] = (URMath::Abs(vplu) + URMath::Abs(vminu))*.5;
	fG2[in-1] = fUp*2 / (fDirin[in-1]*fDirin[in-1]);
	gsmin = fEpsma2*8*URMath::Abs(fX[in-1]);
	fGstep[in-1] = URMath::Max(gsmin,fDirin[in-1]*.1);
	if (fAmin != fUndefi) {
	    a_small      = URMath::Sqrt(fEpsma2*(fAmin + fUp) / fUp);
	    fGstep[in-1] = URMath::Max(gsmin,a_small*fDirin[in-1]);
	}
	fGrd[in-1] = fG2[in-1]*fDirin[in-1];
//*-*-                  if parameter has limits
//        if (fNvarl[k-1] > 1) {
	if (m_userParameterFlag[k] > 1) {
	    if (fGstep[in-1] > .5) fGstep[in-1] = .5;
	    fGstep[in-1] = -fGstep[in-1];
	}
    }
    if (ktofix > 0) {
//       kinfix = fNiofex[ktofix-1];
	kinfix = m_userParameterIdToInternalId[ktofix];
	if (kinfix > 0) mnfixp(kinfix-1, ierr);
	if (ierr > 0)   goto L800;
    }
    ierflg = 0;
    return;
    *m_logStream << endl;
//*-*-                  error on input, unable to implement request  . . . .
L800:
    ierflg = 1;
    *m_logStream << endl;
} /* mnparm_ */

//______________________________________________________________________________
void URMinuit::mnpars(const string &crdbuf, Int_urt &icondn)
{
//*-*-*-*-*-*-*-*Implements one parameter definition*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*            =========== =======================
//*-*        Called from MNREAD and user-callable
//*-*    Implements one parameter definition, that is:
//*-*       parses the string CRDBUF and calls MNPARM
//*-*
//*-* output conditions:
//*-*        ICONDN = 0    all OK
//*-*        ICONDN = 1    error, attempt to define parameter is ignored
//*-*        ICONDN = 2    end of parameter definitions
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt a=0, b=0, fk=0, uk=0, wk=0, xk=0;
    Int_urt ierr, kapo1, kapo2;
    Int_urt k, llist, ibegin, lenbuf, istart, lnc, icy;
    string cnamk, comand, celmnt, ctemp;
    char stmp[128];

    lenbuf = strlen(crdbuf.c_str());
//*-*-                    find out whether fixed or free-field format
    kapo1 = strspn(crdbuf.c_str(), "'");
    if (kapo1 == 0) goto L150;
    kapo2 = strspn(crdbuf.c_str() + kapo1, "'");
    if (kapo2 == 0) goto L150;
//*-*-         new (free-field) format
    kapo2 += kapo1;
//*-*-                            skip leading blanks if any
    for (istart = 1; istart <= kapo1-1; ++istart) {
	if (crdbuf.substr(istart-1,1) != " ") goto L120;
    }
    goto L210;
L120:
//*-*-                              parameter number integer
    celmnt = crdbuf.substr(istart-1, kapo1-istart);
    scanf(celmnt.c_str(),fk);
    k = Int_urt(fk);
    if (k <= 0) goto L210;
    cnamk = "PARAM " + celmnt;
    if (kapo2 - kapo1 > 1) {
	cnamk = crdbuf.substr(kapo1, kapo2-1-kapo1);
    }
//*-*  special handling if comma or blanks and a comma follow 'name'
    for (icy = kapo2 + 1; icy <= lenbuf; ++icy) {
	if (crdbuf.substr(icy-1,1) == ",") goto L139;
	if (crdbuf.substr(icy-1,1) != " ") goto L140;
    }
    uk = 0;
    wk = 0;
    a  = 0;
    b  = 0;
    goto L170;
L139:
    ++icy;
L140:
    ibegin = icy;
    ctemp = crdbuf.substr(ibegin-1,lenbuf-ibegin);
    mncrck(ctemp, 20, comand, lnc, fMaxpar, fPARSplist, llist, ierr, fIsyswr);
    if (ierr > 0) goto L180;
    uk = fPARSplist[0];
    wk = 0;
    if (llist >= 2) wk = fPARSplist[1];
    a = 0;
    if (llist >= 3) a = fPARSplist[2];
    b = 0;
    if (llist >= 4) b = fPARSplist[3];
    goto L170;
//*-*-         old (fixed-field) format
L150:
    scanf(crdbuf.c_str(),xk,stmp,uk,wk,a,b);
    cnamk = stmp;
    k = Int_urt(xk);
    if (k == 0)    goto L210;
//*-*-         parameter format cracked, implement parameter definition
L170:
//       mnparm(k-1, cnamk, uk, wk, a, b, ierr);
    mnparm(k, cnamk, uk, wk, a, b, ierr);
    icondn = ierr;
    return;
//*-*-         format or other error
L180:
    icondn = 1;
    return;
//*-*-       end of data
L210:
    icondn = 2;
} /* mnpars_ */

//______________________________________________________________________________
void URMinuit::mnpfit(Double_urt *parx2p, Double_urt *pary2p, Int_urt npar2p, Double_urt *coef2p, Double_urt &sdev2p)
{
//*-*-*-*-*-*-*-*-*-*To fit a parabola to npar2p points*-*-*-*-*-*-*-*-*-*-*
//*-*                ==================================
//*-*   npar2p   no. of points
//*-*   parx2p(i)   x value of point i
//*-*   pary2p(i)   y value of point i
//*-*
//*-*   coef2p(1...3)  coefficients of the fitted parabola
//*-*   y=coef2p(1) + coef2p(2)*x + coef2p(3)*x**2
//*-*   sdev2p= variance
//*-*   method : chi**2 = min equation solved explicitly
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt a, f, s, t, y, s2, x2, x3, x4, y2, cz[3], xm, xy, x2y;
    x2 = x3 = 0;
    Int_urt i;

    /* Parameter adjustments */
    --coef2p;
    --pary2p;
    --parx2p;

    /* Function Body */
    for (i = 1; i <= 3; ++i) { cz[i-1] = 0; }
    sdev2p = 0;
    if (npar2p < 3) goto L10;
    f = (Double_urt) (npar2p);
//*-* --- center x values for reasons of machine precision
    xm  = 0;
    for (i = 1; i <= npar2p; ++i) { xm += parx2p[i]; }
    xm /= f;
    x2  = 0;
    x3  = 0;
    x4  = 0;
    y   = 0;
    y2  = 0;
    xy  = 0;
    x2y = 0;
    for (i = 1; i <= npar2p; ++i) {
	s    = parx2p[i] - xm;
	t    = pary2p[i];
	s2   = s*s;
	x2  += s2;
	x3  += s*s2;
	x4  += s2*s2;
	y   += t;
	y2  += t*t;
	xy  += s*t;
	x2y += s2*t;
    }
    a = (f*x4 - x2*x2)*x2 - f*(x3*x3);
    if (a == 0) goto L10;
    cz[2] = (x2*(f*x2y - x2*y) - f*x3*xy) / a;
    cz[1] = (xy - x3*cz[2]) / x2;
    cz[0] = (y - x2*cz[2]) / f;
    if (npar2p == 3) goto L6;
    sdev2p = y2 - (cz[0]*y + cz[1]*xy + cz[2]*x2y);
    if (sdev2p < 0) sdev2p = 0;
    sdev2p /= f - 3;
L6:
    cz[0] += xm*(xm*cz[2] - cz[1]);
    cz[1] -= xm*2*cz[2];
L10:
    for (i = 1; i <= 3; ++i) { coef2p[i] = cz[i-1]; }
} /* mnpfit_ */

//______________________________________________________________________________
void URMinuit::mnpint(Double_urt &pexti, Int_urt i1, Double_urt &pinti)
{
//*-*-*-*-*-*-*Calculates the internal parameter value PINTI*-*-*-*-*-*-*-*
//*-*          =============================================
//*-*        corresponding  to the external value PEXTI for parameter I.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt a, alimi, blimi, yy, yy2;
    Int_urt igo;
    string chbuf2, chbufi;

//    Int_urt i = i1+1;
    Int_urt i = i1;
    pinti   = pexti;
//    igo     = fNvarl[i-1];
    igo     = m_userParameterFlag[i];
    if (igo == 4) {
//*-* --                          there are two limits
	alimi = fAlim[i-1];
	blimi = fBlim[i-1];
	yy = (pexti - alimi)*2 / (blimi - alimi) - 1;
	yy2 = yy*yy;
	if (yy2 >= 1 - fEpsma2) {
	    if (yy < 0) {
		a      = fVlimlo;
		chbuf2 = " IS AT ITS LOWER ALLOWED LIMIT.";
	    } else {
		a      = fVlimhi;
		chbuf2 = " IS AT ITS UPPER ALLOWED LIMIT.";
	    }
	    pinti   = a;
	    pexti   = alimi + (blimi - alimi)*.5*(URMath::Sin(a) + 1);
	    fLimset = kurTRUE;
	    if (yy2 > 1) chbuf2 = " BROUGHT BACK INSIDE LIMITS.";
            ostringstream warning;
            warning << "VARIABLE " << i << " " << chbuf2;
	    mnwarn("W", fCfrom.c_str(), warning.str().c_str());
	} else {
	    pinti = URMath::ASin(yy);
	}
    }
} /* mnpint_ */

//______________________________________________________________________________
void URMinuit::mnplot(Double_urt *xpt, Double_urt *ypt, char* chpt, Int_urt nxypt, Int_urt npagwd, Int_urt npagln)
{
 //*-*-*-*Plots points in array xypt onto one page with labelled axes*-*-*-*-*
 //*-*    ===========================================================
 //*-*        NXYPT is the number of points to be plotted
 //*-*        XPT(I) = x-coord. of ith point
 //*-*        YPT(I) = y-coord. of ith point
 //*-*        CHPT(I) = character to be plotted at this position
 //*-*        the input point arrays XPT, YPT, CHPT are destroyed.
 //*-*
 //*-*
 //*-*   
 //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
     
     static string cdot   = ".";
     static string cslash = "/";

     /* Local variables */
     Double_urt xmin, ymin, xmax, ymax, savx, savy, yprt;
     Double_urt bwidx, bwidy, xbest, ybest, ax, ay, bx, by;
     Double_urt xvalus[12], any, dxx, dyy;
     Int_urt iten, i, j, k, maxnx, maxny, iquit, ni, linodd;
     Int_urt nxbest, nybest, km1, ibk, isp1, nx, ny, ks, ix;
     string chmess, ctemp;
     Bool_urt overpr;
     char cline[120];
     char chsav, chbest;

     /* Function Body */
 //*-*  Computing MIN
     maxnx = URMath::Min(npagwd-20,100);
     if (maxnx < 10) maxnx = 10;
     maxny = npagln;
     if (maxny < 10) maxny = 10;
     if (nxypt <= 1) return;
     xbest  = xpt[0];
     ybest  = ypt[0];
     chbest = chpt[0];
 //*-*-        order the points by decreasing y
     km1 = nxypt - 1;
     for (i = 1; i <= km1; ++i) {
 	iquit = 0;
 	ni    = nxypt - i;
 	for (j = 1; j <= ni; ++j) {
 	    if (ypt[j-1] > ypt[j]) continue;
 	    savx     = xpt[j-1];
 	    xpt[j-1] = xpt[j];
 	    xpt[j]   = savx;
 	    savy     = ypt[j-1];
 	    ypt[j-1] = ypt[j];
 	    ypt[j]   = savy;
 	    chsav    = chpt[j-1];
 	    chpt[j-1]= chpt[j];
 	    chpt[j]  = chsav;
 	    iquit    = 1;
 	}
 	if (iquit == 0) break;
     }
 //*-*-        find extreme values
     xmax = xpt[0];
     xmin = xmax;
     for (i = 1; i <= nxypt; ++i) {
 	if (xpt[i-1] > xmax) xmax = xpt[i-1];
 	if (xpt[i-1] < xmin) xmin = xpt[i-1];
     }
     dxx   = (xmax - xmin)*.001;
     xmax += dxx;
     xmin -= dxx;
     mnbins(xmin, xmax, maxnx, xmin, xmax, nx, bwidx);
     ymax = ypt[0];
     ymin = ypt[nxypt-1];
     if (ymax == ymin) ymax = ymin + 1;
     dyy   = (ymax - ymin)*.001;
     ymax += dyy;
     ymin -= dyy;
     mnbins(ymin, ymax, maxny, ymin, ymax, ny, bwidy);
     any = (Double_urt) ny;
 //*-*-        if first point is blank, it is an 'origin'
     if (chbest == ' ') goto L50;
     xbest = (xmax + xmin)*.5;
     ybest = (ymax + ymin)*.5;
 L50:
 //*-*-        find scale constants
     ax = 1 / bwidx;
     ay = 1 / bwidy;
     bx = -ax*xmin + 2;
     by = -ay*ymin - 2;
 //*-*-        convert points to grid positions
     for (i = 1; i <= nxypt; ++i) {
 	xpt[i-1] = ax*xpt[i-1] + bx;
 	ypt[i-1] = any - ay*ypt[i-1] - by;
     }
     nxbest = Int_urt((ax*xbest + bx));
     nybest = Int_urt((any - ay*ybest - by));
 //*-*-        print the points
     ny += 2;
     nx += 2;
     isp1 = 1;
     linodd = 1;
     overpr = kurFALSE;
     for (i = 1; i <= ny; ++i) {
 	for (ibk = 1; ibk <= nx; ++ibk) { cline[ibk-1] = ' '; }
 	cline[nx] = '\0';
 	cline[nx+1] = '\0';
 	cline[0]        = '.';
 	cline[nx-1]     = '.';
 	cline[nxbest-1] = '.';
 	if (i != 1 && i != nybest && i != ny) goto L320;
 	for (j = 1; j <= nx; ++j) { cline[j-1] = '.'; }
 L320:
 	yprt = ymax - Double_urt(i-1)*bwidy;
 	if (isp1 > nxypt) goto L350;
 //*-*-        find the points to be plotted on this line
 	for (k = isp1; k <= nxypt; ++k) {
 	    ks = Int_urt(ypt[k-1]);
 	    if (ks > i) goto L345;
 	    ix = Int_urt(xpt[k-1]);
 	    if (cline[ix-1] == '.')   goto L340;
 	    if (cline[ix-1] == ' ') goto L340;
 	    if (cline[ix-1] == chpt[k-1])   continue;
 	    overpr = kurTRUE;
 //*-*-        OVERPR is true if one or more positions contains more than
 //*-*-           one point
 	    cline[ix-1] = '&';
 	    continue;
 L340:
 	    cline[ix-1] = chpt[k-1];
 	}
 	isp1 = nxypt + 1;
 	goto L350;
 L345:
 	isp1 = k;
 L350:
 	if (linodd == 1 || i == ny) goto L380;
 	linodd = 1;
 	ctemp  = cline;
 	std::printf("                  %s\n",ctemp.c_str());
 	goto L400;
 L380:
 	ctemp = cline;
 	std::printf(" %14.7g ..%s\n",yprt,ctemp.c_str());
 	linodd = 0;
 L400:
 	;
     }
 //*-*-        print labels on x-axis every ten columns
     for (ibk = 1; ibk <= nx; ++ibk) {
 	cline[ibk-1] = ' ';
 	if (ibk % 10 == 1) cline[ibk-1] = '/';
     }
     std::printf("                  %s\n",cline);

     for (ibk = 1; ibk <= 12; ++ibk) {
 	xvalus[ibk-1] = xmin + Double_urt(ibk-1)*10*bwidx;
     }
     std::printf("           ");
     iten = (nx + 9) / 10;
     for (ibk = 1; ibk <= iten; ++ibk) {
         std::printf(" %9.4g", xvalus[ibk-1]);
     }
     chmess = " ";
     if (overpr) chmess = "   Overprint character is &";
     std::printf("                         ONE COLUMN=%13.7g%s\n",bwidx,chmess.c_str());
} /* mnplot_ */

//______________________________________________________________________________
void URMinuit::mnpout(Int_urt iuext1,  string &chnam, Double_urt &val, Double_urt &err, Double_urt &xlolim, Double_urt &xuplim, Int_urt &iuint) const
{
//*-*-*-*Provides the user with information concerning the current status*-*-*
//*-*    ================================================================
//*-*          of parameter number IUEXT. Namely, it returns:
//*-*        CHNAM: the name of the parameter
//*-*        VAL: the current (external) value of the parameter
//*-*        ERR: the current estimate of the parameter uncertainty
//*-*        XLOLIM: the lower bound (or zero if no limits)
//*-*        XUPLIM: the upper bound (or zero if no limits)
//*-*        IUINT: the internal parameter number (or zero if not variable,
//*-*           or negative if undefined).
//*-*  Note also:  If IUEXT is negative, then it is -internal parameter
//*-*           number, and IUINT is returned as the EXTERNAL number.
//*-*     Except for IUINT, this is exactly the inverse of MNPARM
//*-*     User-called
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Int_urt iint, iext, nvl;

//    Int_urt iuext = iuext1 + 1;
    Int_urt iuext = iuext1;
    xlolim = 0;
    xuplim = 0;
    err    = 0;
    if (iuext == 0) goto L100;
    if (iuext < 0) {
//*-*-                  internal parameter number specified
	iint  = -(iuext);
	if (iint > fNpar) goto L100;
	iext  = fNexofi[iint-1];
	iuint = iext;
    } else {
//*-*-                   external parameter number specified
	iext = iuext;
	if (iext == 0) goto L100;
	if (iext > fNu) goto L100;
//	iint  = fNiofex[iext-1];
        iint  = m_userParameterIdToInternalId[iext];
	iuint = iint;
    }
//*-*-                    in both cases
//    nvl = fNvarl[iext-1];
    nvl = m_userParameterFlag[iext];
    if (nvl < 0) goto L100;
//    chnam = fCpnam[iext-1];
//    val   = fU[iext-1];
    chnam = m_userParameterName[iext];
    val   = m_userParameterValue[iext];
    if (iint > 0) err = fWerr[iint-1];
    if (nvl == 4) {
	xlolim = fAlim[iext-1];
	xuplim = fBlim[iext-1];
    }
    return;
//*-*-               parameter is undefined
L100:
    iuint = -1;
    chnam = "undefined";
    val = 0;
} /* mnpout_ */

//______________________________________________________________________________
void URMinuit::mnprin(Int_urt inkode, Double_urt fval)
{
//*-*-*-*Prints the values of the parameters at the time of the call*-*-*-*-*
//*-*    ===========================================================
//*-*        also prints other relevant information such as function value,
//*-*        estimated distance to minimum, parameter errors, step sizes.
//*-*
//*-*         According to the value of IKODE, the printout is:/
//*-*    IKODE=INKODE= 0    only info about function value
//*-*                  1    parameter values, errors, limits
//*-*                  2    values, errors, step sizes, internal values
//*-*                  3    values, errors, step sizes, first derivs.
//*-*                  4    values, parabolic errors, MINOS errors
//*-*    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    static string cblank = "           ";
    static string cnambf = "           ";

    /* Local variables */
    Double_urt dcmax, x1, x2, x3, dc;
    x2 = x3 = 0;
    Int_urt nadd, i, k, l, m, ikode, ic, nc, ntrail, lbl;
    string chedm;
    string colhdl[6], colhdu[6], cx2, cx3, cheval;

    if (fNu == 0) {
	//std::printf(" THERE ARE CURRENTLY NO PARAMETERS DEFINED\n");
        *m_logStream << " THERE ARE CURRENTLY NO PARAMETERS DEFINED" << endl;
	return;
    }
//*-*-                 get value of IKODE based in INKODE, ISW(2)
    ikode = inkode;
    if (inkode == 5) {
	ikode = fISW[1] + 1;
	if (ikode > 3) ikode = 3;
    }
//*-*-                 set 'default' column headings
    for (k = 1; k <= 6; ++k) {
	colhdu[k-1] = "UNDEFINED";
	colhdl[k-1] = "COLUMN HEAD";
    }
//*-*-             print title if Minos errors, and title exists.
    if (ikode == 4 && fCtitl != fCundef) {
	//std::printf(" MINUIT TASK: %s\n",fCtitl.c_str());
        *m_logStream << " MINUIT TASK: " << fCtitl << '\n';
    }
//*-*-             report function value and status
    if (fval == fUndefi) cheval = " unknown       ";
   else {
      ostringstream och;
      och << fval;
      cheval = och.str();
//      = Form("%g",fval);
   }
    if (fEDM == fBigedm) chedm = " unknown  ";
   else  {
      ostringstream och;
      och << fEDM;
      chedm = och.str();
//       chedm = Form("%g",fEDM);
   }
    nc = fNfcn - fNfcnfr;
    //std::printf(" FCN=%s FROM %8s  STATUS=%10s %6d CALLS   %9d TOTAL\n"
    //           ,cheval.c_str()
    //           ,fCfrom.c_str()
    //           ,fCstatu.c_str(),nc,fNfcn);
    *m_logStream << " FCN=" << cheval << " FROM " << setw(8) << fCfrom << "   STATUS="
       << setw(10) << fCstatu << setw(7) << nc << "   " << setw(9) << fNfcn << " TOTAL\n" << endl;
    m = fISW[1];
    if (m == 0 || m == 2 || fDcovar == 0) {
	// std::printf("                     EDM=%s    STRATEGY=%2d      %s\n"
        //              ,chedm.c_str(),fIstrat
        //              ,fCovmes[m].c_str());
       *m_logStream << "                     EDM=" << chedm << "    STRATEGY=" << setw(2)
            << fIstrat << "      " << fCovmes[m] << '\n';
    } else {
	dcmax = 1;
	dc    = URMath::Min(fDcovar,dcmax)*100;
	//std::printf("                     EDM=%s    STRATEGY=%2d  ERROR MATRIX UNCERTAINTY %5.1f per cent\n"
        //              ,chedm.c_str(),fIstrat,dc);
	*m_logStream << "                     EDM=" << chedm
             << "    STRATEGY=" << setw(2) << fIstrat << "  ERROR MATRIX UNCERTAINTY "
             << setw(5) << setprecision(1) << dc << " per cent\n";
    }

    if (ikode == 0) return;
//*-*-              find longest name (for Rene!)
    ntrail = 10;
    for (i = 1; i <= fNu; ++i) {
///       if (fNvarl[i-1] < 0) continue;
	if (m_userParameterFlag[i] < 0) continue;
//        uint lastChar = fCpnam[i-1].size();
        uint lastChar = m_userParameterName[i].size();
        uint icMax = (lastChar < 10 ) ? lastChar : 10;
	for (ic = icMax; ic >= 1; --ic) {
//           if (fCpnam[i-1].substr(ic-1,1) != " ") goto L16;
	    if (m_userParameterName[i].substr(ic-1,1) != " ") goto L16;
	}
	ic = 1;
L16:
	lbl = 10 - ic;
	if (lbl < ntrail) ntrail = lbl;
    }
    nadd = ntrail / 2 + 1;
    if (ikode == 1) {
	colhdu[0] = "              ";
	colhdl[0] = "      ERROR   ";
	colhdu[1] = "      PHYSICAL";
	colhdu[2] = " LIMITS       ";
	colhdl[1] = "    NEGATIVE  ";
	colhdl[2] = "    POSITIVE  ";
    }
    if (ikode == 2) {
	colhdu[0] = "              ";
	colhdl[0] = "      ERROR   ";
	colhdu[1] = "    INTERNAL  ";
	colhdl[1] = "    STEP SIZE ";
	colhdu[2] = "    INTERNAL  ";
	colhdl[2] = "      VALUE   ";
    }
    if (ikode == 3) {
	colhdu[0] = "              ";
	colhdl[0] = "      ERROR   ";
	colhdu[1] = "       STEP   ";
	colhdl[1] = "       SIZE   ";
	colhdu[2] = "      FIRST   ";
	colhdl[2] = "   DERIVATIVE ";
    }
    if (ikode == 4) {
	colhdu[0] = "    PARABOLIC ";
	colhdl[0] = "      ERROR   ";
	colhdu[1] = "        MINOS ";
	colhdu[2] = "ERRORS        ";
	colhdl[1] = "   NEGATIVE   ";
	colhdl[2] = "   POSITIVE   ";
    }

    if (ikode != 4) {
	if (fISW[1] < 3) colhdu[0] = "  APPROXIMATE ";
	if (fISW[1] < 1) colhdu[0] = " CURRENT GUESS";
    }
    // get the maximum parameter name length
    int maxNameSize = 0;
    for ( vector<string>::iterator name = m_userParameterName.begin();
          name != m_userParameterName.end();
          ++name ) {
       int nameSize = name->size();
       maxNameSize = max(maxNameSize,nameSize);
    }
    maxNameSize = max(maxNameSize,9);
//    *m_logStream << "maxNameSize = " << maxNameSize << endl;
    
    *m_logStream << "  EXT " << setw(maxNameSize) << left << "PARAMETER" << "           "
         << setw(14) << colhdu[0].c_str() << setw(14) << colhdu[1].c_str() << setw(14) << colhdu[2].c_str() 
         << endl;
    *m_logStream << "  NO. " << setw(maxNameSize) << left << "  NAME   " << "  VALUE    "
         << setw(14) << colhdl[0].c_str() << setw(14) << colhdl[1].c_str() << setw(14) << colhdl[2].c_str() 
         << endl;
//    std::printf("  EXT PARAMETER              %-14s%-14s%-14s\n",colhdu[0].c_str()
//                                                    ,colhdu[1].c_str()
//                                                    ,colhdu[2].c_str());
//    std::printf("  NO.   NAME      VALUE      %-14s%-14s%-14s\n",colhdl[0].c_str()
//                                                    ,colhdl[1].c_str()
//                                                    ,colhdl[2].c_str());
//*-*-                                       . . . loop over parameters . .
    for (i = 1; i <= fNu; ++i) {
//       if (fNvarl[i-1] < 0) continue;
	if (m_userParameterFlag[i] < 0) continue;
//       l = fNiofex[i-1];
	l = m_userParameterIdToInternalId[i];
//        cnambf = cblank.substr(0,nadd) + fCpnam[i-1];
        cnambf = cblank.substr(0,nadd) + m_userParameterName[i];
	if (l == 0) goto L55;
//*-*-             variable parameter.
	x1  = fWerr[l-1];
	cx2 = "PLEASE GET X..";
	cx3 = "PLEASE GET X..";
	if (ikode == 1) {
//           if (fNvarl[i-1] <= 1) {
	    if (m_userParameterFlag[i] <= 1) {
//                std::printf("%4d %-11s%14.5e%14.5e\n",i,cnambf.c_str(),fU[i-1],x1);
//                std::printf("%4d %-11s%14.5e%14.5e\n",i,cnambf.c_str(),m_userParameterValue[i],x1);
                *m_logStream << setw(5) << right << i << " " << setw(maxNameSize) << m_userParameterName[i]
                     << setw(14) << setprecision(5) << m_userParameterValue[i]
                     << setw(14) << setprecision(5) << x1 << endl;
               continue;
	    } else {
		x2 = fAlim[i-1];
		x3 = fBlim[i-1];
	    }
	}
	if (ikode == 2) {
	    x2 = fDirin[l-1];
	    x3 = fX[l-1];
	}
	if (ikode == 3) {
	    x2 = fDirin[l-1];
	    x3 = fGrd[l-1];
//            if (fNvarl[i-1] > 1 && URMath::Abs(URMath::Cos(fX[l-1])) < .001) {
	    if (m_userParameterFlag[i] > 1 && URMath::Abs(URMath::Cos(fX[l-1])) < .001) {
		cx3 = "** at limit **";
	    }
	}
	if (ikode == 4) {
	    x2 = fErn[l-1];
	    if (x2 == 0)        cx2 = " ";
	    if (x2 == fUndefi)  cx2 = "   at limit   ";
	    x3 = fErp[l-1];
	    if (x3 == 0)        cx3 = " ";
	    if (x3 == fUndefi) 	cx3 = "   at limit   ";
	}
        if (cx2 == "PLEASE GET X..")  {
           ostringstream ocx;
           ocx << setw(14) << setprecision(5) << x2;
           cx2 = ocx.str();
//           cx2 = Form("%14.5e",x2);
        }
        if (cx3 == "PLEASE GET X..")  {
           ostringstream ocx;
           ocx << setw(14) << setprecision(5) << x3;
           cx3 = ocx.str();
           // cx3 = Form("%14.5e",x3);
        }
//        std::printf("%4d %-11s%14.5e%14.5e%-14s%-14s\n",i
//            ,cnambf.c_str(),fU[i-1],x1
//            ,cx2.c_str(),cx3.c_str());
//        std::printf("%4d %-11s%14.5e%14.5e%-14s%-14s\n",i
//                   ,cnambf.c_str(),m_userParameterValue[i],x1
//                   ,cx2.c_str(),cx3.c_str());
        *m_logStream << setw(5) << right << i << " " << setw(maxNameSize) << m_userParameterName[i]
             << setw(14) << setprecision(5) << m_userParameterValue[i]
             << setw(14) << setprecision(5) << x1 
             << cx2 << cx3 << endl;
        
//*-*-              check if parameter is at limit
//        if (fNvarl[i-1] <= 1 || ikode == 3) continue;
	if (m_userParameterFlag[i] <= 1 || ikode == 3) continue;
	if (URMath::Abs(URMath::Cos(fX[l-1])) < .001) {
	    //std::printf("                                 WARNING -   - ABOVE PARAMETER IS AT LIMIT.\n");
           *m_logStream << "                                 WARNING -   - ABOVE PARAMETER IS AT LIMIT.\n";
	}
	continue;

//*-*-                               print constant or fixed parameter.
L55:
	colhdu[0] = "   constant   ";
//        if (fNvarl[i-1] > 0)  colhdu[0] = "     fixed    ";
	if (m_userParameterFlag[i] > 0)  colhdu[0] = "     fixed    ";
//        if (fNvarl[i-1] == 4 && ikode == 1) {
	if (m_userParameterFlag[i] == 4 && ikode == 1) {
           //            std::printf("%4d %-11s%14.5e%-14s%14.5e%14.5e\n",i
           //              ,cnambf.c_str(),fU[i-1]
           //              ,colhdu[0].c_str(),fAlim[i-1],fBlim[i-1]);
           //std::printf("%4d %-11s%14.5e%-14s%14.5e%14.5e\n",i
           //   ,cnambf.c_str(),m_userParameterValue[i]
           //   ,colhdu[0].c_str(),fAlim[i-1],fBlim[i-1]);
           *m_logStream << setw(5) << right << i << " " << setw(maxNameSize) << m_userParameterName[i]
                << setw(14) << setprecision(5) << m_userParameterValue[i]
                << setw(14) << colhdu[0].c_str()
                << setw(14) << setprecision(5) << fAlim[i-1]
                << setw(14) << setprecision(5) << fBlim[i-1]
                << endl;
        } else {
//            std::printf("%4d %-11s%14.5e%s\n",i
//                   ,cnambf.c_str(),fU[i-1],colhdu[0].c_str());
//           std::printf("%4d %-11s%14.5e%s\n",i
//                   ,cnambf.c_str(),m_userParameterValue[i],colhdu[0].c_str());
           *m_logStream << setw(5) << right << i << " " << setw(maxNameSize) << m_userParameterName[i]
                << setw(14) << setprecision(5) << m_userParameterValue[i]
                << setw(14) << colhdu[0].c_str()
                << endl;
	}
    }

    if (fUp != fUpdflt) {
        std::printf("                               ERR DEF= %g\n",fUp);
       *m_logStream << "                               ERR DEF= " << fUp;
    }
    *m_logStream << endl;
    return;
} /* mnprin_ */

//______________________________________________________________________________
void URMinuit::mnpsdf()
{
//*-*-*-*-*-*Calculates the eigenvalues of v to see if positive-def*-*-*-*-*
//*-*        ======================================================
//*-*        if not, adds constant along diagonal to make positive.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dgmin, padd, pmin, pmax, dg, epspdf, epsmin;
    Int_urt ndex, i, j, ndexd, ip, ifault;
    string chbuff, ctemp;

    epsmin = 1e-6;
    epspdf = URMath::Max(epsmin,fEpsma2);
    dgmin  = fVhmat[0];
//*-*-                       Check if negative or zero on diagonal
    for (i = 1; i <= fNpar; ++i) {
	ndex = i*(i + 1) / 2;
	if (fVhmat[ndex-1] <= 0) {
           ostringstream warning;
           warning << "Negative diagonal element " << i << " in Error Matrix";
	    mnwarn("W", fCfrom.c_str(), warning.str().c_str());
	}
	if (fVhmat[ndex-1] < dgmin) dgmin = fVhmat[ndex-1];
    }
    if (dgmin <= 0) {
	dg    = epspdf + 1 - dgmin;
       ostringstream warning;
       warning << dg << " added to diagonal of error matrix";
	mnwarn("W", fCfrom.c_str(), warning.str().c_str());
    } else {
	dg = 0;
    }
//*-*-                   Store VHMAT in P, make sure diagonal pos.
    for (i = 1; i <= fNpar; ++i) {
	ndex  = i*(i-1) / 2;
	ndexd = ndex + i;
	fVhmat[ndexd-1] += dg;
	fPSDFs[i-1] = 1 / URMath::Sqrt(fVhmat[ndexd-1]);
	for (j = 1; j <= i; ++j) {
	    ++ndex;
	    fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[ndex-1]*fPSDFs[i-1]*fPSDFs[j-1];
	}
    }
//*-*-     call eigen (p,p,maxint,npar,pstar,-npar)
    mneig(fP, fMaxint, fNpar, fMaxint, fPstar, epspdf, ifault);
    pmin = fPstar[0];
    pmax = fPstar[0];
    for (ip = 2; ip <= fNpar; ++ip) {
	if (fPstar[ip-1] < pmin) pmin = fPstar[ip-1];
	if (fPstar[ip-1] > pmax) pmax = fPstar[ip-1];
    }
    pmax = URMath::Max(URMath::Abs(pmax),Double_urt(1));
    if ( ( pmin <= 0 && fLwarn ) || fISW[4] >= 2) {
	std::printf(" EIGENVALUES OF SECOND-DERIVATIVE MATRIX:\n");
        ctemp = "       ";
	for (ip = 1; ip <= fNpar; ++ip) {
           ostringstream localTemp;
           localTemp << fPstar[ip-1];
           ctemp += localTemp.str(); //Form(" %11.4e",);
	}
        std::printf("%s",ctemp.c_str());
    }
    if (pmin > epspdf*pmax) return;
    if (fISW[1] == 3) fISW[1] = 2;
    padd = pmax*.001 - pmin;
    for (ip = 1; ip <= fNpar; ++ip) {
	ndex = ip*(ip + 1) / 2;
	fVhmat[ndex-1] *= padd + 1;
    }
    fCstatu = "NOT POSDEF";
    ostringstream warning;
    warning << "MATRIX FORCED POS-DEF BY ADDING %f TO DIAGONAL. " << padd;
    mnwarn("W", fCfrom.c_str(), warning.str().c_str());

} /* mnpsdf_ */

//______________________________________________________________________________
void URMinuit::mnrazz(Double_urt ynew, Double_urt *pnew, Double_urt *y, Int_urt &jh, Int_urt &jl)
{
//*-*-*-*-*Called only by MNSIMP (and MNIMPR) to add a new point*-*-*-*-*-*-*
//*-*      =====================================================
//*-*        and remove an old one from the current simplex, and get the
//*-*        estimated distance to minimum.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt pbig, plit;
    Int_urt i, j, nparp1;

    /* Function Body */
    for (i = 1; i <= fNpar; ++i) { fP[i + jh*fMaxpar - fMaxpar-1] = pnew[i-1]; }
    y[jh-1] = ynew;
    if (ynew < fAmin) {
	for (i = 1; i <= fNpar; ++i) { fX[i-1] = pnew[i-1]; }
	mninex(fX);
	fAmin   = ynew;
	fCstatu = "PROGRESS  ";
	jl      = jh;
    }
    jh     = 1;
    nparp1 = fNpar + 1;
    for (j = 2; j <= nparp1; ++j) { if (y[j-1] > y[jh-1]) jh = j; }
    fEDM = y[jh-1] - y[jl-1];
    if (fEDM <= 0) goto L45;
    for (i = 1; i <= fNpar; ++i) {
	pbig = fP[i-1];
	plit = pbig;
	for (j = 2; j <= nparp1; ++j) {
	    if (fP[i + j*fMaxpar - fMaxpar-1] > pbig) pbig = fP[i + j*fMaxpar - fMaxpar-1];
	    if (fP[i + j*fMaxpar - fMaxpar-1] < plit) plit = fP[i + j*fMaxpar - fMaxpar-1];
	}
	fDirin[i-1] = pbig - plit;
    }
L40:
    return;
L45:
    std::printf("   FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE%d VARIABLE PARAMETERS.\n",fNpar);
    std::printf("          VERIFY THAT STEP SIZES ARE BIG ENOUGH AND CHECK FCN LOGIC.\n");
    std::printf(" *******************************************************************************\n");
    std::printf(" *******************************************************************************\n");
    goto L40;
} /* mnrazz_ */

//______________________________________________________________________________
void URMinuit::mnrn15(Double_urt &val, Int_urt &inseed)
{
//*-*-*-*-*-*-*This is a super-portable random number generator*-*-*-*-*-*-*
//*-*          ================================================
//*-*         It should not overflow on any 32-bit machine.
//*-*         The cycle is only ~10**9, so use with care!
//*-*         Note especially that VAL must not be undefined on input.
//*-*                    Set Default Starting Seed
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    static Int_urt iseed = 12345;

    Int_urt k;

    if (val == 3) goto L100;
    inseed = iseed;
    k      = iseed / 53668;
    iseed  = (iseed - k*53668)*40014 - k*12211;
    if (iseed < 0) iseed += 2147483563;
    val = Double_urt(iseed*4.656613e-10);
    return;
//*-*               "entry" to set seed, flag is VAL=3
L100:
    iseed = inseed;
} /* mnrn15_ */

//______________________________________________________________________________
void URMinuit::mnrset(Int_urt iopt)
{
//*-*-*-*-*-*-*-*Resets function value and errors to UNDEFINED*-*-*-*-*-*-*-*
//*-*            =============================================
//*-*    If IOPT=1,
//*-*    If IOPT=0, sets only MINOS errors to undefined
//*-*        Called from MNCLER and whenever problem changes, for example
//*-*        after SET LIMITS, SET PARAM, CALL FCN 6
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt iext, i;

    fCstatu = "RESET     ";
    if (iopt >= 1) {
	fAmin   = fUndefi;
	fFval3  = URMath::Abs(fAmin)*2 + 1;
	fEDM    = fBigedm;
	fISW[3] = 0;
	fISW[1] = 0;
	fDcovar = 1;
	fISW[0] = 0;
    }
    fLnolim = kurTRUE;
    for (i = 1; i <= fNpar; ++i) {
	iext = fNexofi[i-1];
//       if (fNvarl[iext-1] >= 4) fLnolim = kurFALSE;
	if (m_userParameterFlag[iext] >= 4) fLnolim = kurFALSE;
	fErp[i-1] = 0;
	fErn[i-1] = 0;
	fGlobcc[i-1] = 0;
    }
    if (fISW[1] >= 1) {
	fISW[1] = 1;
	fDcovar = URMath::Max(fDcovar,.5);
    }
} /* mnrset_ */

//______________________________________________________________________________
void URMinuit::mnsave()
{
//*-*-*-*Writes current parameter values and step sizes onto file ISYSSA*-*-*
//*-*    ===============================================================
//*-*          in format which can be reread by Minuit for restarting.
//*-*       The covariance matrix is also output if it exists.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

   std::printf("mnsave is dummy in the base class URMinuit: Use URMinuitOld\n");

} /* mnsave_ */

//______________________________________________________________________________
void URMinuit::mnscan()
{
//*-*-*-*-*Scans the values of FCN as a function of one parameter*-*-*-*-*-*
//*-*      ======================================================
//*-*        and plots the resulting values as a curve using MNPLOT.
//*-*        It may be called to scan one parameter or all parameters.
//*-*        retains the best function and parameter values found.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt step, uhigh, xhreq, xlreq, ubest, fnext, unext, xh, xl;
    Int_urt ipar, iint, icall, ncall, nbins, nparx;
    Int_urt nxypt, nccall, iparwd;

    xlreq = URMath::Min(fWord7[2],fWord7[3]);
    xhreq = URMath::Max(fWord7[2],fWord7[3]);
    ncall = Int_urt((fWord7[1] + .01));
    if (ncall <= 1)  ncall = 41;
    if (ncall > 101) ncall = 101;
    nccall = ncall;
    if (fAmin == fUndefi) mnamin();
    iparwd  = Int_urt((fWord7[0] + .1));
    ipar    = URMath::Max(iparwd,0);
//    iint    = fNiofex[ipar-1];
    iint    = m_userParameterIdToInternalId[ipar];
    fCstatu = "NO CHANGE";
    if (iparwd > 0) goto L200;

//*-*-        equivalent to a loop over parameters requested
L100:
    ++ipar;
    if (ipar > fNu) goto L900;
//    iint = fNiofex[ipar-1];
    iint = m_userParameterIdToInternalId[ipar];
    if (iint <= 0) goto L100;
//*-*-        set up range for parameter IPAR
L200:
//    ubest    = fU[ipar-1];
    ubest    = m_userParameterValue[ipar];
    fXpt[0]  = ubest;
    fYpt[0]  = fAmin;
    fChpt[0] = ' ';
    fXpt[1]  = ubest;
    fYpt[1]  = fAmin;
    fChpt[1] = 'X';
    nxypt    = 2;
//    if (fNvarl[ipar-1] > 1) goto L300;
    if (m_userParameterFlag[ipar] > 1) goto L300;

//*-*-        no limits on parameter
    if (xlreq == xhreq) goto L250;
    unext = xlreq;
    step = (xhreq - xlreq) / Double_urt(ncall-1);
    goto L500;
L250:
    xl = ubest - fWerr[iint-1];
    xh = ubest + fWerr[iint-1];
    mnbins(xl, xh, ncall, unext, uhigh, nbins, step);
    nccall = nbins + 1;
    goto L500;
//*-*-        limits on parameter
L300:
    if (xlreq == xhreq) goto L350;
//*-*  Computing MAX
    xl = URMath::Max(xlreq,fAlim[ipar-1]);
//*-*  Computing MIN
    xh = URMath::Min(xhreq,fBlim[ipar-1]);
    if (xl >= xh) goto L700;
    unext = xl;
    step  = (xh - xl) / Double_urt(ncall-1);
    goto L500;
L350:
    unext = fAlim[ipar-1];
    step = (fBlim[ipar-1] - fAlim[ipar-1]) / Double_urt(ncall-1);
//*-*-        main scanning loop over parameter IPAR
L500:
    for (icall = 1; icall <= nccall; ++icall) {
//       fU[ipar-1] = unext;
	m_userParameterValue[ipar] = unext;
	nparx = fNpar;
//        Eval(nparx, fGin, fnext, fU, 4);	++fNfcn;
	Eval(nparx, fGin, fnext, m_userParameterValue, 4);	++fNfcn;
	++nxypt;
	fXpt[nxypt-1]  = unext;
	fYpt[nxypt-1]  = fnext;
	fChpt[nxypt-1] = '*';
	if (fnext < fAmin) {
	    fAmin   = fnext;
	    ubest   = unext;
	    fCstatu = "IMPROVED  ";
	}
	unext += step;
    }
    fChpt[nccall] = 0;

//*-*-        finished with scan of parameter IPAR
//    fU[ipar-1] = ubest;
    m_userParameterValue[ipar] = ubest;
    mnexin(fX);
//    std::printf("%dSCAN OF PARAMETER NO. %d,  %s\n"
//         ,fNewpag,ipar,fCpnam[ipar-1].c_str());
    std::printf("%dSCAN OF PARAMETER NO. %d,  %s\n"
            ,fNewpag,ipar,m_userParameterName[ipar].c_str());
    mnplot(fXpt, fYpt, fChpt, nxypt, fNpagwd, fNpagln);
    goto L800;
L700:
    std::printf(" REQUESTED RANGE OUTSIDE LIMITS FOR PARAMETER  %d\n",ipar);
L800:
    if (iparwd <= 0) goto L100;
//*-*-        finished with all parameters
L900:
    mnprin(5, fAmin);
} /* mnscan_ */

//______________________________________________________________________________
void URMinuit::mnseek()
{
//*-*-*-*Performs a rough (but global) minimization by monte carlo search*-*
//*-*    ================================================================
//*-*        Each time a new minimum is found, the search area is shifted
//*-*        to be centered at the best value.  Random points are chosen
//*-*        uniformly over a hypercube determined by current step sizes.
//*-*   The Metropolis algorithm accepts a worse point with probability
//*-*      exp(-d/UP), where d is the degradation.  Improved points
//*-*      are of course always accepted.  Actual steps are random
//*-*      multiples of the nominal steps (DIRIN).
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Local variables */
    Double_urt dxdi, rnum, ftry, rnum1, rnum2, alpha;
    Double_urt flast, bar;
    Int_urt ipar, iext, j, ifail, iseed, nparx, istep, ib, mxfail, mxstep;

    mxfail = Int_urt(fWord7[0]);
    if (mxfail <= 0) mxfail = fNpar*20 + 100;
    mxstep = mxfail*10;
    if (fAmin == fUndefi) mnamin();
    alpha = fWord7[1];
    if (alpha <= 0) alpha = 3;
    if (fISW[4] >= 1) {
	std::printf(" MNSEEK: MONTE CARLO MINIMIZATION USING METROPOLIS ALGORITHM\n");
        std::printf(" TO STOP AFTER %6d SUCCESSIVE FAILURES, OR %7d STEPS\n",mxfail,mxstep);
        std::printf(" MAXIMUM STEP SIZE IS %9.3f ERROR BARS.\n",alpha);
    }
    fCstatu = "INITIAL  ";
    if (fISW[4] >= 2) mnprin(2, fAmin);
    fCstatu = "UNCHANGED ";
    ifail   = 0;
    rnum    = 0;
    rnum1   = 0;
    rnum2   = 0;
    nparx   = fNpar;
    flast   = fAmin;
//*-*-             set up step sizes, starting values
    for (ipar = 1; ipar <= fNpar; ++ipar) {
	iext = fNexofi[ipar-1];
	fDirin[ipar-1] = alpha*2*fWerr[ipar-1];
//        if (fNvarl[iext-1] > 1) {
	if (m_userParameterFlag[iext] > 1) {
//*-*-             parameter with limits
	    mndxdi(fX[ipar-1], ipar-1, dxdi);
	    if (dxdi == 0) dxdi = 1;
	    fDirin[ipar-1] = alpha*2*fWerr[ipar-1] / dxdi;
	    if (URMath::Abs(fDirin[ipar-1]) > 6.2831859999999997) {
		fDirin[ipar-1] = 6.2831859999999997;
	    }
	}
	fSEEKxmid[ipar-1]  = fX[ipar-1];
	fSEEKxbest[ipar-1] = fX[ipar-1];
    }
//*-*-                             search loop
    for (istep = 1; istep <= mxstep; ++istep) {
	if (ifail >= mxfail) break;
	for (ipar = 1; ipar <= fNpar; ++ipar) {
	    mnrn15(rnum1, iseed);
	    mnrn15(rnum2, iseed);
	    fX[ipar-1] = fSEEKxmid[ipar-1] + (rnum1 + rnum2 - 1)*.5*fDirin[ipar-1];
	}
	mninex(fX);
//        Eval(nparx, fGin, ftry, fU, 4);	++fNfcn;
	Eval(nparx, fGin, ftry, m_userParameterValue, 4);	++fNfcn;
	if (ftry < flast) {
	    if (ftry < fAmin) {
		fCstatu = "IMPROVEMNT";
		fAmin = ftry;
		for (ib = 1; ib <= fNpar; ++ib) { fSEEKxbest[ib-1] = fX[ib-1]; }
		ifail = 0;
		if (fISW[4] >= 2) mnprin(2, fAmin);
	    }
	    goto L300;
	} else {
	    ++ifail;
//*-*-                  Metropolis algorithm
	    bar = (fAmin - ftry) / fUp;
	    mnrn15(rnum, iseed);
	    if (bar < URMath::Log(rnum)) continue;
	}
//*-*-                   Accept new point, move there
L300:
	for (j = 1; j <= fNpar; ++j) { fSEEKxmid[j-1] = fX[j-1];	}
	flast = ftry;
    }
//*-*-                              end search loop
    if (fISW[4] > 1) {
        std::printf(" MNSEEK: %5d SUCCESSIVE UNSUCCESSFUL TRIALS.\n",ifail);
    }
    for (ib = 1; ib <= fNpar; ++ib) { fX[ib-1] = fSEEKxbest[ib-1]; }
    mninex(fX);
    if (fISW[4] >= 1) mnprin(2, fAmin);
    if (fISW[4] == 0) mnprin(0, fAmin);
} /* mnseek_ */

//______________________________________________________________________________
void URMinuit::mnset()
{
//*-*-*-*-*Interprets the commands that start with SET and SHOW*-*-*-*-*-*-*
//*-*      ====================================================
//*-*        Called from MNEXCM
//*-*        file characteristics for SET INPUT
//*-*       'SET ' or 'SHOW',  'ON ' or 'OFF', 'SUPPRESSED' or 'REPORTED  '
//*-*        explanation of print level numbers -1:3  and strategies 0:2
//*-*        identification of debug options
//*-*        things that can be set or shown
//*-*        options not intended for normal users
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    const char *cname[30] = {
      "FCN value ",
      "PARameters",
      "LIMits    ",
      "COVariance",
      "CORrelatio",
      "PRInt levl",
      "NOGradient",
      "GRAdient  ",
      "ERRor def ",
      "INPut file",
      "WIDth page",
      "LINes page",
      "NOWarnings",
      "WARnings  ",
      "RANdom gen",
      "TITle     ",
      "STRategy  ",
      "EIGenvalue",
      "PAGe throw",
      "MINos errs",
      "EPSmachine",
      "OUTputfile",
      "BATch     ",
      "INTeractiv",
      "VERsion   ",
      "reserve   ",
      "NODebug   ",
      "DEBug     ",
      "SHOw      ",
      "SET       "};

    static Int_urt nname = 25;
    static Int_urt nntot = 30;
    static string cprlev[5] = {
      "-1: NO OUTPUT EXCEPT FROM SHOW    ",
      " 0: REDUCED OUTPUT                ",
      " 1: NORMAL OUTPUT                 ",
      " 2: EXTRA OUTPUT FOR PROBLEM CASES",
      " 3: MAXIMUM OUTPUT                "};

    static string cstrat[3] = {
      " 0: MINIMIZE THE NUMBER OF CALLS TO FUNCTION",
      " 1: TRY TO BALANCE SPEED AGAINST RELIABILITY",
      " 2: MAKE SURE MINIMUM TRUE, ERRORS CORRECT  "};

    static string cdbopt[7] = {
      "REPORT ALL EXCEPTIONAL CONDITIONS      ",
      "MNLINE: LINE SEARCH MINIMIZATION       ",
      "MNDERI: FIRST DERIVATIVE CALCULATIONS  ",
      "MNHESS: SECOND DERIVATIVE CALCULATIONS ",
      "MNMIGR: COVARIANCE MATRIX UPDATES      ",
      "MNHES1: FIRST DERIVATIVE UNCERTAINTIES ",
      "MNCONT: MNCONTOUR PLOT (MNCROS SEARCH) "};

    /* System generated locals */
//    Int_urt f_inqu();

    /* Local variables */
    Double_urt val;
    Int_urt iset, iprm, i, jseed, kname, iseed, iunit, id, ii, kk;
    Int_urt ikseed, idbopt, igrain, iswsav, isw2;
    string  cfname, cmode, ckind,  cwarn, copt, ctemp, ctemp2;
    Bool_urt lname=kurFALSE;

    for (i = 1; i <= nntot; ++i) {
	ctemp  = cname[i-1];
	ckind  = ctemp.substr(0,3);
        ctemp2 = fCword.substr(4,6);
	if (strstr(ctemp2.c_str(),ckind.c_str())) goto L5;
    }
    i = 0;
L5:
    kname = i;

//*-*-          Command could be SET xxx, SHOW xxx,  HELP SET or HELP SHOW
    ctemp2 = fCword.substr(0,3);
    if ( ctemp2.find("HEL") != string::npos)  goto L2000;
    if ( ctemp2.find("SHO") != string::npos)  goto L1000;
    if ( ctemp2.find("SET") == string::npos)  goto L1900;
//*-*-                          ---
    ckind = "SET ";
//*-*-                                       . . . . . . . . . . set unknown
    if (kname <= 0) goto L1900;
//*-*-                                       . . . . . . . . . . set known
    switch ((int)kname) {
	case 1:  goto L3000;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
	case 5:  goto L3000;
	case 6:  goto L60;
	case 7:  goto L70;
	case 8:  goto L80;
	case 9:  goto L90;
	case 10:  goto L100;
	case 11:  goto L110;
	case 12:  goto L120;
	case 13:  goto L130;
	case 14:  goto L140;
	case 15:  goto L150;
	case 16:  goto L160;
	case 17:  goto L170;
	case 18:  goto L3000;
	case 19:  goto L190;
	case 20:  goto L3000;
	case 21:  goto L210;
	case 22:  goto L220;
	case 23:  goto L230;
	case 24:  goto L240;
	case 25:  goto L3000;
	case 26:  goto L1900;
	case 27:  goto L270;
	case 28:  goto L280;
	case 29:  goto L290;
	case 30:  goto L300;
    }

//*-*-                                       . . . . . . . . . . set param
L20:
    iprm = Int_urt(fWord7[0]);
    if (iprm > fNu) goto L25;
    if (iprm <= 0) goto L25;
//    if (fNvarl[iprm-1] < 0) goto L25;
    if (m_userParameterFlag[iprm] < 0) goto L25;
//    fU[iprm-1] = fWord7[1];
    m_userParameterValue[iprm] = fWord7[1];
    mnexin(fX);
    isw2 = fISW[1];
    mnrset(1);
//*-*-       Keep approximate covariance matrix, even if new param value
    fISW[1] = URMath::Min(isw2,1);
    fCfrom  = "SET PARM";
    fNfcnfr = fNfcn;
    fCstatu = "NEW VALUES";
    return;
L25:
    std::printf(" UNDEFINED PARAMETER NUMBER.  IGNORED.\n");
    return;
//*-*-                                       . . . . . . . . . . set limits
L30:
    mnlims();
    return;
//*-*-                                       . . . . . . . . . . set covar
L40:
//*-*   this command must be handled by MNREAD, and is not Fortran-callable
    goto L3000;
//*-*-                                       . . . . . . . . . . set print
L60:
    fISW[4] = Int_urt(fWord7[0]);
    return;
//*-*-                                       . . . . . . . . . . set nograd
L70:
    fISW[2] = 0;
    return;
//*-*-                                       . . . . . . . . . . set grad
L80:
    mngrad();
    return;
//*-*-                                       . . . . . . . . . . set errdef
L90:
    if (fWord7[0] == fUp) return;
    if (fWord7[0] <= 0) {
	if (fUp == fUpdflt) return;
	fUp = fUpdflt;
    } else {
	fUp = fWord7[0];
    }
    for (i = 1; i <= fNpar; ++i) {
	fErn[i-1] = 0;
	fErp[i-1] = 0;
    }
    mnwerr();
    return;
//*-*-                                       . . . . . . . . . . set input
//*-* This command must be handled by MNREAD. If it gets this far,
//*-*-        it is illegal.
L100:
    goto L3000;
//*-*-                                       . . . . . . . . . . set width
L110:
    fNpagwd = Int_urt(fWord7[0]);
    fNpagwd = URMath::Max(fNpagwd,50);
    return;

L120:
    fNpagln = Int_urt(fWord7[0]);
    return;
//*-*-                                       . . . . . . . . . . set nowarn

L130:
    fLwarn = kurFALSE;
    return;
//*-*-                                       . . . . . . . . . . set warn
L140:
    fLwarn = kurTRUE;
    mnwarn("W", "SHO", "SHO");
    return;
//*-*-                                       . . . . . . . . . . set random
L150:
    jseed = Int_urt(fWord7[0]);
    val = 3;
    mnrn15(val, jseed);
    if (fISW[4] > 0) {
        std::printf(" MINUIT RANDOM NUMBER SEED SET TO %d\n",jseed);
    }
    return;
//*-*-                                       . . . . . . . . . . set title
L160:
//*-*   this command must be handled by MNREAD, and is not Fortran-callable
    goto L3000;
//*-*-                                       . . . . . . . . . set strategy
L170:
    fIstrat = Int_urt(fWord7[0]);
    fIstrat = URMath::Max(fIstrat,0);
    fIstrat = URMath::Min(fIstrat,2);
    if (fISW[4] > 0) goto L1172;
    return;
//*-*-                                      . . . . . . . . . set page throw
L190:
    fNewpag = Int_urt(fWord7[0]);
    goto L1190;
//*-*-                                       . . . . . . . . . . set epsmac
L210:
    if (fWord7[0] > 0 && fWord7[0] < .1) {
	fEpsmac = fWord7[0];
    }
    fEpsma2 = URMath::Sqrt(fEpsmac);
    goto L1210;
//*-*-                                      . . . . . . . . . . set outputfile
L220:
    iunit = Int_urt(fWord7[0]);
    fIsyswr = iunit;
    fIstkwr[0] = iunit;
    if (fISW[4] >= 0) goto L1220;
    return;
//*-*-                                       . . . . . . . . . . set batch
L230:
    fISW[5] = 0;
    if (fISW[4] >= 0) goto L1100;
    return;
//*-*-                                      . . . . . . . . . . set interactive
L240:
    fISW[5] = 1;
    if (fISW[4] >= 0) goto L1100;
    return;
//*-*-                                       . . . . . . . . . . set nodebug
L270:
    iset = 0;
    goto L281;
//*-*-                                       . . . . . . . . . . set debug
L280:
    iset = 1;
L281:
    idbopt = Int_urt(fWord7[0]);
    if (idbopt > 6) goto L288;
    if (idbopt >= 0) {
	fIdbg[idbopt] = iset;
	if (iset == 1) fIdbg[0] = 1;
    } else {
//*-*-            SET DEBUG -1  sets all debug options
	for (id = 0; id <= 6; ++id) { fIdbg[id] = iset; }
    }
    fLrepor = fIdbg[0] >= 1;
    mnwarn("D", "SHO", "SHO");
    return;
L288:
    std::printf(" UNKNOWN DEBUG OPTION %d REQUESTED. IGNORED\n",idbopt);
    return;
//*-*-                                       . . . . . . . . . . set show
L290:
//*-*-                                       . . . . . . . . . . set set
L300:
    goto L3000;
//*-*-               -----------------------------------------------------
L1000:
//*-*-              at this point, CWORD must be 'SHOW'
    ckind = "SHOW";
    if (kname <= 0) goto L1900;

    switch ((int)kname) {
	case 1:  goto L1010;
	case 2:  goto L1020;
	case 3:  goto L1030;
	case 4:  goto L1040;
	case 5:  goto L1050;
	case 6:  goto L1060;
	case 7:  goto L1070;
	case 8:  goto L1070;
	case 9:  goto L1090;
	case 10:  goto L1100;
	case 11:  goto L1110;
	case 12:  goto L1120;
	case 13:  goto L1130;
	case 14:  goto L1130;
	case 15:  goto L1150;
	case 16:  goto L1160;
	case 17:  goto L1170;
	case 18:  goto L1180;
	case 19:  goto L1190;
	case 20:  goto L1200;
	case 21:  goto L1210;
	case 22:  goto L1220;
	case 23:  goto L1100;
	case 24:  goto L1100;
	case 25:  goto L1250;
	case 26:  goto L1900;
	case 27:  goto L1270;
	case 28:  goto L1270;
	case 29:  goto L1290;
	case 30:  goto L1300;
    }

//*-*-                                       . . . . . . . . . . show fcn
L1010:
    if (fAmin == fUndefi) mnamin();
    mnprin(0, fAmin);
    return;
//*-*-                                       . . . . . . . . . . show param
L1020:
    if (fAmin == fUndefi) mnamin();
    mnprin(5, fAmin);
    return;
//*-*-                                       . . . . . . . . . . show limits
L1030:
    if (fAmin == fUndefi) mnamin();
    mnprin(1, fAmin);
    return;
//*-*-                                       . . . . . . . . . . show covar
L1040:
    mnmatu(1);
    return;
//*-*-                                       . . . . . . . . . . show corre
L1050:
    mnmatu(0);
    return;
//*-*-                                       . . . . . . . . . . show print
L1060:
    if (fISW[4] < -1) fISW[4] = -1;
    if (fISW[4] > 3)  fISW[4] = 3;
    std::printf(" ALLOWED PRINT LEVELS ARE:\n");
    std::printf("                           %s\n",cprlev[0].c_str());
    std::printf("                           %s\n",cprlev[1].c_str());
    std::printf("                           %s\n",cprlev[2].c_str());
    std::printf("                           %s\n",cprlev[3].c_str());
    std::printf("                           %s\n",cprlev[4].c_str());
    std::printf(" CURRENT PRINTOUT LEVEL IS %s\n",cprlev[fISW[4]+1].c_str());
    return;
//*-*-                                       . . . . . . . show nograd, grad
L1070:
    if (fISW[2] <= 0) {
	std::printf(" NOGRAD IS SET.  DERIVATIVES NOT COMPUTED IN FCN.\n");
    } else {
	std::printf("   GRAD IS SET.  USER COMPUTES DERIVATIVES IN FCN.\n");
    }
    return;
//*-*-                                      . . . . . . . . . . show errdef
L1090:
    std::printf(" ERRORS CORRESPOND TO FUNCTION CHANGE OF %g\n",fUp);
    return;
//*-*-                                      . . . . . . . . . . show input,
//*-*-                                               batch, or interactive
L1100:
//    ioin__1.inerr = 0;
//    ioin__1.inunit = fIsysrd;
//    ioin__1.infile = 0;
//    ioin__1.inex = 0;
//    ioin__1.inopen = 0;
//    ioin__1.innum = 0;
//    ioin__1.innamed = &lname;
//    ioin__1.innamlen = 64;
//    ioin__1.inname = cfname;
//    ioin__1.inacc = 0;
//    ioin__1.inseq = 0;
//    ioin__1.indir = 0;
//    ioin__1.infmt = 0;
//    ioin__1.inform = 0;
//    ioin__1.inunf = 0;
//    ioin__1.inrecl = 0;
//    ioin__1.innrec = 0;
//    ioin__1.inblank = 0;
//    f_inqu(&ioin__1);
    cmode = "BATCH MODE      ";
    if (fISW[5] == 1) cmode  = "INTERACTIVE MODE";
    if (! lname)      cfname = "unknown";
    std::printf(" INPUT NOW BEING READ IN %s FROM UNIT NO. %d FILENAME: %s\n"
           ,cmode.c_str(),fIsysrd,cfname.c_str());
    return;
//*-*-                                      . . . . . . . . . . show width
L1110:
    std::printf("          PAGE WIDTH IS SET TO %d COLUMNS\n",fNpagwd);
    return;
//*-*-                                      . . . . . . . . . . show lines
L1120:
    std::printf("          PAGE LENGTH IS SET TO %d LINES\n",fNpagln);
    return;
//*-*-                                      . . . . . . .show nowarn, warn
L1130:
    cwarn = "SUPPRESSED";
    if (fLwarn) cwarn = "REPORTED  ";
    std::printf("%s\n",cwarn.c_str());
    if (! fLwarn) mnwarn("W", "SHO", "SHO");
    return;
//*-*-                                     . . . . . . . . . . show random
L1150:
    val = 0;
    mnrn15(val, igrain);
    ikseed = igrain;
    std::printf(" MINUIT RNDM SEED IS CURRENTLY=%d\n",ikseed);
    val   = 3;
    iseed = ikseed;
    mnrn15(val, iseed);
    return;
//*-*-                                       . . . . . . . . . show title
L1160:
    std::printf(" TITLE OF CURRENT TASK IS:%s\n",fCtitl.c_str());
    return;
//*-*-                                       . . . . . . . show strategy
L1170:
    std::printf(" ALLOWED STRATEGIES ARE:\n");
    std::printf("                    %s\n",cstrat[0].c_str());
    std::printf("                    %s\n",cstrat[1].c_str());
    std::printf("                    %s\n",cstrat[2].c_str());
L1172:
    std::printf(" NOW USING STRATEGY %s\n",cstrat[fIstrat].c_str());
    return;
//*-*-                                         . . . . . show eigenvalues
L1180:
    iswsav = fISW[4];
    fISW[4] = 3;
    if (fISW[1] < 1) {
	std::printf("%s\n",fCovmes[0].c_str());
    } else {
	mnpsdf();
    }
    fISW[4] = iswsav;
    return;
//*-*-                                           . . . . . show page throw
L1190:
    std::printf(" PAGE THROW CARRIAGE CONTROL = %d\n",fNewpag);
    if (fNewpag == 0) {
	std::printf(" NO PAGE THROWS IN MINUIT OUTPUT\n");
    }
    return;
//*-*-                                       . . . . . . show minos errors
L1200:
    for (ii = 1; ii <= fNpar; ++ii) {
	if (fErp[ii-1] > 0 || fErn[ii-1] < 0) goto L1204;
    }
    std::printf("       THERE ARE NO MINOS ERRORS CURRENTLY VALID.\n");
    return;
L1204:
    mnprin(4, fAmin);
    return;
//*-*-                                       . . . . . . . . . show epsmac
L1210:
    std::printf(" FLOATING-POINT NUMBERS ASSUMED ACCURATE TO %g\n",fEpsmac);
    return;
//*-*-                                       . . . . . . show outputfiles
L1220:
    std::printf("  MINUIT PRIMARY OUTPUT TO UNIT %d\n",fIsyswr);
    return;
//*-*-                                       . . . . . . show version
L1250:
    std::printf(" THIS IS MINUIT VERSION:%s\n",fCvrsn.c_str());
    return;
//*-*-                                       . . . . . . show nodebug, debug
L1270:
    for (id = 0; id <= 6; ++id) {
	copt = "OFF";
	if (fIdbg[id] >= 1) copt = "ON ";
        std::printf("          DEBUG OPTION %3d IS %3s :%s\n"
             ,id,copt.c_str(),cdbopt[id].c_str());
    }
    if (! fLrepor) mnwarn("D", "SHO", "SHO");
    return;
//*-*-                                       . . . . . . . . . . show show
L1290:
    ckind = "SHOW";
    goto L2100;
//*-*-                                       . . . . . . . . . . show set
L1300:
    ckind = "SET ";
    goto L2100;
//*-*-               -----------------------------------------------------
//*-*-                             UNKNOWN COMMAND
L1900:
    std::printf(" THE COMMAND:%10s IS UNKNOWN.\n",fCword.c_str());
    goto L2100;
//*-*-               -----------------------------------------------------
//*-*-                   HELP SHOW,  HELP SET,  SHOW SET, or SHOW SHOW
L2000:
    ckind = "SET ";
    ctemp2 = fCword.substr(3,7);
    if (strcmp(ctemp2.c_str(), "SHO")) ckind = "SHOW";
L2100:
    std::printf(" THE FORMAT OF THE %4s COMMAND IS:\n",ckind.c_str());
    std::printf(" %s xxx    [numerical arguments if any]\n",ckind.c_str());
    std::printf(" WHERE xxx MAY BE ONE OF THE FOLLOWING:\n");
    for (kk = 1; kk <= nname; ++kk) {
        std::printf(" %s\n",cname[kk-1]);
    }
    return;
//*-*-               -----------------------------------------------------
//*-*-                              ILLEGAL COMMAND
L3000:
    std::printf(" ABOVE COMMAND IS ILLEGAL.   IGNORED\n");

} /* mnset_ */

//______________________________________________________________________________
void URMinuit::mnsimp()
{
//*-*-*-*-*Minimization using the simplex method of Nelder and Mead*-*-*-*-*
//*-*      ========================================================
//*-*        Performs a minimization using the simplex method of Nelder
//*-*        and Mead (ref. -- Comp. J. 7,308 (1965)).
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* Initialized data */

    static Double_urt alpha = 1;
    static Double_urt beta = .5;
    static Double_urt gamma = 2;
    static Double_urt rhomin = 4;
    static Double_urt rhomax = 8;

    /* Local variables */
    Double_urt dmin_, dxdi, yrho, f, ynpp1, aming, ypbar;
    Double_urt bestx, ystar, y1, y2, ystst, pb, wg;
    Double_urt absmin, rho, sig2, rho1, rho2;
    Int_urt npfn, i, j, k, jhold, ncycl, nparx;
    Int_urt nparp1, kg, jh, nf, jl, ns;

    if (fNpar <= 0) return;
    if (fAmin == fUndefi) mnamin();
    fCfrom  = "SIMPLEX ";
    fNfcnfr = fNfcn;
    fCstatu = "UNCHANGED ";
    npfn    = fNfcn;
    nparp1  = fNpar + 1;
    nparx   = fNpar;
    rho1    = alpha + 1;
    rho2    = rho1 + alpha*gamma;
    wg      = 1 / Double_urt(fNpar);
    if (fISW[4] >= 0) {
	std::printf(" START SIMPLEX MINIMIZATION.    CONVERGENCE WHEN EDM .LT. %g\n",fEpsi);
    }
    for (i = 1; i <= fNpar; ++i) {
	fDirin[i-1] = fWerr[i-1];
	mndxdi(fX[i-1], i-1, dxdi);
	if (dxdi != 0) fDirin[i-1] = fWerr[i-1] / dxdi;
	dmin_ = fEpsma2*URMath::Abs(fX[i-1]);
	if (fDirin[i-1] < dmin_) fDirin[i-1] = dmin_;
    }
//*-* **       choose the initial simplex using single-parameter searches
L1:
    ynpp1 = fAmin;
    jl = nparp1;
    fSIMPy[nparp1-1] = fAmin;
    absmin = fAmin;
    for (i = 1; i <= fNpar; ++i) {
	aming      = fAmin;
	fPbar[i-1] = fX[i-1];
	bestx      = fX[i-1];
	kg         = 0;
	ns         = 0;
	nf         = 0;
L4:
	fX[i-1] = bestx + fDirin[i-1];
	mninex(fX);
//	Eval(nparx, fGin, f, fU, 4);	++fNfcn;
        Eval(nparx, fGin, f, m_userParameterValue, 4);	++fNfcn;
	if (f <= aming) goto L6;
//*-*-        failure
	if (kg == 1) goto L8;
	kg = -1;
	++nf;
	fDirin[i-1] *= -.4;
	if (nf < 3) goto L4;
	ns = 6;
//*-*-        success
L6:
	bestx        = fX[i-1];
	fDirin[i-1] *= 3;
	aming        = f;
	fCstatu      = "PROGRESS  ";
	kg           = 1;
	++ns;
	if (ns < 6) goto L4;
//*-*-        local minimum found in ith direction
L8:
	fSIMPy[i-1] = aming;
	if (aming < absmin) jl = i;
	if (aming < absmin) absmin = aming;
	fX[i-1] = bestx;
	for (k = 1; k <= fNpar; ++k) { fP[k + i*fMaxpar - fMaxpar-1] = fX[k-1]; }
    }
    jh    = nparp1;
    fAmin = fSIMPy[jl-1];
    mnrazz(ynpp1, fPbar, fSIMPy, jh, jl);
    for (i = 1; i <= fNpar; ++i) { fX[i-1] = fP[i + jl*fMaxpar - fMaxpar-1]; }
    mninex(fX);
    fCstatu = "PROGRESS  ";
    if (fISW[4] >= 1) mnprin(5, fAmin);
    fEDM  = fBigedm;
    sig2  = fEDM;
    ncycl = 0;
//*-*-                                       . . . . .  start main loop
L50:
    if (sig2 < fEpsi && fEDM < fEpsi) goto L76;
    sig2 = fEDM;
    if (fNfcn - npfn > fNfcnmx) goto L78;
//*-*-        calculate new point * by reflection
    for (i = 1; i <= fNpar; ++i) {
	pb = 0;
	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }
	fPbar[i-1]  = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];
	fPstar[i-1] = (alpha + 1)*fPbar[i-1] - alpha*fP[i + jh*fMaxpar - fMaxpar-1];
    }
    mninex(fPstar);
//    Eval(nparx, fGin, ystar, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, ystar, m_userParameterValue, 4);    ++fNfcn;
    if (ystar >= fAmin) goto L70;
//*-*-        point * better than jl, calculate new point **
    for (i = 1; i <= fNpar; ++i) {
	fPstst[i-1] = gamma*fPstar[i-1] + (1 - gamma)*fPbar[i-1];
    }
    mninex(fPstst);
//    Eval(nparx, fGin, ystst, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, ystst, m_userParameterValue, 4);    ++fNfcn;
//*-*-        try a parabola through ph, pstar, pstst.  min = prho
    y1 = (ystar - fSIMPy[jh-1])*rho2;
    y2 = (ystst - fSIMPy[jh-1])*rho1;
    rho = (rho2*y1 - rho1*y2)*.5 / (y1 - y2);
    if (rho < rhomin) goto L66;
    if (rho > rhomax) rho = rhomax;
    for (i = 1; i <= fNpar; ++i) {
	fPrho[i-1] = rho*fPbar[i-1] + (1 - rho)*fP[i + jh*fMaxpar - fMaxpar-1];
    }
    mninex(fPrho);
//    Eval(nparx, fGin, yrho, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, yrho, m_userParameterValue, 4);    ++fNfcn;
    if (yrho  < fSIMPy[jl-1] && yrho < ystst) goto L65;
    if (ystst < fSIMPy[jl-1]) goto L67;
    if (yrho  > fSIMPy[jl-1]) goto L66;
//*-*-        accept minimum point of parabola, PRHO
L65:
    mnrazz(yrho, fPrho, fSIMPy, jh, jl);
    goto L68;
L66:
    if (ystst < fSIMPy[jl-1]) goto L67;
    mnrazz(ystar, fPstar, fSIMPy, jh, jl);
    goto L68;
L67:
    mnrazz(ystst, fPstst, fSIMPy, jh, jl);
L68:
    ++ncycl;
    if (fISW[4] < 2) goto L50;
    if (fISW[4] >= 3 || ncycl % 10 == 0) {
	mnprin(5, fAmin);
    }
    goto L50;
//*-*-        point * is not as good as jl
L70:
    if (ystar >= fSIMPy[jh-1]) goto L73;
    jhold = jh;
    mnrazz(ystar, fPstar, fSIMPy, jh, jl);
    if (jhold != jh) goto L50;
//*-*-        calculate new point **
L73:
    for (i = 1; i <= fNpar; ++i) {
	fPstst[i-1] = beta*fP[i + jh*fMaxpar - fMaxpar-1] + (1 - beta)*fPbar[i-1];
    }
    mninex(fPstst);
//    Eval(nparx, fGin, ystst, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, ystst, m_userParameterValue, 4);    ++fNfcn;
    if (ystst > fSIMPy[jh-1]) goto L1;
//*-*-    point ** is better than jh
    if (ystst < fAmin) goto L67;
    mnrazz(ystst, fPstst, fSIMPy, jh, jl);
    goto L50;
//*-*-                                       . . . . . .  end main loop
L76:
    if (fISW[4] >= 0) {
	std::printf(" SIMPLEX MINIMIZATION HAS CONVERGED.\n");
    }
    fISW[3] = 1;
    goto L80;
L78:
    if (fISW[4] >= 0) {
	std::printf(" SIMPLEX TERMINATES WITHOUT CONVERGENCE.\n");
    }
    fCstatu = "CALL LIMIT";
    fISW[3] = -1;
    fISW[0] = 1;
L80:
    for (i = 1; i <= fNpar; ++i) {
	pb = 0;
	for (j = 1; j <= nparp1; ++j) { pb += wg*fP[i + j*fMaxpar - fMaxpar-1]; }
	fPbar[i-1] = pb - wg*fP[i + jh*fMaxpar - fMaxpar-1];
    }
    mninex(fPbar);
//    Eval(nparx, fGin, ypbar, fU, 4);    ++fNfcn;
    Eval(nparx, fGin, ypbar, m_userParameterValue, 4);    ++fNfcn;
    if (ypbar < fAmin) 	mnrazz(ypbar, fPbar, fSIMPy, jh, jl);
    mninex(fX);
    if (fNfcnmx + npfn - fNfcn < fNpar*3) goto L90;
    if (fEDM > fEpsi*2) goto L1;
L90:
    if (fISW[4] >= 0) mnprin(5, fAmin);
} /* mnsimp_ */

//______________________________________________________________________________
void URMinuit::mnstat(Double_urt &fmin, Double_urt &fedm, Double_urt &errdef, Int_urt &npari, Int_urt &nparx, Int_urt &istat)
{

//*-*-*-*-*Returns concerning the current status of the minimization*-*-*-*-*
//*-*      =========================================================
//*-*       User-called
//*-*          Namely, it returns:
//*-*        FMIN: the best function value found so far
//*-*        FEDM: the estimated vertical distance remaining to minimum
//*-*        ERRDEF: the value of UP defining parameter uncertainties
//*-*        NPARI: the number of currently variable parameters
//*-*        NPARX: the highest (external) parameter number defined by user
//*-*        ISTAT: a status integer indicating how good is the covariance
//*-*           matrix:  0= not calculated at all
//*-*                    1= approximation only, not accurate
//*-*                    2= full matrix, but forced positive-definite
//*-*                    3= full accurate covariance matrix
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    fmin   = fAmin;
    fedm   = fEDM;
    errdef = fUp;
    npari  = fNpar;
    nparx  = fNu;
    istat  = fISW[1];
    if (fEDM == fBigedm) fedm = fUp;
    if (fAmin == fUndefi) {
	fmin  = 0;
	fedm  = fUp;
	istat = 0;
    }
} /* mnstat_ */

//______________________________________________________________________________
void URMinuit::mntiny(Double_urt epsp1, Double_urt &epsbak)
{
//*-*-*-*-*-*-*-*To find the machine precision*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*            =============================
//*-*        Compares its argument with the value 1.0, and returns
//*-*        the value .TRUE. if they are equal.  To find EPSMAC
//*-*        safely by foiling the Fortran optimizer
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    epsbak = epsp1 - 1;
} /* mntiny_ */

//______________________________________________________________________________
Bool_urt URMinuit::mnunpt(const string &cfname)
{
//*-*-*-*-*-*Returns .TRUE. if CFNAME contains unprintable characters*-*-*-*
//*-*        ========================================================
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Int_urt i, l, ic;
    Bool_urt ret_val;
    static string cpt = " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890./;:[]$%*_!@#&+()";

    ret_val = kurFALSE;
    l       = strlen(cfname.c_str());
    for (i = 1; i <= l; ++i) {
	for (ic = 1; ic <= 80; ++ic) {
	    if (cfname[i-1] == cpt[ic-1]) goto L100;
	}
	return kurTRUE;
L100:
	;
    }
    return ret_val;
} /* mnunpt_ */

//______________________________________________________________________________
void URMinuit::mnvert(Double_urt *a, Int_urt l, Int_urt, Int_urt n, Int_urt &ifail)
{
//*-*-*-*-*-*-*-*-*-*-*-*Inverts a symmetric matrix*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    ==========================
//*-*        inverts a symmetric matrix.   matrix is first scaled to
//*-*        have all ones on the diagonal (equivalent to change of units)
//*-*        but no pivoting is done since matrix is positive-definite.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    /* System generated locals */
    Int_urt a_offset;

    /* Local variables */
    Double_urt si;
    Int_urt i, j, k, kp1, km1;

    /* Parameter adjustments */
    a_offset = l + 1;
    a -= a_offset;

    /* Function Body */
    ifail = 0;
    if (n < 1) goto L100;
    if (n > fMaxint) goto L100;
//*-*-                  scale matrix by sqrt of diag elements
    for (i = 1; i <= n; ++i) {
	si = a[i + i*l];
	if (si <= 0) goto L100;
	fVERTs[i-1] = 1 / URMath::Sqrt(si);
    }
    for (i = 1; i <= n; ++i) {
	for (j = 1; j <= n; ++j) {
	    a[i + j*l] = a[i + j*l]*fVERTs[i-1]*fVERTs[j-1];
	}
    }
//*-*-                                       . . . start main loop . . . .
    for (i = 1; i <= n; ++i) {
	k = i;
//*-*-                  preparation for elimination step1
	if (a[k + k*l] != 0) fVERTq[k-1] = 1 / a[k + k*l];
	else goto L100;
	fVERTpp[k-1] = 1;
	a[k + k*l] = 0;
	kp1 = k + 1;
	km1 = k - 1;
	if (km1 < 0) goto L100;
	else if (km1 == 0) goto L50;
	else               goto L40;
L40:
	for (j = 1; j <= km1; ++j) {
	    fVERTpp[j-1] = a[j + k*l];
	    fVERTq[j-1]  = a[j + k*l]*fVERTq[k-1];
	    a[j + k*l]   = 0;
	}
L50:
	if (k - n < 0) goto L51;
	else if (k - n == 0) goto L60;
	else                goto L100;
L51:
	for (j = kp1; j <= n; ++j) {
	    fVERTpp[j-1] = a[k + j*l];
	    fVERTq[j-1]  = -a[k + j*l]*fVERTq[k-1];
	    a[k + j*l]   = 0;
	}
//*-*-                  elimination proper
L60:
	for (j = 1; j <= n; ++j) {
	    for (k = j; k <= n; ++k) { a[j + k*l] += fVERTpp[j-1]*fVERTq[k-1]; }
	}
    }
//*-*-                  elements of left diagonal and unscaling
    for (j = 1; j <= n; ++j) {
	for (k = 1; k <= j; ++k) {
	    a[k + j*l] = a[k + j*l]*fVERTs[k-1]*fVERTs[j-1];
	    a[j + k*l] = a[k + j*l];
	}
    }
    return;
//*-*-                  failure return
L100:
    ifail = 1;
} /* mnvert_ */

//______________________________________________________________________________
void URMinuit::mnwarn(const char* copt1, const char* corg1, const char* cmes1)
{
//*-*-*-*-*-*-*-*-*-*-*-*Prints Warning messages*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                    =======================
//*-*     If COPT='W', CMES is a WARning message from CORG.
//*-*     If COPT='D', CMES is a DEBug message from CORG.
//*-*         If SET WARnings is in effect (the default), this routine
//*-*             prints the warning message CMES coming from CORG.
//*-*         If SET NOWarnings is in effect, the warning message is
//*-*             stored in a circular buffer of length MAXMES.
//*-*         If called with CORG=CMES='SHO', it prints the messages in
//*-*             the circular buffer, FIFO, and empties the buffer.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    string copt = copt1;
    string corg = corg1;
    string cmes = cmes1;

    const Int_urt MAXMES = 10;
    Int_urt ityp, i, ic, nm;
    string englsh, ctyp;

    if (corg.substr(0,3) != "SHO" || cmes.substr(0,3) != "SHO") {

//*-*-            Either print warning or put in buffer
       if (copt == "W") {
          ityp = 1;
	  if (fLwarn) {
             std::printf(" MINUIT WARNING IN %s\n",corg.c_str());
             std::printf(" ============== %s\n",cmes.c_str());
	     return;
	  }
       } else {
	  ityp = 2;
	  if (fLrepor) {
             std::printf(" MINUIT DEBUG FOR %s\n",corg.c_str());
             std::printf(" =============== %s \n",cmes.c_str());
	     return;
	  }
       }
//*-*-                if appropriate flag is off, fill circular buffer
       if (fNwrmes[ityp-1] == 0) fIcirc[ityp-1] = 0;
       ++fNwrmes[ityp-1];
       ++fIcirc[ityp-1];
       if (fIcirc[ityp-1] > 10) 	fIcirc[ityp-1] = 1;
       ic = fIcirc[ityp-1];
       fOrigin[ic] = corg;
       fWarmes[ic] = cmes;
       fNfcwar[ic] = fNfcn;
       return;
   }

//*-*-            'SHO WARnings', ask if any suppressed mess in buffer
    if (copt == "W") {
	ityp = 1;
	ctyp = "WARNING";
    } else {
	ityp = 2;
	ctyp = "*DEBUG*";
    }
    if (fNwrmes[ityp-1] > 0) {
	englsh = " WAS SUPPRESSED.  ";
	if (fNwrmes[ityp-1] > 1) englsh = "S WERE SUPPRESSED.";
        std::printf(" %5d MINUIT %s MESSAGE%s\n",fNwrmes[ityp-1]
             ,ctyp.c_str(),englsh.c_str());
	nm = fNwrmes[ityp-1];
	ic = 0;
	if (nm > MAXMES) {
	    std::printf(" ONLY THE MOST RECENT 10 WILL BE LISTED BELOW.\n");
	    nm = MAXMES;
	    ic = fIcirc[ityp-1];
	}
        std::printf("  CALLS  ORIGIN         MESSAGE\n");
	for (i = 1; i <= nm; ++i) {
	    ++ic;
	    if (ic > MAXMES) ic = 1;
            std::printf(" %6d  %s  %s\n", fNfcwar[ic],fOrigin[ic].c_str(),fWarmes[ic].c_str());
	}
	fNwrmes[ityp-1] = 0;
	std::printf(" \n");
    }
} /* mnwarn_ */

//______________________________________________________________________________
void URMinuit::mnwerr()
{
//*-*-*-*-*-*-*-*Calculates the WERR, external parameter errors*-*-*-*-*-*-*
//*-*            ==============================================
//*-*      and the global correlation coefficients, to be called
//*-*      whenever a new covariance matrix is available.
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    Double_urt denom, ba, al, dx, du1, du2;
    Int_urt ndex, ierr, i, j, k, l, ndiag, k1, iin;

//*-*-                        calculate external error if v exists
    if (fISW[1] >= 1) {
	for (l = 1; l <= fNpar; ++l) {
	    ndex = l*(l + 1) / 2;
	    dx = URMath::Sqrt(URMath::Abs(fVhmat[ndex-1]*fUp));
	    i = fNexofi[l-1];
//            if (fNvarl[i-1] > 1) {
	    if (m_userParameterFlag[i] > 1) {
		al = fAlim[i-1];
		ba = fBlim[i-1] - al;
//		du1 = al + 0.5*(URMath::Sin(fX[l-1] + dx) + 1)*ba - fU[i-1];
//		du2 = al + 0.5*(URMath::Sin(fX[l-1] - dx) + 1)*ba - fU[i-1];
                du1 = al + 0.5*(URMath::Sin(fX[l-1] + dx) + 1)*ba - m_userParameterValue[i];
                du2 = al + 0.5*(URMath::Sin(fX[l-1] - dx) + 1)*ba - m_userParameterValue[i];
		if (dx > 1) du1 = ba;
		dx = 0.5*(URMath::Abs(du1) + URMath::Abs(du2));
	    }
	    fWerr[l-1] = dx;
	}
    }
//*-*-                         global correlation coefficients
    if (fISW[1] >= 1) {
	for (i = 1; i <= fNpar; ++i) {
	    fGlobcc[i-1] = 0;
	    k1 = i*(i-1) / 2;
	    for (j = 1; j <= i; ++j) {
		k = k1 + j;
		fP[i + j*fMaxpar - fMaxpar-1] = fVhmat[k-1];
		fP[j + i*fMaxpar - fMaxpar-1] = fP[i + j*fMaxpar - fMaxpar-1];
	    }
	}
	mnvert(fP, fMaxint, fMaxint, fNpar, ierr);
	if (ierr == 0) {
	    for (iin = 1; iin <= fNpar; ++iin) {
		ndiag = iin*(iin + 1) / 2;
		denom = fP[iin + iin*fMaxpar - fMaxpar-1]*fVhmat[ndiag-1];
		if (denom <= 1 && denom >= 0) fGlobcc[iin-1] = 0;
		else                          fGlobcc[iin-1] = URMath::Sqrt(1 - 1 / denom);
	    }
	}
    }
} /* mnwerr_ */
