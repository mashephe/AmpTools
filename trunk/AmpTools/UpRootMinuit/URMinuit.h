// @(#)root/minuit:$Name:  $:$Id: URMinuit.h,v 1.2 2009/05/15 23:58:05 hrayr Exp $
// Author: Rene Brun, Frederick James   12/08/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// ---------------------------------- minuit.h



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// URMinuit                                                              //
//                                                                      //
// The MINUIT minimisation package (base class)                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef UPROOT_URMinuit
#define UPROOT_URMinuit

#include <string.h>
#include <vector>
#include <ostream>
#include <algorithm>

#include "UpRootMinuit/URtypes.h"
#include "UpRootMinuit/URFcn.h"

class URMinuit {

private:
   URMinuit(const URMinuit &m); // prevent default

// should become private....
public:
        enum{kMAXWARN=100};
          
        Int_urt        fNpfix;            //Number of fixed parameters
        Int_urt        fEmpty;            //Initialization flag (1 = Minuit initialized)
        Int_urt        fMaxpar;           //Maximum number of parameters
        Int_urt        fMaxint;           //Maximum number of internal parameters
        Int_urt        fNpar;             //Number of free parameters (total number of pars = fNpar + fNfix)
        Int_urt        fMaxext;           //Maximum number of external parameters
        Int_urt        fMaxIterations;    //Maximum number of iterations
        Int_urt        fMaxpar5;          // fMaxpar*(fMaxpar+1)/2
        Int_urt        fMaxcpt;
        Int_urt        fMaxpar2;          // fMaxpar*fMaxpar
        Int_urt        fMaxpar1;          // fMaxpar*(fMaxpar+1)
        
        Double_urt     fAmin;             //Minimum value found for FCN
        Double_urt     fUp;               //FCN+-UP defines errors (for chisquare fits UP=1)
        Double_urt     fEDM;              //Estimated vertical distance to the minimum
        Double_urt     fFval3;            //
        Double_urt     fEpsi;             //
        Double_urt     fApsi;             //
        Double_urt     fDcovar;           //Relative change in covariance matrix
        Double_urt     fEpsmac;           //machine precision for floating points:
        Double_urt     fEpsma2;           //sqrt(fEpsmac)
        Double_urt     fVlimlo;           //
        Double_urt     fVlimhi;           //
        Double_urt     fUndefi;           //Undefined number = -54321
        Double_urt     fBigedm;           //Big EDM = 123456
        Double_urt     fUpdflt;           //
        Double_urt     fXmidcr;           //
        Double_urt     fYmidcr;           //
        Double_urt     fXdircr;           //
        Double_urt     fYdircr;           //
        
//        Double_urt     *fU;               //[fMaxpar2] External (visible to user in FCN) value of parameters
        Double_urt     *fAlim;            //[fMaxpar2] Lower limits for parameters. If zero no limits
        Double_urt     *fBlim;            //[fMaxpar2] Upper limits for parameters
        Double_urt     *fErp;             //[fMaxpar] Positive Minos errors if calculated
        Double_urt     *fErn;             //[fMaxpar] Negative Minos errors if calculated
        Double_urt     *fWerr;            //[fMaxpar] External parameters error (standard deviation, defined by UP)
        Double_urt     *fGlobcc;          //[fMaxpar] Global Correlation Coefficients
        Double_urt     *fX;               //[fMaxpar] Internal parameters values
        Double_urt     *fXt;              //[fMaxpar] Internal parameters values X saved as Xt
        Double_urt     *fDirin;           //[fMaxpar] (Internal) step sizes for current step
        Double_urt     *fXs;              //[fMaxpar] Internal parameters values saved for fixed params
        Double_urt     *fXts;             //[fMaxpar] Internal parameters values X saved as Xt for fixed params
        Double_urt     *fDirins;          //[fMaxpar] (Internal) step sizes for current step for fixed params
        Double_urt     *fGrd;             //[fMaxpar] First derivatives
        Double_urt     *fG2;              //[fMaxpar] 
        Double_urt     *fGstep;           //[fMaxpar] Step sizes
        Double_urt     *fGin;             //[fMaxpar2] 
        Double_urt     *fDgrd;            //[fMaxpar] Uncertainties
        Double_urt     *fGrds;            //[fMaxpar] 
        Double_urt     *fG2s;             //[fMaxpar] 
        Double_urt     *fGsteps;          //[fMaxpar] 
        Double_urt     *fVhmat;           //[fMaxpar5] (Internal) error matrix stored as Half MATrix, since it is symmetric
        Double_urt     *fSecDer;          //[fMaxpar5] matrix of second derivatives filled in HESSE
        Double_urt     *fVthmat;          //[fMaxpar5] VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
        Double_urt     *fP;               //[fMaxpar1] 
        Double_urt     *fPstar;           //[fMaxpar2] 
        Double_urt     *fPstst;           //[fMaxpar] 
        Double_urt     *fPbar;            //[fMaxpar] 
        Double_urt     *fPrho;            //[fMaxpar] Minimum point of parabola
        Double_urt     *fWord7;           //[fMaxpar] 
        Double_urt     *fXpt;             //[fMaxcpt] X array of points for contours
        Double_urt     *fYpt;             //[fMaxcpt] Y array of points for contours
        
        Double_urt     *fCONTgcc;         //[fMaxpar] array used in mncont
        Double_urt     *fCONTw;           //[fMaxpar] array used in mncont
        Double_urt     *fFIXPyy;          //[fMaxpar] array used in mnfixp
        Double_urt     *fGRADgf;          //[fMaxpar] array used in mngrad
        Double_urt     *fHESSyy;          //[fMaxpar] array used in mnhess
        Double_urt     *fIMPRdsav;        //[fMaxpar] array used in mnimpr
        Double_urt     *fIMPRy;           //[fMaxpar] array used in mnimpr
        Double_urt     *fMATUvline;       //[fMaxpar] array used in mnmatu
        Double_urt     *fMIGRflnu;        //[fMaxpar] array used in mnmigr
        Double_urt     *fMIGRstep;        //[fMaxpar] array used in mnmigr
        Double_urt     *fMIGRgs;          //[fMaxpar] array used in mnmigr
        Double_urt     *fMIGRvg;          //[fMaxpar] array used in mnmigr
        Double_urt     *fMIGRxxs;         //[fMaxpar] array used in mnmigr
        Double_urt     *fMNOTxdev;        //[fMaxpar] array used in mnmnot
        Double_urt     *fMNOTw;           //[fMaxpar] array used in mnmnot
        Double_urt     *fMNOTgcc;         //[fMaxpar] array used in mnmnot
        Double_urt     *fPSDFs;           //[fMaxpar] array used in mnpsdf
        Double_urt     *fSEEKxmid;        //[fMaxpar] array used in mnseek
        Double_urt     *fSEEKxbest;       //[fMaxpar] array used in mnseek
        Double_urt     *fSIMPy;           //[fMaxpar] array used in mnsimp
        Double_urt     *fVERTq;           //[fMaxpar] array used in mnvert
        Double_urt     *fVERTs;           //[fMaxpar] array used in mnvert
        Double_urt     *fVERTpp;          //[fMaxpar] array used in mnvert
        Double_urt     *fCOMDplist;       //[fMaxpar] array used in mncomd
        Double_urt     *fPARSplist;       //[fMaxpar] array used in mnpars
        
//        Int_urt        *fNvarl;           //[fMaxpar2] parameters flag (-1=undefined, 0=constant..)
//        Int_urt        *fNiofex;          //[fMaxpar2] Internal parameters number, or zero if not currently variable
        Int_urt        *fNexofi;          //[fMaxpar] External parameters number for currently variable parameters
        Int_urt        *fIpfix;           //[fMaxpar] List of fixed parameters
        Int_urt        fNu;               //
        Int_urt        fIsysrd;           //standardInput unit
        Int_urt        fIsyswr;           //standard output unit
        Int_urt        fIsyssa;           //
        Int_urt        fNpagwd;           //Page width
        Int_urt        fNpagln;           //Number of lines per page
        Int_urt        fNewpag;           //
        Int_urt        fIstkrd[10];       //
        Int_urt        fNstkrd;           //
        Int_urt        fIstkwr[10];       //
        Int_urt        fNstkwr;           //
        Int_urt        fISW[7];           //Array of switches
        Int_urt        fIdbg[11];         //Array of internal debug switches
        Int_urt        fNblock;           //Number of Minuit data blocks
        Int_urt        fIcomnd;           //Number of commands
        Int_urt        fNfcn;             //Number of calls to FCN
        Int_urt        fNfcnmx;           //Maximum number of calls to FCN
        Int_urt        fNfcnlc;           //
        Int_urt        fNfcnfr;           //
        Int_urt        fItaur;            //
        Int_urt        fIstrat;           //
        Int_urt        fNwrmes[2];        //
        Int_urt        fNfcwar[20];       //
        Int_urt        fIcirc[2];         //
        Int_urt        fStatus;           //Status flag for the last called Minuit function
        Int_urt        fKe1cr;            //
        Int_urt        fKe2cr;            //
        Bool_urt       fLwarn;            //true if warning messges are to be put out (default=true)
        Bool_urt       fLrepor;           //true if exceptional conditions are put out (default=false)
        Bool_urt       fLimset;           //true if a parameter is up against limits (for MINOS)
        Bool_urt       fLnolim;           //true if there are no limits on any parameters (not yet used)
        Bool_urt       fLnewmn;           //true if the previous process has unexpectedly improved FCN
        Bool_urt       fLphead;           //true if a heading should be put out for the next parameter definition
        char         *fChpt;            //!Character to be plotted at the X,Y contour positions
//        std::string      *fCpnam;           //[fMaxpar2] Array of parameters names
        std::string      fCfrom;            //
        std::string      fCstatu;           //
        std::string      fCtitl;            //
        std::string      fCword;            //
        std::string      fCundef;           //
        std::string      fCvrsn;            //
        std::string      fCovmes[4];        //
        std::string      fOrigin[kMAXWARN]; //
        std::string      fWarmes[kMAXWARN]; //
        URFcn* fFCN; // function to be minimized         
           // (*fFCN)(Int_urt &npar, Double_urt *gin, Double_urt &f, Double_urt *u, Int_urt flag); //!

// methods performed on URMinuit class
public:
                URMinuit();
                URMinuit(Int_urt maxpar);
 virtual       ~URMinuit();
 virtual void   BuildArrays(Int_urt maxpar=15);
 virtual void   zeroPointers();
 virtual Int_urt  Command(const char *command);
 virtual Int_urt  DefineParameter( Int_urt parNo, 
                                 const std::string& name, 
                                 Double_urt initVal, 
                                 Double_urt initErr, 
                                 Double_urt lowerLimit, 
                                 Double_urt upperLimit );
 virtual void   DeleteArrays();
 virtual Int_urt  Eval(Int_urt npar, Double_urt *grad, Double_urt &fval, const std::vector<Double_urt>& par, Int_urt flag);
 virtual Int_urt  FixParameter( Int_urt parNo );
 bool           parameterFixed( Int_urt parameterNumber );
 Int_urt          GetMaxIterations() const {return fMaxIterations;}
 virtual Int_urt  GetNumFixedPars() const;
 virtual Int_urt  GetNumFreePars() const;
 virtual Int_urt  GetNumPars() const;
 virtual Int_urt  GetParameter( Int_urt parNo, Double_urt &currentValue, Double_urt &currentError ) const;
 Int_urt          GetStatus() const {return fStatus;}
 const std::vector<Double_urt>& GetParameterList() const {return m_userParameterValue;}
 virtual Int_urt  Migrad();
 virtual Int_urt  Minos();
 virtual Int_urt  Hesse();
 void             SetLogStream( std::ostream& aStream ) { m_logStream = &aStream; }
 virtual void   mnamin();
 virtual void   mnbins(Double_urt a1, Double_urt a2, Int_urt naa, Double_urt &bl, Double_urt &bh, Int_urt &nb, Double_urt &bwid);
 virtual void   mncalf(Double_urt *pvec, Double_urt &ycalf);
 virtual void   mncler();
 virtual void   mncntr(Int_urt ke1, Int_urt ke2, Int_urt &ierrf);
 virtual void   mncomd(const std::string& crdbin, Int_urt &icondn);
 virtual void   mncont(Int_urt ke1, Int_urt ke2, Int_urt nptu, Double_urt *xptu, Double_urt *yptu, Int_urt &ierrf);
// virtual void   mncrck(const std::string& crdbuf, Int_urt maxcwd, const std::string &comand, Int_urt &lnc
// header change (mrs43): 
//   seems as if the point of this is to modify comand - should not be a const reference
//   cardbuf also appears to be an input only -- pass by value instead of const reference
 virtual void   mncrck(std::string crdbuf, Int_urt maxcwd, std::string &comand, Int_urt &lnc,
		Int_urt mxp, Double_urt *plist, Int_urt &llist, Int_urt &ierr, Int_urt isyswr);
 virtual void   mncros(Double_urt &aopt, Int_urt &iercr);
 virtual void   mncuve();
 virtual void   mnderi();
 virtual void   mndxdi(Double_urt pint, Int_urt ipar, Double_urt &dxdi);
 virtual void   mneig(Double_urt *a, Int_urt ndima, Int_urt n, Int_urt mits, Double_urt *work, Double_urt precis, Int_urt &ifault);
 virtual void   mnemat(Double_urt *emat, Int_urt ndim);
 virtual void   mnerrs(Int_urt number, Double_urt &eplus, Double_urt &eminus, Double_urt &eparab, Double_urt &gcc);
 virtual void   mneval(Double_urt anext, Double_urt &fnext, Int_urt &ierev);
 virtual void   mnexcm(const std::string& command, Double_urt *plist, Int_urt llist, Int_urt &ierflg) ;
 virtual void   mnexin(Double_urt *pint);
 virtual void   mnfixp(Int_urt iint, Int_urt &ierr);
 virtual void   mnfree(Int_urt k);
 virtual void   mngrad();
 virtual void   mnhelp(const std::string& comd);
 virtual void   mnhelp();
 virtual void   mnhess();
 virtual void   mnhes1();
 virtual void   mnimpr();
 virtual void   mninex(Double_urt *pint);
 virtual void   mninit(Int_urt i1, Int_urt i2, Int_urt i3);
 virtual void   mnlims();
 virtual void   mnline(Double_urt *start, Double_urt fstart, Double_urt *step, Double_urt slope, Double_urt toler);
 virtual void   mnmatu(Int_urt kode);
 virtual void   mnmigr();
 virtual void   mnmnos();
 virtual void   mnmnot(Int_urt ilax, Int_urt ilax2, Double_urt &val2pl, Double_urt &val2mi);
 virtual void   mnparm(Int_urt k, const std::string& cnamj, Double_urt uk, Double_urt wk, Double_urt a, Double_urt b, Int_urt &ierflg);
 virtual void   mnpars(const std::string& crdbuf, Int_urt &icondn);
 virtual void   mnpfit(Double_urt *parx2p, Double_urt *pary2p, Int_urt npar2p, Double_urt *coef2p, Double_urt &sdev2p);
 virtual void   mnpint(Double_urt &pexti, Int_urt i, Double_urt &pinti);
 virtual void   mnplot(Double_urt *xpt, Double_urt *ypt, char *chpt, Int_urt nxypt, Int_urt npagwd, Int_urt npagln);
 virtual void   mnpout(Int_urt iuext, std::string& chnam, Double_urt &val, Double_urt &err, Double_urt &xlolim, Double_urt &xuplim, Int_urt &iuint) const;
 virtual void   mnprin(Int_urt inkode, Double_urt fval);
 virtual void   mnpsdf();
 virtual void   mnrazz(Double_urt ynew, Double_urt *pnew, Double_urt *y, Int_urt &jh, Int_urt &jl);
 virtual void   mnrn15(Double_urt &val, Int_urt &inseed);
 virtual void   mnrset(Int_urt iopt);
 virtual void   mnsave();
 virtual void   mnscan();
 virtual void   mnseek();
 virtual void   mnset();
 virtual void   mnsimp();
 virtual void   mnstat(Double_urt &fmin, Double_urt &fedm, Double_urt &errdef, Int_urt &npari, Int_urt &nparx, Int_urt &istat);
 virtual void   mntiny(Double_urt epsp1, Double_urt &epsbak);
 Bool_urt         mnunpt(const std::string& cfname);
 virtual void   mnvert(Double_urt *a, Int_urt l, Int_urt m, Int_urt n, Int_urt &ifail);
 void           mnwarn( const char* copt1,  const char* corg1,  const char* cmes1);
 virtual void   mnwerr();
 virtual Int_urt  Release( Int_urt parNo );
 virtual Int_urt  SetErrorDef( Double_urt up );
// virtual void   SetFCN(void (*fcn)(Int_urt &, Double_urt *, Double_urt &f, Double_urt *, Int_urt));
 virtual void   SetFCN(URFcn* fcn);
 virtual void   SetMaxIterations(Int_urt maxiter=5000) {fMaxIterations = maxiter;}
 virtual Int_urt  SetPrintLevel( Int_urt printLevel=0 );

   private:
 
    // External (user) parameter information
    std::vector<Int_urt>       m_userParameterIdToInternalId; //Internal parameters number, or zero if not currently variable (was fNiofex)
    std::vector<Double_urt>    m_userParameterValue; //  External (visible to user in FCN) value of parameters (was fU)
    std::vector<std::string> m_userParameterName;  // parameters names (was fCpnam)
    std::vector<Int_urt>       m_userParameterFlag;   // parameter flags (-1=undefined, 0=constant..)  (was fNvarl)
    std::ostream* m_logStream;
};

#endif

