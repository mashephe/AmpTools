// @(#)root/base:$Name:  $:$Id: URMath.h,v 1.1.1.1 2006/11/15 18:01:50 mashephe Exp $
// Author: Fons Rademakers   29/07/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_URMath
#define ROOT_URMath


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// URMath                                                                //
//                                                                      //
// Encapsulate math routines. For the time being avoid templates.       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_URtypes
#include "UpRootMinuit/URtypes.h"
#endif

#include <string>

class URMath {

private:
   static Double_urt GamCf(Double_urt a,Double_urt x);
   static Double_urt GamSer(Double_urt a,Double_urt x);

public:

   // Fundamental constants
   static Double_urt Pi()       { return 3.14159265358979323846; }
   static Double_urt TwoPi()    { return 2.0 * Pi(); }
   static Double_urt PiOver2()  { return Pi() / 2.0; }
   static Double_urt PiOver4()  { return Pi() / 4.0; }
   static Double_urt InvPi()    { return 1.0 / Pi(); }
   static Double_urt RadToDeg() { return 180.0 / Pi(); }
   static Double_urt DegToRad() { return Pi() / 180.0; }

   // e (base of natural log)
   static Double_urt E()        { return 2.71828182845904523536; }

   // natural log of 10 (to convert log to ln)
   static Double_urt Ln10()     { return 2.30258509299404568402; }

   // base-10 log of e  (to convert ln to log)
   static Double_urt LogE()     { return 0.43429448190325182765; }

   // velocity of light
   static Double_urt C()        { return 2.99792458e8; }        // m s^-1
   static Double_urt Ccgs()     { return 100.0 * C(); }         // cm s^-1
   static Double_urt CUncertainty() { return 0.0; }             // exact

   // gravitational constant
   static Double_urt G()        { return 6.673e-11; }           // m^3 kg^-1 s^-2
   static Double_urt Gcgs()     { return G() / 1000.0; }        // cm^3 g^-1 s^-2
   static Double_urt GUncertainty() { return 0.010e-11; }

   // G over h-bar C
   static Double_urt GhbarC()   { return 6.707e-39; }           // (GeV/c^2)^-2
   static Double_urt GhbarCUncertainty() { return 0.010e-39; }

   // standard acceleration of gravity
   static Double_urt Gn()       { return 9.80665; }             // m s^-2
   static Double_urt GnUncertainty() { return 0.0; }            // exact

   // Planck's constant
   static Double_urt H()        { return 6.62606876e-34; }      // J s
   static Double_urt Hcgs()     { return 1.0e7 * H(); }         // erg s
   static Double_urt HUncertainty() { return 0.00000052e-34; }

   // h-bar (h over 2 pi)
   static Double_urt Hbar()     { return 1.054571596e-34; }     // J s
   static Double_urt Hbarcgs()  { return 1.0e7 * Hbar(); }      // erg s
   static Double_urt HbarUncertainty() { return 0.000000082e-34; }

   // hc (h * c)
   static Double_urt HC()       { return H() * C(); }           // J m
   static Double_urt HCcgs()    { return Hcgs() * Ccgs(); }     // erg cm

   // Boltzmann's constant
   static Double_urt K()        { return 1.3806503e-23; }       // J K^-1
   static Double_urt Kcgs()     { return 1.0e7 * K(); }         // erg K^-1
   static Double_urt KUncertainty() { return 0.0000024e-23; }

   // Stefan-Boltzmann constant
   static Double_urt Sigma()    { return 5.6704e-8; }           // W m^-2 K^-4
   static Double_urt SigmaUncertainty() { return 0.000040e-8; }

   // Avogadro constant (Avogadro's Number)
   static Double_urt Na()       { return 6.02214199e+23; }      // mol^-1
   static Double_urt NaUncertainty() { return 0.00000047e+23; }

   // universal gas constant (Na * K)
   // http://scienceworld.wolfram.com/physics/UniversalGasConstant.html
   static Double_urt R()        { return K() * Na(); }          // J K^-1 mol^-1
   static Double_urt RUncertainty() { return R()*((KUncertainty()/K()) + (NaUncertainty()/Na())); }

   // Molecular weight of dry air
   // 1976 US Standard Atmosphere,
   // also see http://atmos.nmsu.edu/jsdap/encyclopediawork.html
   static Double_urt MWair()    { return 28.9644; }             // kg kmol^-1 (or gm mol^-1)

   // Dry Air Gas Constant (R / MWair)
   // http://atmos.nmsu.edu/education_and_outreach/encyclopedia/gas_constant.htm
   static Double_urt Rgair()    { return (1000.0 * R()) / MWair(); }  // J kg^-1 K^-1

   // Elementary charge
   static Double_urt Qe()       { return 1.602176462e-19; }     // C
   static Double_urt QeUncertainty() { return 0.000000063e-19; }

   // Trigo
   static Double_urt Sin(Double_urt);
   static Double_urt Cos(Double_urt);
   static Double_urt Tan(Double_urt);
   static Double_urt SinH(Double_urt);
   static Double_urt CosH(Double_urt);
   static Double_urt TanH(Double_urt);
   static Double_urt ASin(Double_urt);
   static Double_urt ACos(Double_urt);
   static Double_urt ATan(Double_urt);
   static Double_urt ATan2(Double_urt, Double_urt);
   static Double_urt ASinH(Double_urt);
   static Double_urt ACosH(Double_urt);
   static Double_urt ATanH(Double_urt);
   static Double_urt Hypot(Double_urt x, Double_urt y);

   // Misc
   static Double_urt Sqrt(Double_urt x);
   static Double_urt Ceil(Double_urt x);
   static Double_urt Floor(Double_urt x);
   static Double_urt Exp(Double_urt);
   static Double_urt Factorial(Int_urt);
   static Double_urt Power(Double_urt x, Double_urt y);
   static Double_urt Log(Double_urt x);
   static Double_urt Log2(Double_urt x);
   static Double_urt Log10(Double_urt x);
   static Int_urt    Nint(Float_urt x);
   static Int_urt    Nint(Double_urt x);
   static Int_urt    Finite(Double_urt x);
   static Int_urt    IsNaN(Double_urt x);

   // Some integer math
   static Long_urt   NextPrime(Long_urt x);   // Least prime number greater than x
   static Long_urt   Sqrt(Long_urt x);
   static Long_urt   Hypot(Long_urt x, Long_urt y);     // sqrt(px*px + py*py)

   // Abs
   static Short_urt  Abs(Short_urt d);
   static Int_urt    Abs(Int_urt d);
   static Long_urt   Abs(Long_urt d);
   static Float_urt  Abs(Float_urt d);
   static Double_urt Abs(Double_urt d);

   // Even/Odd
   static Bool_urt Even(Long_urt a);
   static Bool_urt Odd(Long_urt a);

   // Sign
   static Short_urt  Sign(Short_urt a, Short_urt b);
   static Int_urt    Sign(Int_urt a, Int_urt b);
   static Long_urt   Sign(Long_urt a, Long_urt b);
   static Float_urt  Sign(Float_urt a, Float_urt b);
   static Double_urt Sign(Double_urt a, Double_urt b);

   // Min
   static Short_urt  Min(Short_urt a, Short_urt b);
   static UShort_urt Min(UShort_urt a, UShort_urt b);
   static Int_urt    Min(Int_urt a, Int_urt b);
   static UInt_urt   Min(UInt_urt a, UInt_urt b);
   static Long_urt   Min(Long_urt a, Long_urt b);
   static ULong_urt  Min(ULong_urt a, ULong_urt b);
   static Float_urt  Min(Float_urt a, Float_urt b);
   static Double_urt Min(Double_urt a, Double_urt b);

   // Max
   static Short_urt  Max(Short_urt a, Short_urt b);
   static UShort_urt Max(UShort_urt a, UShort_urt b);
   static Int_urt    Max(Int_urt a, Int_urt b);
   static UInt_urt   Max(UInt_urt a, UInt_urt b);
   static Long_urt   Max(Long_urt a, Long_urt b);
   static ULong_urt  Max(ULong_urt a, ULong_urt b);
   static Float_urt  Max(Float_urt a, Float_urt b);
   static Double_urt Max(Double_urt a, Double_urt b);

   // Locate Min, Max
   static Int_urt  LocMin(Int_urt n, const Short_urt *a);
   static Int_urt  LocMin(Int_urt n, const Int_urt *a);
   static Int_urt  LocMin(Int_urt n, const Float_urt *a);
   static Int_urt  LocMin(Int_urt n, const Double_urt *a);
   static Int_urt  LocMin(Int_urt n, const Long_urt *a);
   static Int_urt  LocMax(Int_urt n, const Short_urt *a);
   static Int_urt  LocMax(Int_urt n, const Int_urt *a);
   static Int_urt  LocMax(Int_urt n, const Float_urt *a);
   static Int_urt  LocMax(Int_urt n, const Double_urt *a);
   static Int_urt  LocMax(Int_urt n, const Long_urt *a);

   // Range
   static Short_urt  Range(Short_urt lb, Short_urt ub, Short_urt x);
   static Int_urt    Range(Int_urt lb, Int_urt ub, Int_urt x);
   static Long_urt   Range(Long_urt lb, Long_urt ub, Long_urt x);
   static ULong_urt  Range(ULong_urt lb, ULong_urt ub, ULong_urt x);
   static Double_urt Range(Double_urt lb, Double_urt ub, Double_urt x);

   // Binary search
   static Int_urt BinarySearch(Int_urt n, const Short_urt *array, Short_urt value);
   static Int_urt BinarySearch(Int_urt n, const Short_urt **array, Short_urt value);
   static Int_urt BinarySearch(Int_urt n, const Int_urt *array, Int_urt value);
   static Int_urt BinarySearch(Int_urt n, const Int_urt **array, Int_urt value);
   static Int_urt BinarySearch(Int_urt n, const Float_urt *array, Float_urt value);
   static Int_urt BinarySearch(Int_urt n, const Float_urt **array, Float_urt value);
   static Int_urt BinarySearch(Int_urt n, const Double_urt *array, Double_urt value);
   static Int_urt BinarySearch(Int_urt n, const Double_urt **array, Double_urt value);
   static Int_urt BinarySearch(Int_urt n, const Long_urt *array, Long_urt value);
   static Int_urt BinarySearch(Int_urt n, const Long_urt **array, Long_urt value);

   // Hashing
   static ULong_urt Hash(const void *txt, Int_urt ntxt);
   static ULong_urt Hash(const char *str);

   // IsInside
   static Bool_urt IsInside(Double_urt xp, Double_urt yp, Int_urt np, Double_urt *x, Double_urt *y);
   static Bool_urt IsInside(Float_urt xp, Float_urt yp, Int_urt np, Float_urt *x, Float_urt *y);
   static Bool_urt IsInside(Int_urt xp, Int_urt yp, Int_urt np, Int_urt *x, Int_urt *y);

   // Sorting
   static void Sort(Int_urt n, const Short_urt *a,  Int_urt *index, Bool_urt down=kurTRUE);
   static void Sort(Int_urt n, const Int_urt *a,    Int_urt *index, Bool_urt down=kurTRUE);
   static void Sort(Int_urt n, const Float_urt *a,  Int_urt *index, Bool_urt down=kurTRUE);
   static void Sort(Int_urt n, const Double_urt *a, Int_urt *index, Bool_urt down=kurTRUE);
   static void Sort(Int_urt n, const Long_urt *a,   Int_urt *index, Bool_urt down=kurTRUE);
   static void BubbleHigh(Int_urt Narr, Double_urt *arr1, Int_urt *arr2);
   static void BubbleLow (Int_urt Narr, Double_urt *arr1, Int_urt *arr2);

   // Advanced
   static Float_urt *Cross(Float_urt v1[3],Float_urt v2[3],Float_urt out[3]);     // Calculate the Cross Product of two vectors
   static Float_urt  Normalize(Float_urt v[3]);                               // Normalize a vector
   static Float_urt  NormCross(Float_urt v1[3],Float_urt v2[3],Float_urt out[3]); // Calculate the Normalized Cross Product of two vectors
   static Float_urt *Normal2Plane(Float_urt v1[3],Float_urt v2[3],Float_urt v3[3], Float_urt normal[3]); // Calcualte a normal vector of a plane

   static Double_urt *Cross(Double_urt v1[3],Double_urt v2[3],Double_urt out[3]);// Calculate the Cross Product of two vectors
   static Double_urt  Erf(Double_urt x);
   static Double_urt  Erfc(Double_urt x);
   static Double_urt  Freq(Double_urt x);
   static Double_urt  Gamma(Double_urt z);
   static Double_urt  Gamma(Double_urt a,Double_urt x);
   static Double_urt  BreitWigner(Double_urt x, Double_urt mean=0, Double_urt gamma=1);
   static Double_urt  Gaus(Double_urt x, Double_urt mean=0, Double_urt sigma=1, Bool_urt norm=kurFALSE);
   static Double_urt  Landau(Double_urt x, Double_urt mean=0, Double_urt sigma=1);
   static Double_urt  LnGamma(Double_urt z);
   static Double_urt  Normalize(Double_urt v[3]);                             // Normalize a vector
   static Double_urt  NormCross(Double_urt v1[3],Double_urt v2[3],Double_urt out[3]); // Calculate the Normalized Cross Product of two vectors
   static Double_urt *Normal2Plane(Double_urt v1[3],Double_urt v2[3],Double_urt v3[3], Double_urt normal[3]); // Calcualte a normal vector of a plane
   static Double_urt  Poisson(Double_urt x, Double_urt par);
   static Double_urt  Prob(Double_urt chi2,Int_urt ndf);
   static Double_urt  KolmogorovProb(Double_urt z);
   static Double_urt  Voigt(Double_urt x, Double_urt sigma, Double_urt lg, Int_urt R = 4);

   // Bessel functions
   static Double_urt BesselI(Int_urt n,Double_urt x);      // integer order modified Bessel function I_n(x)
   static Double_urt BesselK(Int_urt n,Double_urt x);      // integer order modified Bessel function K_n(x)
   static Double_urt BesselI0(Double_urt x);             // modified Bessel function I_0(x)
   static Double_urt BesselK0(Double_urt x);             // modified Bessel function K_0(x)
   static Double_urt BesselI1(Double_urt x);             // modified Bessel function I_1(x)
   static Double_urt BesselK1(Double_urt x);             // modified Bessel function K_1(x)
   static Double_urt BesselJ0(Double_urt x);             // Bessel function J0(x) for any real x
   static Double_urt BesselJ1(Double_urt x);             // Bessel function J1(x) for any real x
   static Double_urt BesselY0(Double_urt x);             // Bessel function Y0(x) for positive x
   static Double_urt BesselY1(Double_urt x);             // Bessel function Y1(x) for positive x
   static Double_urt StruveH0(Double_urt x);             // Struve functions of order 0
   static Double_urt StruveH1(Double_urt x);             // Struve functions of order 1
   static Double_urt StruveL0(Double_urt x);             // Modified Struve functions of order 0
   static Double_urt StruveL1(Double_urt x);             // Modified Struve functions of order 1

   // kludge for now
   void Error( const std::string& errorMessage, double value );
   void Error( const std::string& errorMessage, int iValue, double value );
};


//---- Even/odd ----------------------------------------------------------------

inline Bool_urt URMath::Even(Long_urt a)
   { return ! (a & 1); }

inline Bool_urt URMath::Odd(Long_urt a)
   { return (a & 1); }

//---- Abs ---------------------------------------------------------------------

inline Short_urt URMath::Abs(Short_urt d)
   { return (d > 0) ? d : -d; }

inline Int_urt URMath::Abs(Int_urt d)
   { return (d > 0) ? d : -d; }

inline Long_urt URMath::Abs(Long_urt d)
   { return (d > 0) ? d : -d; }

inline Float_urt URMath::Abs(Float_urt d)
   { return (d > 0) ? d : -d; }

inline Double_urt URMath::Abs(Double_urt d)
   { return (d > 0) ? d : -d; }

//---- Sign --------------------------------------------------------------------

inline Short_urt URMath::Sign(Short_urt a, Short_urt b)
   { return (b >= 0) ? Abs(a) : -Abs(a); }

inline Int_urt URMath::Sign(Int_urt a, Int_urt b)
   { return (b >= 0) ? Abs(a) : -Abs(a); }

inline Long_urt URMath::Sign(Long_urt a, Long_urt b)
   { return (b >= 0) ? Abs(a) : -Abs(a); }

inline Float_urt URMath::Sign(Float_urt a, Float_urt b)
   { return (b >= 0) ? Abs(a) : -Abs(a); }

inline Double_urt URMath::Sign(Double_urt a, Double_urt b)
   { return (b >= 0) ? Abs(a) : -Abs(a); }

//---- Min ---------------------------------------------------------------------

inline Short_urt URMath::Min(Short_urt a, Short_urt b)
   { return a <= b ? a : b; }

inline UShort_urt URMath::Min(UShort_urt a, UShort_urt b)
   { return a <= b ? a : b; }

inline Int_urt URMath::Min(Int_urt a, Int_urt b)
   { return a <= b ? a : b; }

inline UInt_urt URMath::Min(UInt_urt a, UInt_urt b)
   { return a <= b ? a : b; }

inline Long_urt URMath::Min(Long_urt a, Long_urt b)
   { return a <= b ? a : b; }

inline ULong_urt URMath::Min(ULong_urt a, ULong_urt b)
   { return a <= b ? a : b; }

inline Float_urt URMath::Min(Float_urt a, Float_urt b)
   { return a <= b ? a : b; }

inline Double_urt URMath::Min(Double_urt a, Double_urt b)
   { return a <= b ? a : b; }

//---- Max ---------------------------------------------------------------------

inline Short_urt URMath::Max(Short_urt a, Short_urt b)
   { return a >= b ? a : b; }

inline UShort_urt URMath::Max(UShort_urt a, UShort_urt b)
   { return a >= b ? a : b; }

inline Int_urt URMath::Max(Int_urt a, Int_urt b)
   { return a >= b ? a : b; }

inline UInt_urt URMath::Max(UInt_urt a, UInt_urt b)
   { return a >= b ? a : b; }

inline Long_urt URMath::Max(Long_urt a, Long_urt b)
   { return a >= b ? a : b; }

inline ULong_urt URMath::Max(ULong_urt a, ULong_urt b)
   { return a >= b ? a : b; }

inline Float_urt URMath::Max(Float_urt a, Float_urt b)
   { return a >= b ? a : b; }

inline Double_urt URMath::Max(Double_urt a, Double_urt b)
   { return a >= b ? a : b; }

//---- Range -------------------------------------------------------------------

inline Short_urt URMath::Range(Short_urt lb, Short_urt ub, Short_urt x)
   { return x < lb ? lb : (x > ub ? ub : x); }

inline Int_urt URMath::Range(Int_urt lb, Int_urt ub, Int_urt x)
   { return x < lb ? lb : (x > ub ? ub : x); }

inline Long_urt URMath::Range(Long_urt lb, Long_urt ub, Long_urt x)
   { return x < lb ? lb : (x > ub ? ub : x); }

inline ULong_urt URMath::Range(ULong_urt lb, ULong_urt ub, ULong_urt x)
   { return x < lb ? lb : (x > ub ? ub : x); }

inline Double_urt URMath::Range(Double_urt lb, Double_urt ub, Double_urt x)
   { return x < lb ? lb : (x > ub ? ub : x); }

//---- Trig and other functions ------------------------------------------------


#include <float.h>

#ifdef R__WIN32
#   ifndef finite
#      define finite _finite
#      define isnan  _isnan
#   endif
#endif
#if defined(R__AIX) || defined(R__MAC) || defined(R__SOLARIS_CC50) || \
    defined(R__HPUX11) || defined(R__GLIBC)
// math functions are defined inline so we have to include them here
#   include <math.h>
#   ifdef R__SOLARIS_CC50
       extern "C" { int finite(double); }
#   endif
#   if defined(R__GLIBC) && defined(__STRICT_ANSI__)
#      ifndef finite
#         define finite __finite
#      endif
#      ifndef isnan
#         define isnan  __isnan
#      endif
#   endif
#else
// don't want to include complete <math.h>
extern "C" {
   extern double sin(double);
   extern double cos(double);
   extern double tan(double);
   extern double sinh(double);
   extern double cosh(double);
   extern double tanh(double);
   extern double asin(double);
   extern double acos(double);
   extern double atan(double);
   extern double atan2(double, double);
   extern double sqrt(double);
   extern double exp(double);
   extern double pow(double, double);
   extern double log(double);
   extern double log10(double);
#ifndef R__WIN32
#   if !defined(finite)
       extern int finite(double);
#   endif
#   if !defined(isnan)
       extern int isnan(double);
#   endif
#endif
}
#endif

inline Double_urt URMath::Sin(Double_urt x)
   { return sin(x); }

inline Double_urt URMath::Cos(Double_urt x)
   { return cos(x); }

inline Double_urt URMath::Tan(Double_urt x)
   { return tan(x); }

inline Double_urt URMath::SinH(Double_urt x)
   { return sinh(x); }

inline Double_urt URMath::CosH(Double_urt x)
   { return cosh(x); }

inline Double_urt URMath::TanH(Double_urt x)
   { return tanh(x); }

inline Double_urt URMath::ASin(Double_urt x)
   { return asin(x); }

inline Double_urt URMath::ACos(Double_urt x)
   { return acos(x); }

inline Double_urt URMath::ATan(Double_urt x)
   { return atan(x); }

inline Double_urt URMath::ATan2(Double_urt y, Double_urt x)
   { return x != 0 ? atan2(y, x) : (y > 0 ? Pi()/2 : -Pi()/2); }

inline Double_urt URMath::Sqrt(Double_urt x)
   { return sqrt(x); }

inline Double_urt URMath::Exp(Double_urt x)
   { return exp(x); }

inline Double_urt URMath::Power(Double_urt x, Double_urt y)
   { return pow(x, y); }

inline Double_urt URMath::Log(Double_urt x)
   { return log(x); }

inline Double_urt URMath::Log10(Double_urt x)
   { return log10(x); }

inline Int_urt URMath::Finite(Double_urt x)
#ifdef R__HPUX11
   { return isfinite(x); }
#else
   { return finite(x); }
#endif

inline Int_urt URMath::IsNaN(Double_urt x)
   { return isnan(x); }

//-------- Advanced -------------

inline Float_urt URMath::NormCross(Float_urt v1[3],Float_urt v2[3],Float_urt out[3])
{
   // Calculate the Normalized Cross Product of two vectors
   return Normalize(Cross(v1,v2,out));
}

inline Double_urt URMath::NormCross(Double_urt v1[3],Double_urt v2[3],Double_urt out[3])
{
   // Calculate the Normalized Cross Product of two vectors
   return Normalize(Cross(v1,v2,out));
}


#endif
