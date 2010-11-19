// @(#)root/base:$Name:  $:$Id: URMath.cc,v 1.1.1.1 2006/11/15 18:01:50 mashephe Exp $
// Author: Fons Rademakers   29/07/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// URMath                                                                //
//                                                                      //
// Encapsulate math routines (i.e. provide a kind of namespace).        //
// For the time being avoid templates.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "UpRootMinuit/URMath.h"
//#include "TError.h"
#include <math.h>
#include <string.h>

//const Double_urt
//   URMath::Pi = 3.14159265358979323846,
//   URMath::E  = 2.7182818284590452354;

#if defined(hocuspocus)
//______________________________________________________________________________
#if defined(R__MAC) || defined(R__KCC)
Double_urt hypot(Double_urt x, Double_urt y) {
   return sqrt(x*x+y*y);
}
#endif

//______________________________________________________________________________
Long_urt URMath::Sqrt(Long_urt x)
{
   return (Long_urt) (sqrt((Double_urt)x) + 0.5);
}

//______________________________________________________________________________
Long_urt URMath::Hypot(Long_urt x, Long_urt y)
{
   return (Long_urt) (hypot((Double_urt)x, (Double_urt)y) + 0.5);
}

//______________________________________________________________________________
Double_urt URMath::Hypot(Double_urt x, Double_urt y)
{
   return hypot(x, y);
}

//______________________________________________________________________________
Double_urt URMath::ASinH(Double_urt x)
{
#if defined(WIN32) || defined(R__KCC)
   return log(x+sqrt(x*x+1));
#else
   return asinh(x);
#endif
}

//______________________________________________________________________________
Double_urt URMath::ACosH(Double_urt x)
{
#if defined(WIN32) || defined(R__KCC)
   return log(x+sqrt(x*x-1));
#else
   return acosh(x);
#endif
}

//______________________________________________________________________________
Double_urt URMath::ATanH(Double_urt x)
{
#if defined(WIN32) || defined(R__KCC)
   return log((1+x)/(1-x))/2;
#else
   return atanh(x);
#endif
}

//______________________________________________________________________________
Double_urt URMath::Ceil(Double_urt x)
{
   return ceil(x);
}

//______________________________________________________________________________
Double_urt URMath::Floor(Double_urt x)
{
   return floor(x);
}

//______________________________________________________________________________
Double_urt URMath::Log2(Double_urt x)
{
   return log(x)/log(2.0);
}

//______________________________________________________________________________
Long_urt URMath::NextPrime(Long_urt x)
{
   // Return next prime number after x, unless x is a prime in which case
   // x is returned.

   if (x < 2)
      return 2;
   if (x == 3)
      return 3;

   if (x % 2 == 0)
      x++;

   Long_urt sqr = (Long_urt) sqrt((Double_urt)x) + 1;

   for (;;) {
      Long_urt n;
      for (n = 3; (n <= sqr) && ((x % n) != 0); n += 2)
         ;
      if (n > sqr)
         return x;
      x += 2;
   }
}

//______________________________________________________________________________
Int_urt URMath::Nint(Float_urt x)
{
   // Round to nearest integer. Rounds half integers to the nearest
   // even integer.

   int i;
   if (x >= 0) {
      i = int(x + 0.5);
      if (x + 0.5 == Float_urt(i) && i & 1) i--;
   } else {
      i = int(x - 0.5);
      if (x - 0.5 == Float_urt(i) && i & 1) i++;

   }
   return i;
}

//______________________________________________________________________________
Int_urt URMath::Nint(Double_urt x)
{
   // Round to nearest integer. Rounds half integers to the nearest
   // even integer.

   int i;
   if (x >= 0) {
      i = int(x + 0.5);
      if (x + 0.5 == Double_urt(i) && i & 1) i--;
   } else {
      i = int(x - 0.5);
      if (x - 0.5 == Double_urt(i) && i & 1) i++;

   }
   return i;
}

//______________________________________________________________________________
Float_urt *URMath::Cross(Float_urt v1[3],Float_urt v2[3],Float_urt out[3])
{
   // Calculate the Cross Product of two vectors:
   //         out = [v1 x v2]

    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return out;
}

//______________________________________________________________________________
Double_urt *URMath::Cross(Double_urt v1[3],Double_urt v2[3],Double_urt out[3])
{
   // Calculate the Cross Product of two vectors:
   //   out = [v1 x v2]

    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return out;
}

//______________________________________________________________________________
Double_urt URMath::Erf(Double_urt x)
{
   // Computation of the error function erf(x).
   // Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x

   //--- NvE 14-nov-1998 UU-SAP Utrecht

   return (1-Erfc(x));
}

//______________________________________________________________________________
Double_urt URMath::Erfc(Double_urt x)
{
   // Compute the complementary error function erfc(x).
   // Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   // The parameters of the Chebyshev fit
   const Double_urt a1 = -1.26551223,   a2 = 1.00002368,
                  a3 =  0.37409196,   a4 = 0.09678418,
                  a5 = -0.18628806,   a6 = 0.27886807,
                  a7 = -1.13520398,   a8 = 1.48851587,
                  a9 = -0.82215223,  a10 = 0.17087277;

   Double_urt v = 1; // The return value
   Double_urt z = Abs(x);

   if (z <= 0) return v; // erfc(0)=1

   Double_urt t = 1/(1+0.5*z);

   v = t*Exp((-z*z) +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));

   if (x < 0) v = 2-v; // erfc(-x)=2-erfc(x)

   return v;
}

//______________________________________________________________________________
Double_urt URMath::Factorial(Int_urt n)
{
   // Compute factorial(n)

  if(n <= 0) return 1.;
  Double_urt x=1;
  Int_urt b=0;
  do {
     b++;
     x *= b;
  } while(b!=n);
  return x;
}

//______________________________________________________________________________
Double_urt URMath::Freq(Double_urt x)
{
   // Computation of the normal frequency function freq(x).
   // Freq(x) = (1/sqrt(2pi)) Integral(exp(-t^2/2))dt between -infinity and x.
   //
   // Translated from CERNLIB C300 by Rene Brun.

   const Double_urt C1 = 0.56418958354775629;
   const Double_urt W2 = 1.41421356237309505;

   const Double_urt p10 = 2.4266795523053175e+2,  q10 = 2.1505887586986120e+2,
                  p11 = 2.1979261618294152e+1,  q11 = 9.1164905404514901e+1,
                  p12 = 6.9963834886191355e+0,  q12 = 1.5082797630407787e+1,
                  p13 =-3.5609843701815385e-2,  q13 = 1;

   const Double_urt p20 = 3.00459261020161601e+2, q20 = 3.00459260956983293e+2,
                  p21 = 4.51918953711872942e+2, q21 = 7.90950925327898027e+2,
                  p22 = 3.39320816734343687e+2, q22 = 9.31354094850609621e+2,
                  p23 = 1.52989285046940404e+2, q23 = 6.38980264465631167e+2,
                  p24 = 4.31622272220567353e+1, q24 = 2.77585444743987643e+2,
                  p25 = 7.21175825088309366e+0, q25 = 7.70001529352294730e+1,
                  p26 = 5.64195517478973971e-1, q26 = 1.27827273196294235e+1,
                  p27 =-1.36864857382716707e-7, q27 = 1;

   const Double_urt p30 =-2.99610707703542174e-3, q30 = 1.06209230528467918e-2,
                  p31 =-4.94730910623250734e-2, q31 = 1.91308926107829841e-1,
                  p32 =-2.26956593539686930e-1, q32 = 1.05167510706793207e+0,
                  p33 =-2.78661308609647788e-1, q33 = 1.98733201817135256e+0,
                  p34 =-2.23192459734184686e-2, q34 = 1;

   Double_urt v  = URMath::Abs(x)/W2;
   Double_urt vv = v*v;
   Double_urt ap, aq, h, hc, y;
   if (v < 0.5) {
      y=vv;
      ap=p13;
      aq=q13;
      ap = p12 +y*ap;
      ap = p11 +y*ap;
      ap = p10 +y*ap;
      aq = q12 +y*aq;
      aq = q11 +y*aq;
      aq = q10 +y*aq;
      h  = v*ap/aq;
      hc = 1-h;
   } else if (v < 4) {
      ap = p27;
      aq = q27;
      ap = p26 +v*ap;
      ap = p25 +v*ap;
      ap = p24 +v*ap;
      ap = p23 +v*ap;
      ap = p22 +v*ap;
      ap = p21 +v*ap;
      ap = p20 +v*ap;
      aq = q26 +v*aq;
      aq = q25 +v*aq;
      aq = q24 +v*aq;
      aq = q23 +v*aq;
      aq = q22 +v*aq;
      aq = q21 +v*aq;
      aq = q20 +v*aq;
      hc = URMath::Exp(-vv)*ap/aq;
      h  = 1-hc;
   } else {
      y  = 1/vv;
      ap = p34;
      aq = q34;
      ap = p33 +y*ap;
      ap = p32 +y*ap;
      ap = p31 +y*ap;
      ap = p30 +y*ap;
      aq = q33 +y*aq;
      aq = q32 +y*aq;
      aq = q31 +y*aq;
      aq = q30 +y*aq;
      hc = URMath::Exp(-vv)*(C1+y*ap/aq)/v;
      h  = 1-hc;
   }
   if (x > 0) return 0.5 +0.5*h;
   else return 0.5*hc;
}

//______________________________________________________________________________
Double_urt URMath::Gamma(Double_urt z)
{
   // Computation of gamma(z) for all z>0.
   //
   // C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   if (z<=0) return 0;

   Double_urt v = LnGamma(z);
   return Exp(v);
}

//______________________________________________________________________________
Double_urt URMath::Gamma(Double_urt a,Double_urt x)
{
   // Computation of the incomplete gamma function P(a,x)
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   if (a <= 0 || x <= 0) return 0;

   if (x < (a+1)) return GamSer(a,x);
   else           return GamCf(a,x);
}

//______________________________________________________________________________
Double_urt URMath::GamCf(Double_urt a,Double_urt x)
{
   // Computation of the incomplete gamma function P(a,x)
   // via its continued fraction representation.
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   Int_urt itmax    = 100;      // Maximum number of iterations
   Double_urt eps   = 3.e-7;    // Relative accuracy
   Double_urt fpmin = 1.e-30;   // Smallest Double_urt value allowed here

   if (a <= 0 || x <= 0) return 0;

   Double_urt gln = LnGamma(a);
   Double_urt b   = x+1-a;
   Double_urt c   = 1/fpmin;
   Double_urt d   = 1/b;
   Double_urt h   = d;
   Double_urt an,del;
   for (Int_urt i=1; i<=itmax; i++) {
      an = Double_urt(-i)*(Double_urt(i)-a);
      b += 2;
      d  = an*d+b;
      if (Abs(d) < fpmin) d = fpmin;
      c = b+an/c;
      if (Abs(c) < fpmin) c = fpmin;
      d   = 1/d;
      del = d*c;
      h   = h*del;
      if (Abs(del-1) < eps) break;
      //if (i==itmax) cout << "*GamCf(a,x)* a too large or itmax too small" << endl;
   }
   Double_urt v = Exp(-x+a*Log(x)-gln)*h;
   return (1-v);
}

//______________________________________________________________________________
Double_urt URMath::GamSer(Double_urt a,Double_urt x)
{
   // Computation of the incomplete gamma function P(a,x)
   // via its series representation.
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   Int_urt itmax  = 100;   // Maximum number of iterations
   Double_urt eps = 3.e-7; // Relative accuracy

   if (a <= 0 || x <= 0) return 0;

   Double_urt gln = LnGamma(a);
   Double_urt ap  = a;
   Double_urt sum = 1/a;
   Double_urt del = sum;
   for (Int_urt n=1; n<=itmax; n++) {
      ap  += 1;
      del  = del*x/ap;
      sum += del;
      if (URMath::Abs(del) < Abs(sum*eps)) break;
      //if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
   }
   Double_urt v = sum*Exp(-x+a*Log(x)-gln);
   return v;
}

//______________________________________________________________________________
Double_urt URMath::BreitWigner(Double_urt x, Double_urt mean, Double_urt gamma)
{
   // Calculate a Breit Wigner function with mean and gamma.

   Double_urt bw = gamma/((x-mean)*(x-mean) + gamma*gamma/4);
   return bw/(2*Pi());
}

//______________________________________________________________________________
Double_urt URMath::Gaus(Double_urt x, Double_urt mean, Double_urt sigma, Bool_urt norm)
{
   // Calculate a gaussian function with mean and sigma.
   // if norm=kurTRUE (default is kurFALSE) the result is divided by sqrt(2*Pi)*sigma
   if (sigma == 0) return 1.e30;
   Double_urt arg = (x-mean)/sigma;
   Double_urt res = URMath::Exp(-0.5*arg*arg);
   if (!norm) return res;
   return res/(2.50662827463100024*sigma); //sqrt(2*Pi)=2.50662827463100024
}

//______________________________________________________________________________
Double_urt URMath::Landau(Double_urt x, Double_urt mpv, Double_urt sigma)
{
   // The LANDAU function with mpv(most probable value) and sigma.
   // This function has been adapted from the CERNLIB routine G110 denlan.

   Double_urt p1[5] = {0.4259894875,-0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
   Double_urt q1[5] = {1.0         ,-0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};

   Double_urt p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
   Double_urt q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};

   Double_urt p3[5] = {0.1788544503, 0.09359161662,0.006325387654, 0.00006611667319,-0.000002031049101};
   Double_urt q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};

   Double_urt p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
   Double_urt q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};

   Double_urt p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
   Double_urt q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};

   Double_urt p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
   Double_urt q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};

   Double_urt a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};

   Double_urt a2[2] = {-1.845568670,-4.284640743};

   if (sigma <= 0) return 0;
   Double_urt v = (x-mpv)/sigma;
   Double_urt u, ue, us, den;
   if (v < -5.5) {
      u   = URMath::Exp(v+1.0);
      ue  = URMath::Exp(-1/u);
      us  = URMath::Sqrt(u);
      den = 0.3989422803*(ue/us)*(1+(a1[0]+(a1[1]+a1[2]*u)*u)*u);
   } else if(v < -1) {
       u   = URMath::Exp(-v-1);
       den = URMath::Exp(-u)*URMath::Sqrt(u)*
             (p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v)*v)*v)*v)/
             (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v)*v)*v)*v);
   } else if(v < 1) {
       den = (p2[0]+(p2[1]+(p2[2]+(p2[3]+p2[4]*v)*v)*v)*v)/
             (q2[0]+(q2[1]+(q2[2]+(q2[3]+q2[4]*v)*v)*v)*v);
   } else if(v < 5) {
       den = (p3[0]+(p3[1]+(p3[2]+(p3[3]+p3[4]*v)*v)*v)*v)/
             (q3[0]+(q3[1]+(q3[2]+(q3[3]+q3[4]*v)*v)*v)*v);
   } else if(v < 12) {
       u   = 1/v;
       den = u*u*(p4[0]+(p4[1]+(p4[2]+(p4[3]+p4[4]*u)*u)*u)*u)/
                 (q4[0]+(q4[1]+(q4[2]+(q4[3]+q4[4]*u)*u)*u)*u);
   } else if(v < 50) {
       u   = 1/v;
       den = u*u*(p5[0]+(p5[1]+(p5[2]+(p5[3]+p5[4]*u)*u)*u)*u)/
                 (q5[0]+(q5[1]+(q5[2]+(q5[3]+q5[4]*u)*u)*u)*u);
   } else if(v < 300) {
       u   = 1/v;
       den = u*u*(p6[0]+(p6[1]+(p6[2]+(p6[3]+p6[4]*u)*u)*u)*u)/
                 (q6[0]+(q6[1]+(q6[2]+(q6[3]+q6[4]*u)*u)*u)*u);
   } else {
       u   = 1/(v-v*URMath::Log(v)/(v+1));
       den = u*u*(1+(a2[0]+a2[1]*u)*u);
   }
   return den;
}

//______________________________________________________________________________
Double_urt URMath::LnGamma(Double_urt z)
{
   // Computation of ln[gamma(z)] for all z>0.
   //
   // C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
   //
   // The accuracy of the result is better than 2e-10.
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   if (z<=0) return 0;

   // Coefficients for the series expansion
   Double_urt c[7] = { 2.5066282746310005, 76.18009172947146, -86.50532032941677
                   ,24.01409824083091,  -1.231739572450155, 0.1208650973866179e-2
                   ,-0.5395239384953e-5};

   Double_urt x   = z;
   Double_urt y   = x;
   Double_urt tmp = x+5.5;
   tmp = (x+0.5)*Log(tmp)-tmp;
   Double_urt ser = 1.000000000190015;
   for (Int_urt i=1; i<7; i++) {
      y   += 1;
      ser += c[i]/y;
   }
   Double_urt v = tmp+Log(c[0]*ser/x);
   return v;
}

//______________________________________________________________________________
Float_urt URMath::Normalize(Float_urt v[3])
{
   // Normalize a vector v in place.
   // Returns the norm of the original vector.

   Float_urt d = Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
   if (d != 0) {
      v[0] /= d;
      v[1] /= d;
      v[2] /= d;
   }
   return d;
}

//______________________________________________________________________________
Double_urt URMath::Normalize(Double_urt v[3])
{
   // Normalize a vector v in place.
   // Returns the norm of the original vector.

    Double_urt d = Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (d != 0)
    {
      v[0] /= d;
      v[1] /= d;
      v[2] /= d;
    }
    return d;
}

//______________________________________________________________________________
Float_urt *URMath::Normal2Plane(Float_urt p1[3],Float_urt p2[3],Float_urt p3[3], Float_urt normal[3])
{
   // Calculate a normal vector of a plane.
   //
   //  Input:
   //     Float_urt *p1,*p2,*p3  -  3 3D points belonged the plane to define it.
   //
   //  Return:
   //     Pointer to 3D normal vector (normalized)

   Float_urt v1[3], v2[3];

   v1[0] = p2[0] - p1[0];
   v1[1] = p2[1] - p1[1];
   v1[2] = p2[2] - p1[2];

   v2[0] = p3[0] - p1[0];
   v2[1] = p3[1] - p1[1];
   v2[2] = p3[2] - p1[2];

   NormCross(v1,v2,normal);
   return normal;
}

//______________________________________________________________________________
Double_urt *URMath::Normal2Plane(Double_urt p1[3],Double_urt p2[3],Double_urt p3[3], Double_urt normal[3])
{
   // Calculate a normal vector of a plane.
   //
   //  Input:
   //     Float_urt *p1,*p2,*p3  -  3 3D points belonged the plane to define it.
   //
   //  Return:
   //     Pointer to 3D normal vector (normalized)

   Double_urt v1[3], v2[3];

   v1[0] = p2[0] - p1[0];
   v1[1] = p2[1] - p1[1];
   v1[2] = p2[2] - p1[2];

   v2[0] = p3[0] - p1[0];
   v2[1] = p3[1] - p1[1];
   v2[2] = p3[2] - p1[2];

   NormCross(v1,v2,normal);
   return normal;
}

//______________________________________________________________________________
Double_urt URMath::Poisson(Double_urt x, Double_urt par)
{
  // compute the Poisson distribution function for (x,par)

  if (x<0) return 0;                                                         
  if (x==0) return URMath::Exp(-par);                                                         
  return URMath::Power(par,x)/URMath::Gamma(x+1)/URMath::Exp(par); 
}                                                                              
    
//______________________________________________________________________________
Double_urt URMath::Prob(Double_urt chi2,Int_urt ndf)
{
   // Computation of the probability for a certain Chi-squared (chi2)
   // and number of degrees of freedom (ndf).
   //
   // Calculations are based on the incomplete gamma function P(a,x),
   // where a=ndf/2 and x=chi2/2.
   //
   // P(a,x) represents the probability that the observed Chi-squared
   // for a correct model should be less than the value chi2.
   //
   // The returned probability corresponds to 1-P(a,x),
   // which denotes the probability that an observed Chi-squared exceeds
   // the value chi2 by chance, even for a correct model.
   //
   //--- NvE 14-nov-1998 UU-SAP Utrecht

   if (ndf <= 0) return 0; // Set CL to zero in case ndf<=0

   if (chi2 <= 0) {
      if (chi2 < 0) return 0;
      else          return 1;
   }

   if (ndf==1) {
      Double_urt v = 1.-Erf(Sqrt(chi2)/Sqrt(2.));
      return v;
   }

   // Gaussian approximation for large ndf
   Double_urt q = Sqrt(2*chi2)-Sqrt(Double_urt(2*ndf-1));
   if (ndf > 30 && q > 5) {
      Double_urt v = 0.5*(1-Erf(q/Sqrt(2.)));
      return v;
   }

   // Evaluate the incomplete gamma function
   return (1-Gamma(0.5*ndf,0.5*chi2));
}

//______________________________________________________________________________
Double_urt URMath::KolmogorovProb(Double_urt z)
{
   // Calculates the Kolmogorov distribution function,
   //Begin_Html
   /*
   <img src="gif/kolmogorov.gif">
   */
   //End_Html
   // which gives the probability that Kolmogorov's test statistic will exceed
   // the value z assuming the null hypothesis. This gives a very powerful
   // test for comparing two one-dimensional distributions.
   // see, for example, Eadie et al, "statistocal Methods in Experimental
   // Physics', pp 269-270).
   //
   // This function returns the confidence level for the null hypothesis, where:
   //   z = dn*sqrt(n), and
   //   dn  is the maximum deviation between a hypothetical distribution
   //       function and an experimental distribution with
   //   n    events
   //
   // NOTE: To compare two experimental distributions with m and n events,
   //       use z = sqrt(m*n/(m+n))*dn
   //
   // Accuracy: The function is far too accurate for any imaginable application.
   //           Probabilities less than 10^-15 are returned as zero.
   //           However, remember that the formula is only valid for "large" n.
   // Theta function inversion formula is used for z <= 1
   //
   // This function was translated by Rene Brun from PROBKL in CERNLIB.

   Double_urt fj[4] = {-2,-8,-18,-32}, r[4];
   const Double_urt w = 2.50662827;
   // c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
   const Double_urt c1 = -1.2337005501361697;
   const Double_urt c2 = -11.103304951225528;
   const Double_urt c3 = -30.842513753404244;
   
   Double_urt u = URMath::Abs(z);
   Double_urt p;
   if (u < 0.2) {
      p = 1;
   } else if (u < 0.755) {
      Double_urt v = 1./(u*u);
      p = 1 - w*(URMath::Exp(c1*v) + URMath::Exp(c2*v) + URMath::Exp(c3*v))/u;
   } else if (u < 6.8116) {
      r[1] = 0;
      r[2] = 0;
      r[3] = 0;
      Double_urt v = u*u;
      Int_urt maxj = URMath::Max(1,URMath::Nint(3./u));
      for (Int_urt j=0; j<maxj;j++) {
         r[j] = URMath::Exp(fj[j]*v);
      }
      p = 2*(r[0] - r[1] +r[2] - r[3]);
   } else {
      p = 0;
   }
   return p;     
   }

//______________________________________________________________________________
Double_urt URMath::Voigt(Double_urt x, Double_urt sigma, Double_urt lg, Int_urt R)
{
   // Computation of Voigt function (normalised).
   // Voigt is a convolution of
   // gauss(x) = 1/(sqrt(2*pi)*sigma) * exp(x*x/(2*sigma*sigma)
   // and
   // lorentz(x) = (1/pi) * (lg/2) / (x*x + g*g/4)
   // functions.
   //
   // The Voigt function is known to be the real part of Faddeeva function also
   // called complex error function [2].
   //
   // The algoritm was developed by J. Humlicek [1].
   // This code is based on fortran code presented by R. J. Wells [2].
   // Translated and adapted by Miha D. Puc
   //
   // To calculate the Faddeeva function with relative error less than 10^(-R).
   // R can be set by the the user subject to the constraints 2 <= R <= 5.
   //
   // [1] J. Humlicek, JQSRT, 21, 437 (1982).
   // [2] R.J. Wells "Rapid Approximation to the Voigt/Faddeeva Function and its
   // Derivatives" JQSRT 62 (1999), pp 29-48.
   // http://www-atm.physics.ox.ac.uk/user/wells/voigt.html

   if ((sigma < 0 || lg < 0) || (sigma==0 && lg==0)) {
      return 0;  // Not meant to be for those who want to be thinner than 0
   }

   if (sigma == 0) {
      return lg * 0.159154943  / (x*x + lg*lg /4); //pure Lorentz
   }

   if (lg == 0) {   //pure gauss
      return 0.39894228 / sigma * URMath::Exp(-x*x / (2*sigma*sigma));
   }

   Double_urt X, Y, K;
   X = x / sigma / 1.41421356;
   Y = lg / 2 / sigma / 1.41421356;

   Double_urt R0, R1;

   if (R < 2) R = 2;
   if (R > 5) R = 5;

   R0=1.51 * exp(1.144 * (Double_urt)R);
   R1=1.60 * exp(0.554 * (Double_urt)R);

   // Constants

   const Double_urt RRTPI = 0.56418958;  // 1/SQRT(pi)

   Double_urt Y0, Y0PY0, Y0Q;                      // for CPF12 algorithm
   Y0 = 1.5;
   Y0PY0 = Y0 + Y0;
   Y0Q = Y0 * Y0;

   Double_urt C[6] = { 1.0117281, -0.75197147, 0.012557727, 0.010022008, -0.00024206814, 0.00000050084806};
   Double_urt S[6] = { 1.393237, 0.23115241, -0.15535147, 0.0062183662, 0.000091908299, -0.00000062752596};
   Double_urt T[6] = { 0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249};

   // Local variables

   int J;					 // Loop variables
   int RG1, RG2, RG3; 				 // y polynomial flags
   Double_urt ABX, XQ, YQ, YRRTPI; 		 // --x--, x^2, y^2, y/SQRT(pi)
   Double_urt XLIM0, XLIM1, XLIM2, XLIM3, XLIM4; 	 // --x-- on region boundaries
   Double_urt A0=0, D0=0, D2=0, E0=0, E2=0, E4=0, H0=0, H2=0, H4=0, H6=0;// W4 temporary variables
   Double_urt P0=0, P2=0, P4=0, P6=0, P8=0, Z0=0, Z2=0, Z4=0, Z6=0, Z8=0;
   Double_urt XP[6], XM[6], YP[6], YM[6];          // CPF12 temporary values
   Double_urt MQ[6], PQ[6], MF[6], PF[6];
   Double_urt D, YF, YPY0, YPY0Q;

   //***** Start of executable code *****************************************

   RG1 = 1;  // Set flags
   RG2 = 1;
   RG3 = 1;
   YQ = Y * Y;  // y^2
   YRRTPI = Y * RRTPI;  // y/SQRT(pi)

   // Region boundaries when both K and L are required or when R<>4

   XLIM0 = R0 - Y;
   XLIM1 = R1 - Y;
   XLIM3 = 3.097 * Y - 0.45;
   XLIM2 = 6.8 - Y;
   XLIM4 = 18.1 * Y + 1.65;
   if ( Y <= 1e-6 ) { 	                   // When y<10^-6 avoid W4 algorithm.
      XLIM1 = XLIM0; 
      XLIM2 = XLIM0;
   }

   ABX = fabs(X); 				 // |x|
   XQ = ABX * ABX;			         // x^2
   if ( ABX > XLIM0 ) { 			 // Region 0 algorithm
      K = YRRTPI / (XQ + YQ);
   } else if ( ABX > XLIM1 ) { 			 // Humlicek W4 Region 1
      if ( RG1 != 0 ) { 			 // First point in Region 1
	 RG1 = 0;
	 A0 = YQ + 0.5; 			 // Region 1 y-dependents
	 D0 = A0*A0;
	 D2 = YQ + YQ - 1.0;
      }
      D = RRTPI / (D0 + XQ*(D2 + XQ));
      K = D * Y * (A0 + XQ);
   } else if ( ABX > XLIM2 ) { 			 // Humlicek W4 Region 2
      if ( RG2 != 0 ) { 			 // First point in Region 2
	 RG2 = 0;
	 H0 = 0.5625 + YQ * (4.5 + YQ * (10.5 + YQ * (6.0 + YQ)));
	                                         // Region 2 y-dependents
	 H2 = -4.5 + YQ * (9.0 + YQ * ( 6.0 + YQ * 4.0));
	 H4 = 10.5 - YQ * (6.0 - YQ * 6.0);
	 H6 = -6.0 + YQ * 4.0;
	 E0 = 1.875 + YQ * (8.25 + YQ * (5.5 + YQ));
	 E2 = 5.25 + YQ * (1.0 + YQ * 3.0);
	 E4 = 0.75 * H6;
      }
      D = RRTPI / (H0 + XQ * (H2 + XQ * (H4 + XQ * (H6 + XQ))));
      K = D * Y * (E0 + XQ * (E2 + XQ * (E4 + XQ)));
   } else if ( ABX < XLIM3 ) { 			 // Humlicek W4 Region 3
      if ( RG3 != 0 ) { 			 // First point in Region 3
	 RG3 = 0;
	 Z0 = 272.1014 + Y * (1280.829 + Y *
			      (2802.870 + Y *
			       (3764.966 + Y *
				(3447.629 + Y *
				 (2256.981 + Y *
				  (1074.409 + Y *
				   (369.1989  + Y *
				    (88.26741 + Y *
				     (13.39880 + Y)
				     ))))))));   // Region 3 y-dependents
	 Z2 = 211.678 + Y * (902.3066 + Y *
			     (1758.336 + Y *
			      (2037.310 + Y *
			       (1549.675 + Y *
				(793.4273 + Y *
				 (266.2987 + Y *
				  (53.59518 + Y * 5.0)
				  ))))));
	 Z4 = 78.86585 + Y * (308.1852 + Y *
			      (497.3014 + Y *
			       (479.2576 + Y *
				(269.2916 + Y *
				 (80.39278 + Y * 10.0)
				 ))));
	 Z6 = 22.03523 + Y * (55.02933 + Y *
			      (92.75679 + Y *
			       (53.59518 + Y * 10.0)
			       ));
	 Z8 = 1.496460 + Y * (13.39880 + Y * 5.0);
	 P0 = 153.5168 + Y * (549.3954 + Y *
			      (919.4955 + Y *
			       (946.8970 + Y *
				(662.8097 + Y *
				 (328.2151 + Y *
				  (115.3772 + Y *
				   (27.93941 + Y *
				    (4.264678 + Y * 0.3183291)
				    )))))));
	 P2 = -34.16955 + Y * (-1.322256+ Y *
			       (124.5975 + Y *
				(189.7730 + Y *
				 (139.4665 + Y *
				  (56.81652 + Y *
				   (12.79458 + Y * 1.2733163)
				   )))));
	 P4 = 2.584042 + Y * (10.46332 + Y *
			      (24.01655 + Y *
			       (29.81482 + Y *
				(12.79568 + Y * 1.9099744)
				)));
	 P6 = -0.07272979 + Y * (0.9377051 + Y *
				 (4.266322 + Y * 1.273316));
	 P8 = 0.0005480304 + Y * 0.3183291;
      }
      D = 1.7724538 / (Z0 + XQ * (Z2 + XQ * (Z4 + XQ * (Z6 + XQ * (Z8 + XQ)))));
      K = D * (P0 + XQ * (P2 + XQ * (P4 + XQ * (P6 + XQ * P8))));
   } else { 					 // Humlicek CPF12 algorithm
      YPY0 = Y + Y0;
      YPY0Q = YPY0 * YPY0;
      K = 0.0;
      for (J = 0; J <= 5; J++) {
	 D = X - T[J];
	 MQ[J] = D * D;
	 MF[J] = 1.0 / (MQ[J] + YPY0Q);
	 XM[J] = MF[J] * D;
	 YM[J] = MF[J] * YPY0;
	 D = X + T[J];
	 PQ[J] = D * D;
	 PF[J] = 1.0 / (PQ[J] + YPY0Q);
	 XP[J] = PF[J] * D;
	 YP[J] = PF[J] * YPY0;
      }
      if ( ABX <= XLIM4 ) {           		 // Humlicek CPF12 Region I
	 for (J = 0; J <= 5; J++) {
	    K = K + C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]) ;
	 }
      } else { 					 // Humlicek CPF12 Region II
	 YF = Y + Y0PY0;
	 for ( J = 0; J <= 5; J++) {
	    K = K + (C[J] *
		     (MQ[J] * MF[J] - Y0 * YM[J])
		     + S[J] * YF * XM[J]) / (MQ[J]+Y0Q)
	       + (C[J] * (PQ[J] * PF[J] - Y0 * YP[J])
		  - S[J] * YF * XP[J]) / (PQ[J]+Y0Q);
	 }
	 K = Y * K + exp( -XQ );
      }
   }
   return K / 2.506628 / sigma; // Normalize by dividing by sqrt(2*pi)*sigma.
}

//______________________________________________________________________________
Int_urt URMath::LocMin(Int_urt n, const Short_urt *a)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;
  Short_urt xmin = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMin(Int_urt n, const Int_urt *a)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;
  Int_urt xmin = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMin(Int_urt n, const Float_urt *a)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;
  Float_urt xmin = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMin(Int_urt n, const Double_urt *a)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;
  Double_urt xmin = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMin(Int_urt n, const Long_urt *a)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;
  Long_urt xmin = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMax(Int_urt n, const Short_urt *a)
{
   // Return index of array with the maximum element.
   // If more than one element is maximum returns first found.

  if  (n <= 0) return -1;
  Short_urt xmax = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMax(Int_urt n, const Int_urt *a)
{
   // Return index of array with the maximum element.
   // If more than one element is maximum returns first found.

  if  (n <= 0) return -1;
  Int_urt xmax = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMax(Int_urt n, const Float_urt *a)
{
   // Return index of array with the maximum element.
   // If more than one element is maximum returns first found.

  if  (n <= 0) return -1;
  Float_urt xmax = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMax(Int_urt n, const Double_urt *a)
{
   // Return index of array with the maximum element.
   // If more than one element is maximum returns first found.

  if  (n <= 0) return -1;
  Double_urt xmax = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::LocMax(Int_urt n, const Long_urt *a)
{
   // Return index of array with the maximum element.
   // If more than one element is maximum returns first found.

  if  (n <= 0) return -1;
  Long_urt xmax = a[0];
  Int_urt loc = 0;
  for  (Int_urt i = 1; i < n; i++) {
     if (xmax < a[i])  {
         xmax = a[i];
         loc = i;
     }
  }
  return loc;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Short_urt *array, Short_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == array[middle-1]) return middle-1;
      if (value  < array[middle-1]) nabove = middle;
      else                          nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Short_urt **array, Short_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == *array[middle-1]) return middle-1;
      if (value  < *array[middle-1]) nabove = middle;
      else                           nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Int_urt *array, Int_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == array[middle-1]) return middle-1;
      if (value  < array[middle-1]) nabove = middle;
      else                          nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Int_urt **array, Int_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == *array[middle-1]) return middle-1;
      if (value  < *array[middle-1]) nabove = middle;
      else                           nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Float_urt *array, Float_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == array[middle-1]) return middle-1;
      if (value  < array[middle-1]) nabove = middle;
      else                          nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Float_urt **array, Float_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == *array[middle-1]) return middle-1;
      if (value  < *array[middle-1]) nabove = middle;
      else                           nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Double_urt *array, Double_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == array[middle-1]) return middle-1;
      if (value  < array[middle-1]) nabove = middle;
      else                          nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Double_urt **array, Double_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == *array[middle-1]) return middle-1;
      if (value  < *array[middle-1]) nabove = middle;
      else                           nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Long_urt *array, Long_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == array[middle-1]) return middle-1;
      if (value  < array[middle-1]) nabove = middle;
      else                          nbelow = middle;
   }
   return nbelow-1;
}

//______________________________________________________________________________
Int_urt URMath::BinarySearch(Int_urt n, const Long_urt **array, Long_urt value)
{
   // Binary search in an array of n values to locate value.
   //
   // Array is supposed  to be sorted prior to this call.
   // If match is found, function returns position of element.
   // If no match found, function gives nearest element smaller than value.

   Int_urt nabove, nbelow, middle;
   nabove = n+1;
   nbelow = 0;
   while(nabove-nbelow > 1) {
      middle = (nabove+nbelow)/2;
      if (value == *array[middle-1]) return middle-1;
      if (value  < *array[middle-1]) nabove = middle;
      else                           nbelow = middle;
   }
   return nbelow-1;
}

//_____________________________________________________________________________
Bool_urt URMath::IsInside(Double_urt xp, Double_urt yp, Int_urt np, Double_urt *x, Double_urt *y)
{
   // Function which returns kurTRUE if point xp,yp lies inside the
   // polygon defined by the np points in arrays x and y, kurFALSE otherwise
   // NOTE that the polygon must be a closed polygon (1st and last point
   // must be identical)
   Double_urt xint;
   Int_urt i;
   Int_urt inter = 0;
   for (i=0;i<np-1;i++) {
      if (y[i] == y[i+1]) continue;
      if (yp <= y[i] && yp <= y[i+1]) continue;
      if (y[i] < yp && y[i+1] < yp) continue;
      xint = x[i] + (yp-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]);
      if (xp < xint) inter++;
   }
   if (inter%2) return kurTRUE;
   return kurFALSE;
}

//_____________________________________________________________________________
Bool_urt URMath::IsInside(Float_urt xp, Float_urt yp, Int_urt np, Float_urt *x, Float_urt *y)
{
   // Function which returns kurTRUE if point xp,yp lies inside the
   // polygon defined by the np points in arrays x and y, kurFALSE otherwise
   // NOTE that the polygon must be a closed polygon (1st and last point
   // must be identical)
   Double_urt xint;
   Int_urt i;
   Int_urt inter = 0;
   for (i=0;i<np-1;i++) {
      if (y[i] == y[i+1]) continue;
      if (yp <= y[i] && yp <= y[i+1]) continue;
      if (y[i] < yp && y[i+1] < yp) continue;
      xint = x[i] + (yp-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]);
      if ((Double_urt)xp < xint) inter++;
   }
   if (inter%2) return kurTRUE;
   return kurFALSE;
}

//_____________________________________________________________________________
Bool_urt URMath::IsInside(Int_urt xp, Int_urt yp, Int_urt np, Int_urt *x, Int_urt *y)
{
   // Function which returns kurTRUE if point xp,yp lies inside the
   // polygon defined by the np points in arrays x and y, kurFALSE otherwise
   // NOTE that the polygon must be a closed polygon (1st and last point
   // must be identical)
   Double_urt xint;
   Int_urt i;
   Int_urt inter = 0;
   for (i=0;i<np-1;i++) {
      if (y[i] == y[i+1]) continue;
      if (yp <= y[i] && yp <= y[i+1]) continue;
      if (y[i] < yp && y[i+1] < yp) continue;
      xint = x[i] + (yp-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]);
      if ((Double_urt)xp < xint) inter++;
   }
   if (inter%2) return kurTRUE;
   return kurFALSE;
}

//_____________________________________________________________________________
void URMath::Sort(Int_urt n1, const Short_urt *a, Int_urt *index, Bool_urt down)
{
   // Sort the n1 elements of the Short_urt array a.
   // In output the array index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).
   // This is a translation of the CERNLIB routine sortzv (M101)
   // based on the quicksort algorithm.
   // NOTE that the array index must be created with a length >= n1
   // before calling this function.

   Int_urt i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_urt i22 = 0;
   Short_urt ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   for (i=0;i<n1;i++) index[i]--;
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}

//_____________________________________________________________________________
void URMath::Sort(Int_urt n1, const Int_urt *a, Int_urt *index, Bool_urt down)
{
   // Sort the n1 elements of the Int_urt array a.
   // In output the array index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).
   // This is a translation of the CERNLIB routine sortzv (M101)
   // based on the quicksort algorithm.
   // NOTE that the array index must be created with a length >= n1
   // before calling this function.

   Int_urt i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_urt i22 = 0;
   Int_urt ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   for (i=0;i<n1;i++) index[i]--;
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}

//_____________________________________________________________________________
void URMath::Sort(Int_urt n1, const Float_urt *a, Int_urt *index, Bool_urt down)
{
   // Sort the n1 elements of the Float_urt array a.
   // In output the array index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).
   // This is a translation of the CERNLIB routine sortzv (M101)
   // based on the quicksort algorithm.
   // NOTE that the array index must be created with a length >= n1
   // before calling this function.

   Int_urt i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_urt i22 = 0;
   Float_urt ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   for (i=0;i<n1;i++) index[i]--;
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}

//_____________________________________________________________________________
void URMath::Sort(Int_urt n1, const Double_urt *a, Int_urt *index, Bool_urt down)
{
   // Sort the n1 elements of the Double_urt array a.
   // In output the array index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).
   // This is a translation of the CERNLIB routine sortzv (M101)
   // based on the quicksort algorithm.
   // NOTE that the array index must be created with a length >= n1
   // before calling this function.

   Int_urt i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_urt i22 = 0;
   Double_urt ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   for (i=0;i<n1;i++) index[i]--;
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}

//_____________________________________________________________________________
void URMath::Sort(Int_urt n1, const Long_urt *a, Int_urt *index, Bool_urt down)
{
   // Sort the n1 elements of the Long_urt array a.
   // In output the array index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).
   // This is a translation of the CERNLIB routine sortzv (M101)
   // based on the quicksort algorithm.
   // NOTE that the array index must be created with a length >= n1
   // before calling this function.

   Int_urt i,i1,n,i2,i3,i33,i222,iswap,n2;
   Int_urt i22 = 0;
   Long_urt ai;
   n = n1;
   if (n <= 0) return;
   if (n == 1) {index[0] = 0; return;}
   for (i=0;i<n;i++) index[i] = i+1;
   for (i1=2;i1<=n;i1++) {
      i3 = i1;
      i33 = index[i3-1];
      ai  = a[i33-1];
      while(1) {
         i2 = i3/2;
         if (i2 <= 0) break;
         i22 = index[i2-1];
         if (ai <= a[i22-1]) break;
         index[i3-1] = i22;
         i3 = i2;
      }
      index[i3-1] = i33;
   }

   while(1) {
      i3 = index[n-1];
      index[n-1] = index[0];
      ai = a[i3-1];
      n--;
      if(n-1 < 0) {index[0] = i3; break;}
      i1 = 1;
      while(2) {
         i2 = i1+i1;
         if (i2 <= n) i22 = index[i2-1];
         if (i2-n > 0) {index[i1-1] = i3; break;}
         if (i2-n < 0) {
            i222 = index[i2];
            if (a[i22-1] - a[i222-1] < 0) {
                i2++;
                i22 = i222;
            }
         }
         if (ai - a[i22-1] > 0) {index[i1-1] = i3; break;}
         index[i1-1] = i22;
         i1 = i2;
      }
   }
   for (i=0;i<n1;i++) index[i]--;
   if (!down) return;
   n2 = n1/2;
   for (i=0;i<n2;i++) {
      iswap         = index[i];
      index[i]      = index[n1-i-1];
      index[n1-i-1] = iswap;
   }
}


//______________________________________________________________________________
void URMath::BubbleHigh(Int_urt Narr, Double_urt *arr1, Int_urt *arr2)
{
   // Bubble sort variant to obtain the order of an array's elements into
   // an index in order to do more useful things than the standard built
   // in functions.
   // *arr1 is unchanged;
   // *arr2 is the array of indicies corresponding to the decending value
   // of arr1 with arr2[0] corresponding to the largest arr1 value and
   // arr2[Narr] the smallest.
   //
   //  Author:        Adrian Bevan (bevan@slac.stanford.edu)
   //  Copyright:     Liverpool University, July 2001

   if (Narr <= 0) return;
   double *localArr1 = new double[Narr];
   int    *localArr2 = new int[Narr];
   int iEl;
   int iEl2;

   for(iEl = 0; iEl < Narr; iEl++) {
      localArr1[iEl] = arr1[iEl];
      localArr2[iEl] = iEl;
   }

   for (iEl = 0; iEl < Narr; iEl++) {
      for (iEl2 = Narr-1; iEl2 > iEl; --iEl2) {
         if (localArr1[iEl2-1] < localArr1[iEl2]) {
            double tmp        = localArr1[iEl2-1];
            localArr1[iEl2-1] = localArr1[iEl2];
            localArr1[iEl2]   = tmp;

            int    tmp2       = localArr2[iEl2-1];
            localArr2[iEl2-1] = localArr2[iEl2];
            localArr2[iEl2]   = tmp2;
         }
      }
   }

   for (iEl = 0; iEl < Narr; iEl++) {
      arr2[iEl] = localArr2[iEl];
   }
   delete [] localArr2;
   delete [] localArr1;
}

//______________________________________________________________________________
void URMath::BubbleLow(Int_urt Narr, Double_urt *arr1, Int_urt *arr2)
{
   // Opposite ordering of the array arr2[] to that of BubbleHigh.
   //
   //  Author:        Adrian Bevan (bevan@slac.stanford.edu)
   //  Copyright:     Liverpool University, July 2001

   if (Narr <= 0) return;
   double *localArr1 = new double[Narr];
   int    *localArr2 = new int[Narr];
   int iEl;
   int iEl2;

   for (iEl = 0; iEl < Narr; iEl++) {
      localArr1[iEl] = arr1[iEl];
      localArr2[iEl] = iEl;
   }

   for (iEl = 0; iEl < Narr; iEl++) {
      for (iEl2 = Narr-1; iEl2 > iEl; --iEl2) {
         if (localArr1[iEl2-1] > localArr1[iEl2]) {
            double tmp        = localArr1[iEl2-1];
            localArr1[iEl2-1] = localArr1[iEl2];
            localArr1[iEl2]   = tmp;

            int    tmp2       = localArr2[iEl2-1];
            localArr2[iEl2-1] = localArr2[iEl2];
            localArr2[iEl2]   = tmp2;
         }
      }
   }

   for (iEl = 0; iEl < Narr; iEl++) {
      arr2[iEl] = localArr2[iEl];
   }
   delete [] localArr2;
   delete [] localArr1;
}

#ifdef OLD_HASH

//______________________________________________________________________________
ULong_urt URMath::Hash(const void *txt, Int_urt ntxt)
{
   // Calculates hash index from any char string.
   // Based on precalculated table of 256 specially selected random numbers.
   //
   //   For string:  i = URMath::Hash(string,nstring);
   //   For int:     i = URMath::Hash(&intword,sizeof(int));
   //   For pointer: i = URMath::Hash(&pointer,sizeof(void*));
   //
   //   Limitation: for ntxt>256 calculates hash only from first 256 bytes
   //
   //              V.Perev

   const UChar_urt *uc = (const UChar_urt*) txt;
   ULong_urt  u = 0, uu = 0;

   static ULong_urt utab[256] =
      {0xb93f6fc0,0x553dfc51,0xb22c1e8c,0x638462c0,0x13e81418,0x2836e171,0x7c4abb90,0xda1a4f39
      ,0x38f211d1,0x8c804829,0x95a6602d,0x4c590993,0x1810580a,0x721057d4,0x0f587215,0x9f49ce2a
      ,0xcd5ab255,0xab923a99,0x80890f39,0xbcfa2290,0x16587b52,0x6b6d3f0d,0xea8ff307,0x51542d5c
      ,0x189bf223,0x39643037,0x0e4a326a,0x214eca01,0x47645a9b,0x0f364260,0x8e9b2da4,0x5563ebd9
      ,0x57a31c1c,0xab365854,0xdd63ab1f,0x0b89acbd,0x23d57d33,0x1800a0fd,0x225ac60a,0xd0e51943
      ,0x6c65f669,0xcb966ea0,0xcbafda95,0x2e5c0c5f,0x2988e87e,0xc781cbab,0x3add3dc7,0x693a2c30
      ,0x42d6c23c,0xebf85f26,0x2544987e,0x2e315e3f,0xac88b5b5,0x7ebd2bbb,0xda07c87b,0x20d460f1
      ,0xc61c3f40,0x182046e7,0x3b6c3b66,0x2fc10d4a,0x0780dfbb,0xc437280c,0x0988dd07,0xe1498606
      ,0x8e61d728,0x4f1f3909,0x040a9682,0x49411b29,0x391b0e1c,0xd7905241,0xdd77d95b,0x88426c13
      ,0x33033e58,0xe158e30e,0x7e342647,0x1e09544b,0x4637353d,0x18ea0924,0x39212b08,0x12580ae8
      ,0x269a6f06,0x3e10b73b,0x123db33b,0x085412da,0x3bb5f464,0xd9b2d442,0x103d26bb,0xd0038bab
      ,0x45b6177f,0xfb48f1fe,0x99074c55,0xb545e82e,0x5f79fd0d,0x570f3ae4,0x57e43255,0x037a12ae
      ,0x4357bdb2,0x337c2c4d,0x2982499d,0x2ab72793,0x3900e7d1,0x57a6bb81,0x7503609b,0x3f39c0d0
      ,0x717b389d,0x5748034f,0x4698162b,0x5801b97c,0x1dfd5d7e,0xc1386d1c,0xa387a72a,0x084547e4
      ,0x2e54d8e9,0x2e2f384c,0xe09ccc20,0x8904b71e,0x3e24edc5,0x06a22e16,0x8a2be1df,0x9e5058b2
      ,0xe01a2f16,0x03325eed,0x587ecfe6,0x584d9cd3,0x32926930,0xe943d68c,0xa9442da8,0xf9650560
      ,0xf003871e,0x1109c663,0x7a2f2f89,0x1c2210bb,0x37335787,0xb92b382f,0xea605cb5,0x336bbe38
      ,0x08126bd3,0x1f8c2bd6,0xba6c46f2,0x1a4d1b83,0xc988180d,0xe2582505,0xa8a1b375,0x59a08c49
      ,0x3db54b48,0x44400f35,0x272d4e7f,0x5579f733,0x98eb590e,0x8ee09813,0x12cc9301,0xc85c402d
      ,0x135c1039,0x22318128,0x4063c705,0x87a8a3fa,0xfc14431f,0x6e27bf47,0x2d080a19,0x01dba174
      ,0xe343530b,0xaa1bfced,0x283bb2c8,0x5df250c8,0x4ff9140b,0x045039c1,0xa377780d,0x750f2661
      ,0x2b108918,0x0b152120,0x3cbc251f,0x5e87b350,0x060625bb,0xe068ba3b,0xdb73ebd7,0x66014ff3
      ,0xdb003000,0x161a3a0b,0xdc24e142,0x97ea5575,0x635a3cab,0xa719100a,0x256084db,0xc1f4a1e7
      ,0xe13388f2,0xb8199fc9,0x50c70dc9,0x08154211,0xd60e5220,0xe52c6592,0x584c5fe1,0xfe5e0875
      ,0x21072b30,0x3370d773,0x92608fe2,0x2d013d93,0x53414b3c,0x2c066142,0x64676644,0x0420887c
      ,0x35c01187,0x6822119b,0xf9bfe6df,0x273f4ee4,0x87973149,0x7b41282d,0x635d0d1f,0x5f7ecc1e
      ,0x14c3608a,0x462dfdab,0xc33d8808,0x1dcd995e,0x0fcb11ba,0x11755914,0x5a62044b,0x37f76755
      ,0x345bd058,0x8831c2b5,0x204a8468,0x3b0b1cd2,0x444e56f4,0x97a93e2c,0xd5f15067,0x266a95fa
      ,0xff4f8036,0x6160060d,0x930c472f,0xed922184,0x37120251,0xc0add74f,0x1c0bc89d,0x018d47f2
      ,0xff59ef66,0xd1901a17,0x91f6701b,0x0960082f,0x86f6a8f3,0x1154fecd,0x9867d1de,0x0945482f
      ,0x790ffcac,0xe5610011,0x4765637e,0xa745dbff,0x841fdcb3,0x4f7372a0,0x3c05013d,0xf1ac4ab7
      ,0x3bc5b5cc,0x49a73349,0x356a7f67,0x1174f031,0x11d32634,0x4413d301,0x1dd285c4,0x3fae4800
   };

   if (ntxt > 255) ntxt = 255;

   for ( ; ntxt--; uc++) {
      uu = uu<<1 ^ utab[(*uc) ^ ntxt];
      u ^= uu;
   }
   return u;
}

#else

//______________________________________________________________________________
ULong_urt URMath::Hash(const void *txt, Int_urt ntxt)
{
   // Calculates hash index from any char string.
   // Based on precalculated table of 256 specially selected numbers.
   // These numbers are selected in such a way, that for string
   // length == 4 (integer number) the hash is unambigous, i.e.
   // from hash value we can recalculate input (no degeneration).
   //
   // The quality of hash method is good enough, that
   // "random" numbers made as R = Hash(1), Hash(2), ...Hash(N)
   // tested by <R>, <R*R>, <Ri*Ri+1> gives the same result
   // as for libc rand().
   //
   // For string:  i = URMath::Hash(string,nstring);
   // For int:     i = URMath::Hash(&intword,sizeof(int));
   // For pointer: i = URMath::Hash(&pointer,sizeof(void*));
   //
   //              V.Perev

   static const ULong_urt utab[] = {
       0xdd367647,0x9caf993f,0x3f3cc5ff,0xfde25082,0x4c764b21,0x89affca7,0x5431965c,0xce22eeec
      ,0xc61ab4dc,0x59cc93bd,0xed3107e3,0x0b0a287a,0x4712475a,0xce4a4c71,0x352c8403,0x94cb3cee
      ,0xc3ac509b,0x09f827a2,0xce02e37e,0x7b20bbba,0x76adcedc,0x18c52663,0x19f74103,0x6f30e47b
      ,0x132ea5a1,0xfdd279e0,0xa3d57d00,0xcff9cb40,0x9617f384,0x6411acfa,0xff908678,0x5c796b2c
      ,0x4471b62d,0xd38e3275,0xdb57912d,0x26bf953f,0xfc41b2a5,0xe64bcebd,0x190b7839,0x7e8e6a56
      ,0x9ca22311,0xef28aa60,0xe6b9208e,0xd257fb65,0x45781c2c,0x9a558ac3,0x2743e74d,0x839417a8
      ,0x06b54d5d,0x1a82bcb4,0x06e97a66,0x70abdd03,0xd163f30d,0x222ed322,0x777bfeda,0xab7a2e83
      ,0x8494e0cf,0x2dca2d4f,0x78f94278,0x33f04a09,0x402b6452,0x0cd8b709,0xdb72a39e,0x170e00a2
      ,0x26354faa,0x80e57453,0xcfe8d4e1,0x19e45254,0x04c291c3,0xeb503738,0x425af3bc,0x67836f2a
      ,0xfac22add,0xfafc2b8c,0x59b8c2a0,0x03e806f9,0xcb4938b9,0xccc942af,0xcee3ae2e,0xfbe748fa
      ,0xb223a075,0x85c49b5d,0xe4576ac9,0x0fbd46e2,0xb49f9cf5,0xf3e1e86a,0x7d7927fb,0x711afe12
      ,0xbf61c346,0x157c9956,0x86b6b046,0x2e402146,0xb2a57d8a,0x0d064bb1,0x30ce390c,0x3a3e1eb1
      ,0xbe7f6f8f,0xd8e30f87,0x5be2813c,0x73a3a901,0xa3aaf967,0x59ff092c,0x1705c798,0xf610dd66
      ,0xb17da91e,0x8e59534e,0x2211ea5b,0xa804ba03,0xd890efbb,0xb8b48110,0xff390068,0xc8c325b4
      ,0xf7289c07,0x787e104f,0x3d0df3d0,0x3526796d,0x10548055,0x1d59a42b,0xed1cc5a3,0xdd45372a
      ,0x31c50d57,0x65757cb7,0x3cfb85be,0xa329910d,0x6ad8ce39,0xa2de44de,0x0dd32432,0xd4a5b617
      ,0x8f3107fc,0x96485175,0x7f94d4f3,0x35097634,0xdb3ca782,0x2c0290b8,0x2045300b,0xe0f5d15a
      ,0x0e8cbffa,0xaa1cc38a,0x84008d6f,0xe9a9e794,0x5c602c25,0xfa3658fa,0x98d9d82b,0x3f1497e7
      ,0x84b6f031,0xe381eff9,0xfc7ae252,0xb239e05d,0xe3723d1f,0xcc3bda82,0xe21b1ad3,0x9104f7c8
      ,0x4bb2dfcd,0x4d14a8bc,0x6ba7f28c,0x8f89886c,0xad44c97e,0xb30fd975,0x633cdab1,0xf6c2d514
      ,0x067a49d2,0xdc461ad9,0xebaf9f3f,0x8dc6cac3,0x7a060f16,0xbab063ad,0xf42e25e6,0x60724ca6
      ,0xc7245c2e,0x4e48ea3c,0x9f89a609,0xa1c49890,0x4bb7f116,0xd722865c,0xa8ee3995,0x0ee070b1
      ,0xd9bffcc2,0xe55b64f9,0x25507a5a,0xc7a3e2b5,0x5f395f7e,0xe7957652,0x7381ba6a,0xde3d21f1
      ,0xdf1708dd,0xad0c9d0c,0x00cbc9e5,0x1160e833,0x6779582c,0x29d5d393,0x3f11d7d7,0x826a6b9b
      ,0xe73ff12f,0x8bad3d86,0xee41d3e5,0x7f0c8917,0x8089ef24,0x90c5cb28,0x2f7f8e6b,0x6966418a
      ,0x345453fb,0x7a2f8a68,0xf198593d,0xc079a532,0xc1971e81,0x1ab74e26,0x329ef347,0x7423d3d0
      ,0x942c510b,0x7f6c6382,0x14ae6acc,0x64b59da7,0x2356fa47,0xb6749d9c,0x499de1bb,0x92ffd191
      ,0xe8f2fb75,0x848dc913,0x3e8727d3,0x1dcffe61,0xb6e45245,0x49055738,0x827a6b55,0xb4788887
      ,0x7e680125,0xd19ce7ed,0x6b4b8e30,0xa8cadea2,0x216035d8,0x1c63bc3c,0xe1299056,0x1ad3dff4
      ,0x0aefd13c,0x0e7b921c,0xca0173c6,0x9995782d,0xcccfd494,0xd4b0ac88,0x53d552b1,0x630dae8b
      ,0xa8332dad,0x7139d9a2,0x5d76f2c4,0x7a4f8f1e,0x8d1aef97,0xd1cf285d,0xc8239153,0xce2608a9
      ,0x7b562475,0xe4b4bc83,0xf3db0c3a,0x70a65e48,0x6016b302,0xdebd5046,0x707e786a,0x6f10200c
   };

   static const ULong_urt msk[] = { 0x11111111, 0x33333333, 0x77777777, 0xffffffff };

   const UChar_urt *uc = (const UChar_urt *) txt;
   ULong_urt uu = 0;
   union {
      ULong_urt  u;
      UShort_urt s[2];
   } u;
   u.u = 0;
   Int_urt i, idx;

   for (i = 0; i < ntxt; i++) {
      idx = (uc[i] ^ i) & 255;
      uu  = (uu << 1) ^ (utab[idx] & msk[i & 3]);
      if (i & 3 == 3) u.u ^= uu;
   }
   if (i & 3) u.u ^= uu;

   u.u *= 1879048201;      // prime number
   u.s[0] += u.s[1];
   u.u *= 1979048191;      // prime number
   u.s[1] ^= u.s[0];
   u.u *= 2079048197;      // prime number

   return u.u;
}

#endif

//______________________________________________________________________________
ULong_urt URMath::Hash(const char *txt)
{
   return Hash(txt, Int_urt(strlen(txt)));
}

//______________________________________________________________________________
Double_urt URMath::BesselI0(Double_urt x)
{
   // Compute the modified Bessel function I_0(x) for any real x.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const Double_urt p1=1.0,          p2=3.5156229,    p3=3.0899424,
                  p4=1.2067492,    p5=0.2659732,    p6=3.60768e-2,  p7=4.5813e-3;

   const Double_urt q1= 0.39894228,  q2= 1.328592e-2, q3= 2.25319e-3,
                  q4=-1.57565e-3,  q5= 9.16281e-3,  q6=-2.057706e-2,
                  q7= 2.635537e-2, q8=-1.647633e-2, q9= 3.92377e-3;

   const Double_urt k1 = 3.75;
   Double_urt ax = URMath::Abs(x);

   Double_urt y=0, result=0;

   if (ax < k1) {
      Double_urt xx = x/k1;
      y = xx*xx;
      result = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
   } else {
      y = k1/ax;
      result = (URMath::Exp(ax)/URMath::Sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselK0(Double_urt x)
{
   // Compute the modified Bessel function K_0(x) for positive real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const Double_urt p1=-0.57721566,  p2=0.42278420,   p3=0.23069756,
                  p4= 3.488590e-2, p5=2.62698e-3,   p6=1.0750e-4,    p7=7.4e-5;

   const Double_urt q1= 1.25331414,  q2=-7.832358e-2, q3= 2.189568e-2,
                  q4=-1.062446e-2, q5= 5.87872e-3,  q6=-2.51540e-3,  q7=5.3208e-4;

   if (x <= 0) {
      Error("URMath::BesselK0", "*K0* Invalid argument x = %g\n",x);
      return 0;
   }

   Double_urt y=0, result=0;

   if (x <= 2) {
      y = x*x/4;
      result = (-log(x/2.)*URMath::BesselI0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = 2/x;
      result = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselI1(Double_urt x)
{
   // Compute the modified Bessel function I_1(x) for any real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const Double_urt p1=0.5,          p2=0.87890594,   p3=0.51498869,
                  p4=0.15084934,   p5=2.658733e-2,  p6=3.01532e-3,  p7=3.2411e-4;

   const Double_urt q1= 0.39894228,  q2=-3.988024e-2, q3=-3.62018e-3,
                  q4= 1.63801e-3,  q5=-1.031555e-2, q6= 2.282967e-2,
                  q7=-2.895312e-2, q8= 1.787654e-2, q9=-4.20059e-3;

   const Double_urt k1 = 3.75;
   Double_urt ax = URMath::Abs(x);

   Double_urt y=0, result=0;

   if (ax < k1) {
      Double_urt xx = x/k1;
      y = xx*xx;
      result = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = k1/ax;
      result = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
      if (x < 0) result = -result;
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselK1(Double_urt x)
{
   // Compute the modified Bessel function K_1(x) for positive real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const Double_urt p1= 1.,          p2= 0.15443144,  p3=-0.67278579,
                  p4=-0.18156897,  p5=-1.919402e-2, p6=-1.10404e-3,  p7=-4.686e-5;

   const Double_urt q1= 1.25331414,  q2= 0.23498619,  q3=-3.655620e-2,
                  q4= 1.504268e-2, q5=-7.80353e-3,  q6= 3.25614e-3,  q7=-6.8245e-4;

   if (x <= 0) {
      Error("URMath::BesselK1", "*K1* Invalid argument x = %g\n",x);
      return 0;
   }

   Double_urt y=0,result=0;

   if (x <= 2) {
      y = x*x/4;
      result = (log(x/2.)*URMath::BesselI1(x))+(1./x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = 2/x;
      result = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselK(Int_urt n,Double_urt x)
{
   // Compute the Integer Order Modified Bessel function K_n(x)
   // for n=0,1,2,... and positive real x.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   if (x <= 0 || n < 0) {
      Error("URMath::BesselK", "*K* Invalid argument(s) (n,x) = (%d, %g)\n",n,x);
      return 0;
   }

   if (n==0) return URMath::BesselK0(x);
   if (n==1) return URMath::BesselK1(x);

   // Perform upward recurrence for all x
   Double_urt tox = 2/x;
   Double_urt bkm = URMath::BesselK0(x);
   Double_urt bk  = URMath::BesselK1(x);
   Double_urt bkp = 0;
   for (Int_urt j=1; j<n; j++) {
      bkp = bkm+Double_urt(j)*tox*bk;
      bkm = bk;
      bk  = bkp;
   }
   return bk;
}

//______________________________________________________________________________
Double_urt URMath::BesselI(Int_urt n,Double_urt x)
{
   // Compute the Integer Order Modified Bessel function I_n(x)
   // for n=0,1,2,... and any real x.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   Int_urt iacc = 40; // Increase to enhance accuracy
   const Double_urt kBigPositive = 1.e10;
   const Double_urt kBigNegative = 1.e-10;

   if (n < 0) {
      Error("URMath::BesselI", "*I* Invalid argument (n,x) = (%d, %g)\n",n,x);
      return 0;
   }

   if (n==0) return URMath::BesselI0(x);
   if (n==1) return URMath::BesselI1(x);

   if (URMath::Abs(x) > kBigPositive) return 0;

   Double_urt tox = 2/URMath::Abs(x);
   Double_urt bip = 0, bim = 0;
   Double_urt bi  = 1;
   Double_urt result = 0;
   Int_urt m = 2*((n+Int_urt(sqrt(Float_urt(iacc*n)))));
   for (Int_urt j=m; j>=1; j--) {
      bim = bip+Double_urt(j)*tox*bi;
      bip = bi;
      bi  = bim;
      // Renormalise to prevent overflows
      if (URMath::Abs(bi) > kBigPositive)  {
         result *= kBigNegative;
         bi     *= kBigNegative;
         bip    *= kBigNegative;
      }
      if (j==n) result=bip;
   }

   result *= URMath::BesselI0(x)/bi; // Normalise with BesselI0(x)
   if ((x < 0) && (n%2 == 1)) result = -result;

   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselJ0(Double_urt x)
{
   // Returns the Bessel function J0(x) for any real x.

   Double_urt ax,z;
   Double_urt xx,y,result,result1,result2;
   const Double_urt p1  = 57568490574.0, p2  = -13362590354.0, p3 = 651619640.7;
   const Double_urt p4  = -11214424.18,  p5  = 77392.33017,    p6 = -184.9052456;
   const Double_urt p7  = 57568490411.0, p8  = 1029532985.0,   p9 = 9494680.718;
   const Double_urt p10 = 59272.64853,   p11 = 267.8532712;

   const Double_urt q1  = 0.785398164;
   const Double_urt q2  = -0.1098628627e-2,  q3  = 0.2734510407e-4;
   const Double_urt q4  = -0.2073370639e-5,  q5  = 0.2093887211e-6;
   const Double_urt q6  = -0.1562499995e-1,  q7  = 0.1430488765e-3;
   const Double_urt q8  = -0.6911147651e-5,  q9  = 0.7621095161e-6;
   const Double_urt q10 =  0.934935152e-7,   q11 = 0.636619772;

   if ((ax=fabs(x)) < 8) {
      y=x*x;
      result1 = p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6))));
      result2 = p7 + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
      result  = result1/result2;
   } else {
      z  = 8/ax;
      y  = z*z;
      xx = ax-q1;
      result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
      result2 = q6 + y*(q7 + y*(q8 + y*(q9 - y*q10)));
      result  = sqrt(q11/ax)*(cos(xx)*result1-z*sin(xx)*result2);
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselJ1(Double_urt x)
{
   // Returns the Bessel function J1(x) for any real x.

   Double_urt ax,z;
   Double_urt xx,y,result,result1,result2;
   const Double_urt p1  = 72362614232.0,  p2  = -7895059235.0, p3 = 242396853.1;
   const Double_urt p4  = -2972611.439,   p5  = 15704.48260,   p6 = -30.16036606;
   const Double_urt p7  = 144725228442.0, p8  = 2300535178.0,  p9 = 18583304.74;
   const Double_urt p10 = 99447.43394,    p11 = 376.9991397;

   const Double_urt q1  = 2.356194491;
   const Double_urt q2  = 0.183105e-2,     q3  = -0.3516396496e-4;
   const Double_urt q4  = 0.2457520174e-5, q5  = -0.240337019e-6;
   const Double_urt q6  = 0.04687499995,   q7  = -0.2002690873e-3;
   const Double_urt q8  = 0.8449199096e-5, q9  = -0.88228987e-6;
   const Double_urt q10 = 0.105787412e-6,  q11 = 0.636619772;

   if ((ax=fabs(x)) < 8) {
      y=x*x;
      result1 = x*(p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6)))));
      result2 = p7    + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
      result  = result1/result2;
   } else {
      z  = 8/ax;
      y  = z*z;
      xx = ax-q1;
      result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
      result2 = q6 + y*(q7 + y*(q8 + y*(q9 + y*q10)));
      result  = sqrt(q11/ax)*(cos(xx)*result1-z*sin(xx)*result2);
      if (x < 0) result = -result;
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselY0(Double_urt x)
{
   // Returns the Bessel function Y0(x) for positive x.

   Double_urt z,xx,y,result,result1,result2;
   const Double_urt p1  = -2957821389.,  p2  = 7062834065.0, p3  = -512359803.6;
   const Double_urt p4  = 10879881.29,   p5  = -86327.92757, p6  = 228.4622733;
   const Double_urt p7  = 40076544269.,  p8  = 745249964.8,  p9  = 7189466.438;
   const Double_urt p10 = 47447.26470,   p11 = 226.1030244,  p12 = 0.636619772;

   const Double_urt q1  =  0.785398164;
   const Double_urt q2  = -0.1098628627e-2, q3  = 0.2734510407e-4;
   const Double_urt q4  = -0.2073370639e-5, q5  = 0.2093887211e-6;
   const Double_urt q6  = -0.1562499995e-1, q7  = 0.1430488765e-3;
   const Double_urt q8  = -0.6911147651e-5, q9  = 0.7621095161e-6;
   const Double_urt q10 = -0.934945152e-7,  q11 = 0.636619772;

   if (x < 8) {
      y = x*x;
      result1 = p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6))));
      result2 = p7 + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
      result  = (result1/result2) + p12*URMath::BesselJ0(x)*log(x);
   } else {
      z  = 8/x;
      y  = z*z;
      xx = x-q1;
      result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
      result2 = q6 + y*(q7 + y*(q8 + y*(q9 + y*q10)));
      result  = sqrt(q11/x)*(sin(xx)*result1+z*cos(xx)*result2);
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::BesselY1(Double_urt x)
{
   // Returns the Bessel function Y1(x) for positive x.

   Double_urt z,xx,y,result,result1,result2;
   const Double_urt p1  = -0.4900604943e13, p2  = 0.1275274390e13;
   const Double_urt p3  = -0.5153438139e11, p4  = 0.7349264551e9;
   const Double_urt p5  = -0.4237922726e7,  p6  = 0.8511937935e4;
   const Double_urt p7  =  0.2499580570e14, p8  = 0.4244419664e12;
   const Double_urt p9  =  0.3733650367e10, p10 = 0.2245904002e8;
   const Double_urt p11 =  0.1020426050e6,  p12 = 0.3549632885e3;
   const Double_urt p13 =  0.636619772;
   const Double_urt q1  =  2.356194491;
   const Double_urt q2  =  0.183105e-2,     q3  = -0.3516396496e-4;
   const Double_urt q4  =  0.2457520174e-5, q5  = -0.240337019e-6;
   const Double_urt q6  =  0.04687499995,   q7  = -0.2002690873e-3;
   const Double_urt q8  =  0.8449199096e-5, q9  = -0.88228987e-6;
   const Double_urt q10 =  0.105787412e-6,  q11 =  0.636619772;

   if (x < 8) {
      y=x*x;
      result1 = x*(p1 + y*(p2 + y*(p3 + y*(p4 + y*(p5  + y*p6)))));
      result2 = p7 + y*(p8 + y*(p9 + y*(p10 + y*(p11  + y*(p12+y)))));
      result  = (result1/result2) + p13*(URMath::BesselJ1(x)*log(x)-1/x);
   } else {
      z  = 8/x;
      y  = z*z;
      xx = x-q1;
      result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
      result2 = q6 + y*(q7 + y*(q8 + y*(q9 + y*q10)));
      result  = sqrt(q11/x)*(sin(xx)*result1+z*cos(xx)*result2);
   }
   return result;
}

//______________________________________________________________________________
Double_urt URMath::StruveH0(Double_urt x)
{
   // Struve Functions of Order 0
   //
   // Converted from CERNLIB M342 by Rene Brun.

   const Int_urt n1 = 15;
   const Int_urt n2 = 25;
   const Double_urt c1[16] = { 1.00215845609911981, -1.63969292681309147,
                             1.50236939618292819, -.72485115302121872,
                              .18955327371093136, -.03067052022988,
                              .00337561447375194, -2.6965014312602e-4,
                             1.637461692612e-5,   -7.8244408508e-7,
                             3.021593188e-8,      -9.6326645e-10,
                             2.579337e-11,        -5.8854e-13,
                             1.158e-14,           -2e-16 };
   const Double_urt c2[26] = {  .99283727576423943, -.00696891281138625,
                             1.8205103787037e-4,  -1.063258252844e-5,
                             9.8198294287e-7,     -1.2250645445e-7,
                             1.894083312e-8,      -3.44358226e-9,
                             7.1119102e-10,       -1.6288744e-10,
                             4.065681e-11,        -1.091505e-11,
                             3.12005e-12,         -9.4202e-13,
                             2.9848e-13,          -9.872e-14,
                             3.394e-14,           -1.208e-14,
                             4.44e-15,            -1.68e-15,
                             6.5e-16,             -2.6e-16,
                             1.1e-16,             -4e-17,
                             2e-17,               -1e-17 };

   const Double_urt c0  = 2/URMath::Pi();

   Int_urt i;
   Double_urt alfa, h, r, y, b0, b1, b2;
   Double_urt v = URMath::Abs(x);

   v = URMath::Abs(x);
   if (v < 8) {
      y = v/8;
      h = 2*y*y -1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n1; i >= 0; --i) {
         b0 = c1[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = y*(b0 - h*b2);
   } else {
      r = 1/v;
      h = 128*r*r -1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n2; i >= 0; --i) {
         b0 = c2[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = URMath::BesselY0(v) + r*c0*(b0 - h*b2);
   }
   if (x < 0)  h = -h;
   return h;
}

//______________________________________________________________________________
Double_urt URMath::StruveH1(Double_urt x)
{
   // Struve Functions of Order 1
   //
   // Converted from CERNLIB M342 by Rene Brun.

   const Int_urt n3 = 16;
   const Int_urt n4 = 22;
   const Double_urt c3[17] = { .5578891446481605,   -.11188325726569816,
                            -.16337958125200939,   .32256932072405902,
                            -.14581632367244242,   .03292677399374035,
                            -.00460372142093573,  4.434706163314e-4,
                            -3.142099529341e-5,   1.7123719938e-6,
                            -7.416987005e-8,      2.61837671e-9,
                            -7.685839e-11,        1.9067e-12,
                            -4.052e-14,           7.5e-16,
                            -1e-17 };
   const Double_urt c4[23] = { 1.00757647293865641,  .00750316051248257,
                            -7.043933264519e-5,   2.66205393382e-6,
                            -1.8841157753e-7,     1.949014958e-8,
                            -2.6126199e-9,        4.236269e-10,
                            -7.955156e-11,        1.679973e-11,
                            -3.9072e-12,          9.8543e-13,
                            -2.6636e-13,          7.645e-14,
                            -2.313e-14,           7.33e-15,
                            -2.42e-15,            8.3e-16,
                            -3e-16,               1.1e-16,
                            -4e-17,               2e-17,-1e-17 };

   const Double_urt c0  = 2/URMath::Pi();
   const Double_urt cc  = 2/(3*URMath::Pi());

   Int_urt i, i1;
   Double_urt alfa, h, r, y, b0, b1, b2;
   Double_urt v = URMath::Abs(x);
   
   if (v == 0) {
      h = 0;
   } else if (v <= 0.3) {
      y = v*v;
      r = 1;
      h = 1;
      i1 = (Int_urt)(-8. / URMath::Log10(v));
      for (i = 1; i <= i1; ++i) {
         h = -h*y / ((2*i+ 1)*(2*i + 3));
         r += h;
      }
      h = cc*y*r;
   } else if (v < 8) {
      h = v*v/32 -1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n3; i >= 0; --i) {
         b0 = c3[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = b0 - h*b2;
   } else {
      h = 128/(v*v) -1;
      alfa = h + h;
      b0 = 0;
      b1 = 0;
      b2 = 0;
      for (i = n4; i >= 0; --i) {
         b0 = c4[i] + alfa*b1 - b2;
         b2 = b1;
         b1 = b0;
      }
      h = URMath::BesselY1(v) + c0*(b0 - h*b2);
   }
   return h;
}


//______________________________________________________________________________
Double_urt URMath::StruveL0(Double_urt x)
{
  // Modified Struve Function of Order 0
  //  (from Kirill Filimonov)
  
  const Double_urt pi=URMath::Pi();
  
  Double_urt s=1.0;
  Double_urt r=1.0;
  
  Double_urt a0,sl0,a1,bi0;
  
  Int_urt km,i;
  
  if (x<=20.) {
    a0=2.0*x/pi;
    for (int i=1; i<=60;i++) {
      r*=(x/(2*i+1))*(x/(2*i+1));
      s+=r;
      if(URMath::Abs(r/s)<1.e-12) break;
    }
    sl0=a0*s;
  } else {
    km=int(5*(x+1.0));
    if(x>=50.0)km=25;
    for (i=1; i<=km; i++) {
      r*=(2*i-1)*(2*i-1)/x/x;
      s+=r;
      if(URMath::Abs(r/s)<1.0e-12) break;
    }
    a1=URMath::Exp(x)/URMath::Sqrt(2*pi*x);
    r=1.0;
    bi0=1.0;
    for (i=1; i<=16; i++) {
      r=0.125*r*(2.0*i-1.0)*(2.0*i-1.0)/(i*x);
      bi0+=r;
      if(URMath::Abs(r/bi0)<1.0e-12) break;
    }
    
    bi0=a1*bi0;
    sl0=-2.0/(pi*x)*s+bi0;
  }  
  return sl0;  
}

//______________________________________________________________________________
Double_urt URMath::StruveL1(Double_urt x)
{
  // Modified Struve Function of Order 1
  //  (from Kirill Filimonov)

  const Double_urt pi=URMath::Pi();
  Double_urt a1,sl1,bi1,s;
  Double_urt r=1.0;
  Int_urt km,i;
  
  if (x<=20.) {
    s=0.0;
    for (i=1; i<=60;i++) {
      r*=x*x/(4.0*i*i-1.0);
      s+=r;
      if(URMath::Abs(r)<URMath::Abs(s)*1.e-12) break;
    }
    sl1=2.0/pi*s;
  } else {
    s=1.0;
    km=int(0.5*x);
    if(x>50.0)km=25;
    for (i=1; i<=km; i++) {
      r*=(2*i+3)*(2*i+1)/x/x;
      s+=r;
      if(URMath::Abs(r/s)<1.0e-12) break;
    }
    sl1=2.0/pi*(-1.0+1.0/(x*x)+3.0*s/(x*x*x*x));
    a1=URMath::Exp(x)/URMath::Sqrt(2*pi*x);
    r=1.0;
    bi1=1.0;
    for (i=1; i<=16; i++) {
      r=-0.125*r*(4.0-(2.0*i-1.0)*(2.0*i-1.0))/(i*x);
      bi1+=r;
      if(URMath::Abs(r/bi1)<1.0e-12) break;
    }
    sl1+=a1*bi1;
  }  
  return sl1;  
}

void
URMath::Error( const string& errorMessage, double value ) {
   cout << errorMessage << ' ' << value << endl;
}

void
URMath::Error( const string& errorMessage, int iValue, double value ) {
   cout << errorMessage << ' ' << intValue << ' ' << value << endl;
}

#endif 

