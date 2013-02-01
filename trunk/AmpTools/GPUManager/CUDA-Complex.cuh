//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
// 
// Copyright Trustees of Indiana University 2010, all rights reserved
// 
// This software written by Matthew Shepherd, Ryan Mitchell, and 
//                  Hrayr Matevosyan at Indiana University, Bloomington
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
// obtained from the Trustees of Indiana University.
// 
// INDIANA UNIVERSITY AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES, 
// EXPRESS OR IMPLIED.  By way of example, but not limitation, INDIANA 
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR 
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR 
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS, 
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be 
// held liable for any liability with respect to any claim by the user or 
// any other party arising from use of the program.
//******************************************************************************

#ifndef __CUDA_COMPLEX__H__
#define __CUDA_COMPLEX__H__

#include "GPUManager/GPUCustomTypes.h"

#ifdef __CUDACC__

	#define HOSTCPU __host__
	#define DEVICE __forceinline__ __device__
	#define HOSTDEVICE __host__ __device__

#else

	#define HOSTCPU  inline
	#define DEVICE  inline
	#define HOSTDEVICE  inline

#endif

// Struct alignment is handled differently between the CUDA compiler and other
// compilers (e.g. GCC, MS Visual C++ .NET)
#ifdef __CUDACC__

	#define ALIGN(x)  __align__(x)
	
#else

	#if defined(_MSC_VER) && (_MSC_VER >= 1300)
	// Visual C++ .NET and later
		#define ALIGN(x) __declspec(align(x))
		
	#else
	
		#if defined(__GNUC__)
		// GCC
			#define ALIGN(x)  __attribute__ ((aligned (x)))
		#else
		// all other compilers
			#define ALIGN(x)
			
		#endif
		
	#endif
	
#endif


#if !defined(__DEVICE_EMULATION__) || (defined(_MSC_VER) && (_MSC_VER >= 1300))
	#define REF(x) &x
	#define ARRAYREF(x,y) (&x)[y]
#else
	#define REF(x) x
	#define ARRAYREF(x,y) x[y]
#endif

//############################################ COMPLEX NUMBERS ################################################

typedef struct ALIGN(8) WCUComplexST
{
	//Data
	GDouble m_dRe;
	GDouble m_dIm;
	
	DEVICE WCUComplexST& operator=(const int& iA) {
		m_dRe = (GDouble)iA; m_dIm = 0;
		return *this;
	};
	
	DEVICE WCUComplexST& operator=(const GDouble& fA) {
		m_dRe = (GDouble)fA; m_dIm = 0;
		return *this;
	};
	
/*	DEVICE WCUComplexST& operator=(const double& dA) {
		m_dRe = (GDouble)dA; m_dIm = 0;
		return *this;
	};
*/	
	// assignment of a pair of floats to complex
	DEVICE WCUComplexST& operator=(const GDouble ARRAYREF(a,2)) {
		m_dRe = a[0]; m_dIm = a[1];
		return *this;
	};
	
	DEVICE GDouble Re(){ return m_dRe;};
	
	DEVICE GDouble Im(){ return m_dIm;};
	
	DEVICE WCUComplexST& Set(GDouble dRe, GDouble dIm=0){ m_dRe=dRe; m_dIm=dIm; return *this;};
	
	DEVICE WCUComplexST& operator+=(const WCUComplexST& cA) {m_dRe+=cA.m_dRe; m_dIm+=cA.m_dIm;return *this;};
	
	DEVICE WCUComplexST& operator-=(const WCUComplexST& cA) {m_dRe-=cA.m_dRe; m_dIm-=cA.m_dIm;return *this;};
	
	DEVICE WCUComplexST& operator/=(const GDouble& dA) {m_dRe/=dA; m_dIm/=dA;return *this;};
	
	DEVICE WCUComplexST& operator*=(const GDouble& dA) {m_dRe*=dA; m_dIm*=dA;return *this;};
	
	DEVICE WCUComplexST& operator*=(const WCUComplexST& cA) {GDouble tRe=m_dRe,tIm=m_dIm; 
					m_dRe=tRe*cA.m_dRe-tIm*cA.m_dIm; m_dIm=tRe*cA.m_dIm+tIm*cA.m_dRe; return *this;};
	
	// assignment of a pair of floats to complex
	DEVICE WCUComplexST Conjugate() {WCUComplexST cRes ={m_dRe, -m_dIm}; return cRes;};	
		
	DEVICE GDouble ModSq(){return (m_dRe*m_dRe+m_dIm*m_dIm);};	
	
	DEVICE GDouble Mod(){return G_SQRT(m_dRe*m_dRe+m_dIm*m_dIm);};
	
	DEVICE WCUComplexST SqrtP(){ //G_SQRT according to Peng
			WCUComplexST cRes;
			if(m_dRe>=0)
			{
				cRes.m_dRe=G_SQRT(m_dRe);
				if(m_dIm<0)	cRes.m_dRe*=-1;
				cRes.m_dIm=D_EPS;//G_FABS(m_dIm);//
			}
			else
			{
				cRes.m_dRe =D_EPS;//m_dIm;//
				if(m_dIm<0)	cRes.m_dRe*=-1;
				cRes.m_dIm=G_SQRT(-m_dRe);
			}
			return cRes;
		};
		
	DEVICE	WCUComplexST Log(){ //Complex log that stays in the first Reiman sheet
			WCUComplexST cRes;
			cRes.m_dRe= G_LOG(G_SQRT(SQ(m_dRe)+SQ(m_dIm)));
			cRes.m_dIm= G_ATAN2(m_dIm,m_dRe);
			return cRes;
		};

} WCUComplex;

static DEVICE WCUComplex operator+(const WCUComplex& cA, const WCUComplex& cB){ 

  WCUComplex cRes={cA.m_dRe+cB.m_dRe, cA.m_dIm+cB.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator+(const WCUComplex& cA, const GDouble& dB) { 
  
  WCUComplex cRes={cA.m_dRe+dB,cA.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator+(const GDouble& dB, const WCUComplex& cA) { 

  WCUComplex cRes={cA.m_dRe+dB,cA.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator-(const WCUComplex& cA, const WCUComplex& cB){ 
  
  WCUComplex cRes={cA.m_dRe-cB.m_dRe, cA.m_dIm-cB.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator-(const WCUComplex& cA, const GDouble& dB) { 

  WCUComplex cRes={cA.m_dRe-dB, cA.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator-(const GDouble& dB, const WCUComplex& cA) { 

  WCUComplex cRes={dB-cA.m_dRe,-cA.m_dIm};
  return cRes;
}

static DEVICE WCUComplex operator*(const WCUComplex& cA,const  WCUComplex& cB){

	WCUComplex cRes={cA.m_dRe*cB.m_dRe-cA.m_dIm*cB.m_dIm, cA.m_dRe*cB.m_dIm+cA.m_dIm*cB.m_dRe}; 
  return cRes;
}

static DEVICE WCUComplex operator*(const WCUComplex& cA, const GDouble& dB){ 
  
  WCUComplex cRes={cA.m_dRe*dB,cA.m_dIm*dB};
  return cRes;
}

static DEVICE WCUComplex operator*(const GDouble& dB, const WCUComplex& cA){ 
  
  WCUComplex cRes={cA.m_dRe*dB,cA.m_dIm*dB};
  return cRes;
}

static DEVICE WCUComplex operator/(const WCUComplex& cA,const WCUComplex& cB){

			GDouble tDen=cB.m_dRe*cB.m_dRe+cB.m_dIm*cB.m_dIm;
			//Skip the dumb division by zero check not to loose speed.
			WCUComplex cRes={(cA.m_dRe*cB.m_dRe+cA.m_dIm*cB.m_dIm)/tDen, (-cA.m_dRe*cB.m_dIm+cA.m_dIm*cB.m_dRe)/tDen};
			return cRes;
}

static DEVICE WCUComplex operator/( const GDouble& dA, const WCUComplex& cB) {
			GDouble tDen=cB.m_dRe*cB.m_dRe+cB.m_dIm*cB.m_dIm;
			WCUComplex cRes={dA*cB.m_dRe/tDen, -dA*cB.m_dIm/tDen};
			return cRes;
}

static DEVICE WCUComplex operator/( const WCUComplex& cA, const GDouble& dB ) {

  WCUComplex cRes={cA.m_dRe/dB, cA.m_dIm/dB};
  return cRes;
}

static DEVICE WCUComplex operator~( const WCUComplex& cA ){
 
  WCUComplex cRes = {cA.m_dRe, -cA.m_dIm}; 
  return cRes;
}
	
static DEVICE WCUComplex Conjugate( const WCUComplex& cA ) { 
  
    WCUComplex cRes = {cA.m_dRe, -cA.m_dIm}; 
    return cRes;
}

static DEVICE WCUComplex TimesI(const WCUComplex& cA) { 

  WCUComplex cRes = {-cA.m_dIm,cA.m_dRe};
  return cRes;
}

static DEVICE WCUComplex TimesI(const GDouble& dA){ 
  
    WCUComplex cRes = { 0, dA };
    return cRes;
}

static DEVICE GDouble ModSq(const WCUComplex& cA){

  return (cA.m_dRe*cA.m_dRe+cA.m_dIm*cA.m_dIm);
}

static DEVICE GDouble Mod(const WCUComplex& cA){

  return G_SQRT(cA.m_dRe*cA.m_dRe+cA.m_dIm*cA.m_dIm);
}


//######################################## END OF COMPLEX NUMBERS ############################################

#endif //__CUDA_COMPLEX__H__
