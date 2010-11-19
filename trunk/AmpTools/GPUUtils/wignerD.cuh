
#ifndef WIGNERD
#define WIGNERD

#include "GPUManager/GPUCustomTypes.h"
#include "GPUManager/CUDA-Complex.cuh"

static __device__ WCUComplex wignerD( int l, int m, int n, GDouble cosTheta, GDouble phi );
static __device__ GDouble    wignerDSmall( GDouble aj, GDouble am, GDouble an, GDouble beta );
static __device__ WCUComplex Y( int l, int m, GDouble cosTheta, GDouble phi );

static __device__ GDouble
wignerDSmall( GDouble aj, GDouble am, GDouble an, GDouble beta ){

	// Calculates the beta-term
	//                         d j mn (beta)
	// in the matrix element of the finite rotation operator
	// (Wigner's D-function), according to formula 4.3.1(3) in
	// D.A. Varshalovich, A.N. Moskalev, and V.K. Khersonskii,
	// Quantum Theory of Angular Momentum, World Scientific,
	// Singapore 1988.
	// CERNLIB DDJMNB function translated from Fortran to C++ by Rene Brun
	  
	GDouble f = 8.72664625997164788e-3;    
  
  GDouble fcl[51] = { 0 , 0 ,
		6.93147180559945309e-1 ,1.79175946922805500e00,
		3.17805383034794562e00 ,4.78749174278204599e00,
		6.57925121201010100e00 ,8.52516136106541430e00,
		1.06046029027452502e01 ,1.28018274800814696e01,
		1.51044125730755153e01 ,1.75023078458738858e01,
		1.99872144956618861e01 ,2.25521638531234229e01,
		2.51912211827386815e01 ,2.78992713838408916e01,
		3.06718601060806728e01 ,3.35050734501368889e01,
		3.63954452080330536e01 ,3.93398841871994940e01,
		4.23356164607534850e01 ,4.53801388984769080e01,
		4.84711813518352239e01 ,5.16066755677643736e01,
		5.47847293981123192e01 ,5.80036052229805199e01,
		6.12617017610020020e01 ,6.45575386270063311e01,
		6.78897431371815350e01 ,7.12570389671680090e01,
		7.46582363488301644e01 ,7.80922235533153106e01,
		8.15579594561150372e01 ,8.50544670175815174e01,
		8.85808275421976788e01 ,9.21361756036870925e01,
		9.57196945421432025e01 ,9.93306124547874269e01,
		1.02968198614513813e02 ,1.06631760260643459e02,
		1.10320639714757395e02 ,1.14034211781461703e02,
		1.17771881399745072e02 ,1.21533081515438634e02,
		1.25317271149356895e02 ,1.29123933639127215e02,
		1.32952575035616310e02 ,1.36802722637326368e02,
		1.40673923648234259e02 ,1.44565743946344886e02,
		1.48477766951773032e02};
	
	int jpm = int(aj+am);
	int jpn = int(aj+an);
	int jmm = int(aj-am);	
	
	int jmn = int(aj-an);
	int mpn = int(am+an);
	
	GDouble r = 0;
	if (beta == 0) 
	{
		if (jpm == jpn) r = 1;
  } 
	else if (beta == 180) 
	{
		if (jpm == jmn) 
		{
			r = 1;
			if ( (jpm > 0 ? jpm : -jpm ) % 2 == 1 ) r = -1;
		}
  } 
	else if (beta == 360)
	{
		if (jpm == jpn)
		{
			r = 1;
      if ( (jpm > 0 ? jpm : -jpm ) % 2 == 1 ) r = -1;
		}
  } 
	else
	{
		GDouble b  = f*beta;
		GDouble s  = G_LOG(G_SIN(b));
		GDouble c  = G_LOG(G_FABS(G_COS(b)));
		GDouble rt = 0.5*(fcl[jpm]+fcl[jmm]+fcl[jpn]+fcl[jmn]);
		int k0    = ( 0 > mpn ? 0 : mpn ); //max( 0 , mpn )
		int kq    = k0+jpm;
		if (beta > 180) kq += mpn;
		GDouble q  = 1;
		if (kq%2 == 1) q = -1;
		kq = k0+k0;
		GDouble cx = kq-mpn;
		GDouble sx = jpm+jpn-kq;
		for( int k = k0 ; k <= ( jpm < jpn ? jpm : jpn ); k++ ) 
		{
			r  += q*G_EXP(rt-fcl[k]-fcl[jpm-k]-fcl[jpn-k]-fcl[k-mpn]+ cx*c+sx*s);
			cx += 2;
			sx -= 2;
			q   = -q;
		}
	}
  
	return r;
}

static __device__ WCUComplex
wignerD( int l, int m, int n, GDouble cosTheta, GDouble phi ){
	
  GDouble dpart = wignerDSmall( l, m, n, G_ACOS( cosTheta ) * 180.0 / PI );
	
  WCUComplex res = { G_COS( -1.0 * m * phi ) * dpart, 
                     G_SIN( -1.0 * m * phi ) * dpart  };

  return res;
}

static __device__ WCUComplex
Y( int l, int m, GDouble cosTheta, GDouble phi ){

  return G_SQRT( (2*l+1 ) / (4*PI) ) * 
    Conjugate( wignerD( l, m, 0, cosTheta, phi ) );
}

#endif
