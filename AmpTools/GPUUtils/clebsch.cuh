
static __device__ GDouble dfact( GDouble __x );
static __device__ GDouble clebsch( GDouble __j1, GDouble __m1, GDouble __j2,
                                   GDouble __m2, GDouble __J,  GDouble __M );

static __device__ 
GDouble dfact( GDouble __x ){

  if((__x < 0.00001) && (__x >= 0.0)) return 1.;
  if(__x < 0) return 0.;

/*
  return __x*dfact(__x - 1.);
*/

 // recursive calls not supported on GPU yet
  GDouble result = 1;
  
  while( __x >= 1.00001 ){
  
    result *= __x;
    __x -= 1;
  }

  return result;
}

/** 
 * Calcultates Clebsch-Gordon coeficents on the GPU.
 *
 * <b> Returns </b>
 *
 *  \f$\left(j_1,m_1;j_2,m_2|J,M\right) \f$ 
 *
 * Note: This function was copied from one written by Mike Williams which
 * is derived from a function written by Denis Weygand.
 */
 
static __device__ 
GDouble clebsch( GDouble __j1, GDouble __m1, GDouble __j2,
        	       GDouble __m2, GDouble __J,  GDouble __M ){

  // convert to pure integers (each 2*spin)
  int j1 = (int)(2.*__j1);
  int m1 = (int)(2.*__m1);
  int j2 = (int)(2.*__j2);
  int m2 = (int)(2.*__m2);
  int J = (int)(2.*__J);
  int M = (int)(2.*__M);

  if((__m1 + __m2) != __M) return 0.;

  GDouble n0,n1,n2,n3,n4,n5,d0,d1,d2,d3,d4,A,exp;
  int nu = 0;
  
  GDouble sum = 0;
  while(((d3=(j1-j2-M)/2+nu) < 0)||((n2=(j1-m1)/2+nu) < 0 )) { nu++;}
  while (((d1=(J-j1+j2)/2-nu) >= 0) && ((d2=(J+M)/2-nu) >= 0) 
	 &&((n1=(j2+J+m1)/2-nu) >= 0 )){
    d3=((j1-j2-M)/2+nu);
    n2=((j1-m1)/2+nu);
    d0=dfact((GDouble) nu);
    exp=nu+(j2+m2)/2;
    n0 = (GDouble) G_POW(-1.,exp);
    sum += ((n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3)));
    nu++;
  }

  if (sum == 0) return 0;

  n0 = J+1;
  n1 = dfact((GDouble) (J+j1-j2)/2);
  n2 = dfact((GDouble) (J-j1+j2)/2);
  n3 = dfact((GDouble) (j1+j2-J)/2);
  n4 = dfact((GDouble) (J+M)/2);
  n5 = dfact((J-M)/2);
  
  d0 = dfact((GDouble) (j1+j2+J)/2+1);
  d1 = dfact((GDouble) (j1-m1)/2);
  d2 = dfact((GDouble) (j1+m1)/2);
  d3 = dfact((GDouble) (j2-m2)/2);
  d4 = dfact((GDouble) (j2+m2)/2);
  
  A = ((GDouble) (n0*n1*n2*n3*n4*n5))/((GDouble) (d0*d1*d2*d3*d4));
	
  return G_SQRT( A )*sum;		
}