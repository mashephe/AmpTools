#if !defined(WIGNERD)
#define WIGNERD

#include <complex>

#include "GPUManager/GPUCustomTypes.h"

using std::complex;

GDouble wignerDSmall( GDouble aj, GDouble am, GDouble an, GDouble beta );
complex< GDouble > wignerD( int l, int m, int n, GDouble cosTheta, GDouble phi );
complex< GDouble > wignerD_1( int l, int m, int n, GDouble theta, GDouble phi );
complex< GDouble > wignerD( int l, int m, int n, GDouble cosTheta, GDouble phi1, GDouble phi2 );
complex< GDouble > Y( int l, int m, GDouble cosTheta, GDouble phi );

#endif
