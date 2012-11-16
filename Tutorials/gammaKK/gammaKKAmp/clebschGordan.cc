
#include "./clebschGordan.h"

#include <math.h>

/* Name: s3j
**       Evaluates 3j symbol
**
** Author: Riccardo Gusmeroli (web.address@libero.it)
**
** Notes: 
**     - defining S3J_TEST enables the compilation of a very small test suite.
**     - the maximum allowed factorial is S3J_MAX_FACT (currently 25!). 
**
**
** This program is free software; you can redistribute it and/or  
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2  
** of the License, or (at your option) any later version.
** This program is distributed in the hope that it will be useful, 
** but WITHOUT ANY WARRANTY; without even the implied warranty of 
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
** GNU General Public License for more details.
** You should have received a copy of the GNU General Public License 
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. */


#define S3J_0		1e-10

#define S3J_MAX_FACT	25
#define S3J_TEST

#define S3J_EQUAL(a,b)		(fabs((a)-(b))<S3J_0)
#define S3J_MAX(a,b,c,ris)	(((a)>(b)?(ris=(a)):(ris=(b)))>(c)?ris:(ris=(c)))
#define S3J_MIN(a,b,c,ris)	(((a)<(b)?(ris=(a)):(ris=(b)))<(c)?ris:(ris=(c)))


double s3j(double j1, double j2, double j3, 
		   double m1, double m2, double m3) {
	
	/*  ( j1 j2 j3 )
	(            ) = delta(m1+m2+m3,0) * (-1)^(j1-j2-m3) * 
    ( m1 m2 m3 )
	
	+-
	|  (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)! 
    * | -------------------------------------- ...
	|
	+-
	-+ 1/2
	(j1-m1)! (j1+m1)! (j2-m2)! (j2+m2)! (j3-m3)! (j3+m3)!  |
	... ------------------------------------------------------- |     * 
	(j1+j2+j3+1)!                       |
	-+
	
	+---
	\                       (-1)^k
	*    |   ---------------------------------------------------------------------
	/      k! (j1+j2-j3-k)! (j1-m1-k)! (j2+m2-k)! (j3-j2+m1+k)! (j3-j1-m2+k)!
	+---
	k
	
	Where factorials must have non-negative integral values:
	
	j1+j2-j3   >= 0		j1-j2+j3   >= 0		-j1+j2+j3 >= 0		j1+j2+j3+1 >= 0
	k          >= 0		j1+j2-j3-k >= 0		j1-m1-k   >= 0		j2+m2-k    >= 0
	j3-j2+m1+k >= 0		j3-j1-m2+k >= 0
	
	The 3j symbol is therefore non-null if
	
				j1+j2 >= j3		(1)
				j1+j3 >= j2		(2)
				j2+j3 >= j1		(3)
	
	and k values in the sum must be such that
	
				k <= j1+j2-j3		(4)			k >= 0				(7)
				k <= j1-m1			(5)			k >= -j3+j2-m1		(8)
				k <= j2+m2			(6)			k >= -j3+j1+m2		(9)
	
	If no values of k satisfy the (4) to (9), the result is null because the sum is null,
	otherwise one can find   kmin < kmax   such that
	
	kmin <= k <= kmax
	
	(4) to (6) => kmin=MAX(j1+j2-j3,  j1-m1,       j2+m2     )
	(7) to (9) => kmax=MIN(0,         -j3+j2-m1,   -j3+j1+m2 )
	
	The condition kmin < kmax includes (1) to (3) because
	
	(4) and (7)    =>    (1)
	(5) and (8)    =>    (2)
	(6) and (9)    =>    (3)
	
	Once the values of kmin and kmax are found, the only "selection rule" is kmin<kmax.
	*/
	
	int k, kmin, kmax;
	int jpm1, jmm1, jpm2, jmm2, jpm3, jmm3;
	int j1pj2mj3, j3mj2pm1, j3mj1mm2;
	double ris, mult, f[S3J_MAX_FACT];
	
	f[0]=1.0;
	mult=1.0;
	for (k=1; k<S3J_MAX_FACT; ++k) {
		
		f[k]=f[k-1]*mult;
		mult+=1.0;
	}
	
	jpm1=(int)(j1+m1);
	if (!S3J_EQUAL(jpm1,j1+m1)) return 0.0;
	
	jpm2=(int)(j2+m2);
	if (!S3J_EQUAL(jpm2,j2+m2)) return 0.0;
	
	jpm3=(int)(j3+m3);
	if (!S3J_EQUAL(jpm3,j3+m3)) return 0.0;
	
	jmm1=(int)(j1-m1);
	if (!S3J_EQUAL(jmm1,j1-m1)) return 0.0;
	
	jmm2=(int)(j2-m2);
	if (!S3J_EQUAL(jmm2,j2-m2)) return 0.0;
	
	jmm3=(int)(j3-m3);
	if (!S3J_EQUAL(jmm3,j3-m3)) return 0.0;
	
	/* delta(m1+m2+m3,0) */
	if ((jpm1-jmm1+jpm2-jmm2+jpm3-jmm3)!=0) return 0.0;
	
	/* j1+j2-j3 = (j1+j2-j3) + (m1+m2+m3) = jpm1+jpm2-jmm3 */
	j1pj2mj3=jpm1+jpm2-jmm3;
	
	/* j3-j2+m1 = (j3-j2+m1) - (m1+m2+m3) = jmm3-jpm2 */
	j3mj2pm1=jmm3-jpm2;
	
	/* j3-j1-m2 = (j3-j1-m2) + (m1+m2+m3) = jpm3-jmm1 */
	j3mj1mm2=jpm3-jmm1;
	
	
	S3J_MAX(-j3mj2pm1, -j3mj1mm2,   0,         kmin);
	S3J_MIN(j1pj2mj3,  jmm1,       jpm2,       kmax);
	if (kmin>kmax) return 0.0;	
	
	ris=0.0;
	if (kmin%2==0) mult=1.0;
	else mult=-1.0;
	for (k=kmin; k<=kmax; ++k) {
		
		ris+=mult/(f[k]*f[j1pj2mj3-k]*f[jmm1-k]*f[jpm2-k]*f[j3mj2pm1+k]*f[j3mj1mm2+k]);
		
		mult=-mult;
	}
	
	/* (-1)^(j1-j2-m3)=(-1)^(j1-j2-m3+m1+m2+m3)=(-1)^(jpm1-jmm2) */
	if ((jpm1-jmm2)%2!=0) ris=-ris;
	
	ris*=sqrt(f[j1pj2mj3]*f[jpm1-jmm2+jpm3]*f[-jmm1+jpm2+jpm3]*
			  f[jpm1]*f[jpm2]*f[jpm3]*f[jmm1]*f[jmm2]*f[jmm3]/
			  f[jpm1+jpm2+jpm3+1]);
	
	return ris;
}


double clebschGordan(int ij1, int ij2, int im1, int im2, int ij, int im) {
	
	int esp;
	double cgris;
	double j1,j2,m1,m2,j,m;
	
	j1 = (double) ij1;
	j2 = (double) ij2;
	m1 = (double) im1;
	m2 = (double) im2;
	j  = (double) ij;
	m  = (double) im;
	
	esp=(int)(j1-j2+m);
	if (!S3J_EQUAL(esp,j1-j2+m)) return 0;
	
	if (esp%2==0) cgris=1.0;
	else cgris=-1.0;
	
	cgris*=sqrt(2*j+1)*s3j(j1,j2,j,m1,m2,-m);
		
	return cgris;
}


