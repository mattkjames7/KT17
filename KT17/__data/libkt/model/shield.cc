#include "shield.h"

/***********************************************************************
 * NAME : 		void ShieldField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the magnetopause shielding field required
 * 					to contain a magnetospheric source.
 * 
 * INPUTS : 
 * 			double x		x-MSM coordinate (in Rm).
 * 			double y		y-MSM coordinate (in Rm).
 * 			double z		z-MSM coordinate (in Rm).
 * 			int n			dimension length(s) of A and P.
 *			double A[][n]	Shielding coefficients.
 * 			double P[]		Shielding coefficients.
 * 
 * OUTPUTS :
 * 			double *Bx	x-component of magnetopause field (in nT).
 * 			double *By	y-component of magnetopause field (in nT).
 * 			double *Bz	z-component of magnetopause field (in nT).
 *
 * 
 * ********************************************************************/
void ShieldField(	double x, double y, double z,
					int n, const double A[], const double P[],
					double *Bx, double *By, double *Bz) {
	
	/* set the outputs equal to zero */
	Bx[0] = 0.0;
	By[0] = 0.0;
	Bz[0] = 0.0;

	/* some temporaray variables */
	int i, k, iA;
	double pik, Aexpx, cosPy, sinPz, cosPz, sinPy, Pi2;
	
	/* loop through sum (equation 12 of Korth et al 2015) */
	for (i=0;i<n;i++) {
		Pi2 = P[i]*P[i];
		for (k=0;k<n;k++) {
			iA = i*n + k;
			pik = sqrt(Pi2 + P[k]*P[k]);
			Aexpx = A[iA]*exp(pik*x);
			cosPy = cos(P[i]*y);
			cosPz = cos(P[k]*z);
			sinPy = sin(P[i]*y);
			sinPz = sin(P[k]*z);
			
			Bx[0] += Aexpx*pik*cosPy*sinPz;
			By[0] +=-Aexpx*P[i]*sinPy*sinPz;
			Bz[0] += Aexpx*P[k]*cosPy*cosPz;
		}
	}
}
