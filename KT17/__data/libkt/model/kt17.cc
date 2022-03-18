#include "kt17.h"

/***********************************************************************
 * NAME : 		void KT17(x,y,z,Rsun,DistIndex,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the total KT17 magnetic field.
 * 
 * INPUTS : 
 * 		double x		x-MSM coordinate (in Rm).
 * 		double y		y-MSM coordinate (in Rm).
 * 		double z		z-MSM coordinate (in Rm).
 * 		double Rsun		Radial distance of Mercury from the Sun (AU).
 * 		double DistIndex	Anderson et al 2013 disturbance index in the
 * 						range 0-97 doi: 10.1002/ggge.20242.	
 *
 * OUTPUTS :
 * 			double *Bx	x-component of the KT17 field (in nT).
 * 			double *By	y-component of the KT17 field (in nT).
 * 			double *Bz	z-component of the KT17 field (in nT).
 *
 * 
 * ********************************************************************/
void KT17(	double x, double y, double z,
			double Rsun, double DistIndex, 
			double *Bx, double *By, double *Bz) {
			
	/* calculate the parameters for the KT14 model */
	double f, Rsm, t1, t2;
	f = 2.06873 - 0.00279*DistIndex;
	Rsm = f*pow(Rsun,1.0/3.0);
	t1 = 6.495 + 0.0229*DistIndex;
	t2 = 1.6245 + 0.0088*DistIndex;
	
	/* call the KT14 model */
	KT14(x,y,z,Rsm,t1,t2,Bx,By,Bz);	
	
}
