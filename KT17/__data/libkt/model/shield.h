#ifndef __SHIELD_H__
#define __SHIELD_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#endif


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
					double *Bx, double *By, double *Bz);
