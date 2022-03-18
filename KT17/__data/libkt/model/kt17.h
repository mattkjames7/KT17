#ifndef __KT17_H__
#define __KT17_H__
#include <stdio.h>
#include <stdlib.h>
#include "kt14.h"

#endif


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
			double *Bx, double *By, double *Bz);
