#ifndef __KT17_H__
#define __KT17_H__
#include <stdio.h>
#include <stdlib.h>
#include "kt14.h"

#endif


/***********************************************************************
 * NAME : 		void KT14Params(Rsun,DistIndex,Rsm,t1,t2)
 * 
 * DESCRIPTION : 	Converts KT17 parameters to KT14 ones.
 * 
 * INPUTS : 
 * 		double Rsun		Radial distance of Mercury from the Sun (AU).
 * 		double DistIndex	Anderson et al 2013 disturbance index in the
 * 						range 0-97 doi: 10.1002/ggge.20242.	
 *
 * OUTPUTS :
 * 		double Rsm	Subsolar standoff distance of the magnetopause
 * 					in Rm.
 * 		double t1	Scale factor for the disk field contribution to 
 * 					the model.
 * 		double t2	Scale factor for the quasi-harris sheet field 
 * 					contribution to the model.
 *
 * 
 * ********************************************************************/
void KT14Params(	double Rsun, double DistIndex,
					double *Rsm, double *t1, double *t2);
					

/***********************************************************************
 * NAME : 		void KT17Params(Rsm,t1,t2,Rsun,DistIndex)
 * 
 * DESCRIPTION : 	Converts KT14 parameters to KT17 ones.
 * 
 * INPUTS : 
 * 		double Rsm	Subsolar standoff distance of the magnetopause
 * 					in Rm.
 * 		double t1	Scale factor for the disk field contribution to 
 * 					the model.
 * 		double t2	Scale factor for the quasi-harris sheet field 
 * 					contribution to the model.
 *
 * OUTPUTS :
 * 		double Rsun		Radial distance of Mercury from the Sun (AU).
 * 		double DistIndex	Anderson et al 2013 disturbance index in the
 * 						range 0-97 doi: 10.1002/ggge.20242.	
 *
 * 
 * ********************************************************************/
void KT17Params(	double Rsm, double t1, double t2,
					double *Rsun, double *DistIndex);					
					
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
