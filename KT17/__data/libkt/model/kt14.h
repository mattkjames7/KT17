#ifndef __KT14_H__
#define __KT14_H__
#include <stdio.h>
#include <stdlib.h>
#include "dipole.h"
#include "disk.h"
#include "qhsheet.h"

#endif


/***********************************************************************
 * NAME : 		void KT14(x,y,z,Rsm,t1,t2,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the total KT14 magnetic field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 * 			double Rsm	Subsolar standoff distance of the magnetopause
 * 						in Rm.
 * 			double t1	Scale factor for the disk field contribution to 
 * 						the model.
 * 			double t2	Scale factor for the quasi-harris sheet field 
 * 						contribution to the model.
 *
 * OUTPUTS :
 * 			double *Bx	x-component of the KT14 field (in nT).
 * 			double *By	y-component of the KT14 field (in nT).
 * 			double *Bz	z-component of the KT14 field (in nT).
 *
 * 
 * ********************************************************************/
void KT14(	double x, double y, double z,
			double Rsm, double t1, double t2, 
			double *Bx, double *By, double *Bz);
