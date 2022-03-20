#ifndef __DIPOLE_H__
#define __DIPOLE_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shield.h"
#endif

/***********************************************************************
 * NAME : 		void DipoleField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's dipole field at a given 
 * 					postion.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's magnetic field (in nT).
 * 			double *By	y-component of Mercury's magnetic field (in nT).
 * 			double *Bz	z-component of Mercury's magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void DipoleField(	double x, double y, double z,
					double *Bx, double *By, double *Bz);
					
/***********************************************************************
 * NAME : 		void DipoleShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the magnetopause shielding field required
 * 					to contain the dipole field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of magnetopause field (in nT).
 * 			double *By	y-component of magnetopause field (in nT).
 * 			double *Bz	z-component of magnetopause field (in nT).
 *
 * 
 * ********************************************************************/
void DipoleShield(	double x, double y, double z,
					double *Bx, double *By, double *Bz);
					
					
/***********************************************************************
 * NAME : 		void Dipole(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's dipole field at a given 
 * 					postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's magnetic field (in nT).
 * 			double *By	y-component of Mercury's magnetic field (in nT).
 * 			double *Bz	z-component of Mercury's magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void Dipole(	double x, double y, double z, 
				double *Bx, double *By, double *Bz);
