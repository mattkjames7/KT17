#ifndef __DISK_H__
#define __DISK_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shield.h"
#endif


/***********************************************************************
 * NAME : 		double DiskThickness(x,y)
 * 
 * DESCRIPTION : 	Calculates the disk thickness in Rm.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 *
 * RETURNS :
 * 			double d	Current disk thickness (Rm).
 *
 * 
 * ********************************************************************/
double DiskThickness(double x, double y);

/***********************************************************************
 * NAME : 		void DiskField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the disk field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of magnetic field (in nT).
 * 			double *By	y-component of magnetic field (in nT).
 * 			double *Bz	z-component of magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void DiskField(	double x, double y, double z,
				double *Bx, double *By, double *Bz);

/***********************************************************************
 * NAME : 		void DiskShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the magnetopause shielding field required
 * 					to contain the disk field.
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
void DiskShield(	double x, double y, double z,
					double *Bx, double *By, double *Bz);				

/***********************************************************************
 * NAME : 		void Disk(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's disk field at a given 
 * 					postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's disk field (in nT).
 * 			double *By	y-component of Mercury's disk field (in nT).
 * 			double *Bz	z-component of Mercury's disk field (in nT).
 *
 * 
 * ********************************************************************/
void Disk(	double x, double y, double z, double t1,
			double *Bx, double *By, double *Bz);
