#ifndef __QHSHEET_H__
#define __QHSHEET_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shield.h"

#endif

/***********************************************************************
 * NAME : 		double QHThickness(x,y)
 * 
 * DESCRIPTION : 	Calculates the thickness of the quasi-harris sheet.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 *
 * RETURNS :
 * 			double d	Thickness of the current sheet in Rm.
 *
 * 
 * ********************************************************************/
double QHThickness(double x, double y);

/***********************************************************************
 * NAME : 		void QHSheetField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's quasi-harris sheet field at a 
 * 					given postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet field (in nT).
 * 			double *By	y-component of Mercury's QH sheet field (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet field (in nT).
 *
 * 
 * ********************************************************************/
void QHSheetField(	double x, double y, double z, 
					double *Bx, double *By, double *Bz);
					
					/***********************************************************************
 * NAME : 		void QHSheetShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the shielding field to contain the 
 * 					quasi-harris sheet field within the magnetopause.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet shield (in nT).
 * 			double *By	y-component of Mercury's QH sheet shield (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet shield (in nT).
 *
 * 
 * ********************************************************************/
void QHSheetShield(	double x, double y, double z, 
					double *Bx, double *By, double *Bz);
					
					
/***********************************************************************
 * NAME : 		void QHSheetShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the shielding field to contain the 
 * 					quasi-harris sheet field within the magnetopause.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet shield (in nT).
 * 			double *By	y-component of Mercury's QH sheet shield (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet shield (in nT).
 *
 * 
 * ********************************************************************/
void QHSheetShield(	double x, double y, double z, 
					double *Bx, double *By, double *Bz);
					
/***********************************************************************
 * NAME : 		void QHSheet(x,y,z,t2,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's quasi-harris sheet field at a 
 * 					given postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 * 			double t2	Scale factor for the field contribution to the
 * 						model.
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet field (in nT).
 * 			double *By	y-component of Mercury's QH sheet field (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet field (in nT).
 *
 * 
 * ********************************************************************/
void QHSheet(	double x, double y, double z, double t2,
				double *Bx, double *By, double *Bz);
