#ifndef __MODEL_H__
#define __MODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include "ktmodel.h"

#endif

/* this object will be configured for each trace */
extern KTModel ktmodel;

/***********************************************************************
 * NAME : void ModelField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Wrapper for the KT14/17 model field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of the KT14/17 field (in nT).
 * 			double *By	y-component of the KT14/17 field (in nT).
 * 			double *Bz	z-component of the KT14/17 field (in nT).
 *
 * 
 * ********************************************************************/
void ModelField(double x, double y, double z,	
				double *Bx, double *By, double *Bz);
