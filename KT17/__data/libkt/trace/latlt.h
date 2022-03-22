#ifndef __LATLT_H__
#define __LATLT_H__
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#endif


/***********************************************************************
 * NAME : void LatLT(x,y,z,Lat,LT)
 * 
 * DESCRIPTION : Calculates the latitude and local time of a position.
 * 
 * INPUTS : 
 *		double x	x-position
 *		double y	y-position
 *		double z	z-position
 * 
 * OUPUTS :
 *		double Lat	Latitude in degrees.
 * 		double LT 	Local time in hours.
 * 	  
 * ********************************************************************/
void LatLT(double x, double y, double z, double *Lat, double *LT);

/***********************************************************************
 * NAME : void LshellMLT(x,y,z,L,MLT)
 * 
 * DESCRIPTION : Calculates the L-shell and local time of a position
 * 				on the magnetic equatorial plane.
 * 
 * INPUTS : 
 *		double x	x-position
 *		double y	y-position
 *		double z	z-position
 * 
 * OUPUTS :
 *		double L	L-shell.
 * 		double MLT 	Local time in hours.
 *		  
 * ********************************************************************/
void LshellMLT(double x, double y, double z, double *L, double *MLT);
