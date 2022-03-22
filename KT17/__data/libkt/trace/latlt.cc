#include "latlt.h"

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
void LatLT(double x, double y, double z, double *Lat, double *LT) {
	
	double rho;
	if (isnan(x)) {
		Lat[0] = NAN;
		LT[0] = NAN;
	} else {
		rho = sqrt(x*x + y*y);
		Lat[0] = atan2(z,rho)*180.0/M_PI;
		LT[0] = fmod(atan2(-y,-x)*12.0/M_PI + 24.0,24.0);
	}
}

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
void LshellMLT(double x, double y, double z, double *L, double *MLT) {
	
	if (isnan(x)) {
		L[0] = NAN;
		MLT[0] = NAN;
	} else {
		L[0] = sqrt(x*x + y*y + z*z);
		MLT[0] = fmod(atan2(-y,-x)*12.0/M_PI+24.0,24.0);
	}
}
