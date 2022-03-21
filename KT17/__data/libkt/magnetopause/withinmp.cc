#include "withinmp.h"

/***********************************************************************
 * NAME : bool WithinMP(r,theta,Rsm)
 * 
 * DESCRIPTION : Checks whether a position is within the magnetopause.
 * 
 * INPUTS : 
 * 		double r		Radial coordinate.
 * 		double theta	Angle between y/z and x coordinates (rad).
 * 		double Rsm		Subsolar magnetopause distance (Rm).
 *
 * RETURNS :
 *		bool wmp		true if within magnetopause.
 * 
 * ********************************************************************/
bool WithinMPRT(double r, double theta, double Rsm) {
	
	double rmp;
	
	rmp = MagnetopauseR(theta,Rsm);
	
	return (r <= rmp);
	
}

/***********************************************************************
 * NAME : bool WithinMP(x,y,z,Rsm)
 * 
 * DESCRIPTION : Checks whether a position is within the magnetopause.
 * 
 * INPUTS : 
 * 		double x		x-coordinate (Rm).
 * 		double y		y-coordinate (Rm).
 * 		double z		z-coordinate (Rm).
 * 		double Rsm		Subsolar magnetopause distance (Rm).
 *
 * RETURNS :
 *		bool wmp		true if within magnetopause.
 * 
 * ********************************************************************/
bool WithinMP(double x, double y, double z, double Rsm) {
	
	double theta, r, rho;
	rho = y*y + z*z;
	r = x*x + rho;
	rho = sqrt(rho);
	r = sqrt(r);
	theta = atan2(rho,x);
	
	return WithinMPRT(r,theta,Rsm);

}
