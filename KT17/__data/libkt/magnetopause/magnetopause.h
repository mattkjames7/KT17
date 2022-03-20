#ifndef __MAGNETOPAUSE_H__
#define __MAGNETOPAUSE_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#endif


/***********************************************************************
 * NAME : double MagnetopauseR(theta,Rsm,alpha)
 * 
 * DESCRIPTION : Calculate the radius of the magnetopause at some angle.
 * 
 * INPUTS : 
 * 		double theta	Angle between y/z and x coordinates (rad).
 * 		double Rsm		Subsolar magnetopause distance (Rm).
 * 		double alpha	Flaring parameter - default is 0.5.
 *
 * RETURNS :
 *		double Rmp		Magnetopause radial coordinate at theta.
 * 
 * ********************************************************************/
double MagnetopauseR( double theta, double Rsm, double alpha);

/***********************************************************************
 * NAME : double MagnetopauseR(theta,Rsm)
 * 
 * DESCRIPTION : Calculate the radius of the magnetopause at some angle. 
 * 
 * INPUTS : 
 * 		double theta	Angle between y/z and x coordinates (rad).
 * 		double Rsm		Subsolar magnetopause distance (Rm).
 *
 * RETURNS :
 *		double Rmp		Magnetopause radial coordinate at theta.
 * 
 * ********************************************************************/
double MagnetopauseR( double theta, double Rsm);
