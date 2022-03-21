#ifndef __WITHINMP_H__
#define __WITHINMP_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "magnetopause.h"

#endif

extern "C" {
	/***********************************************************************
	 * NAME : bool WithinMPRT(r,theta,Rsm)
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
	bool WithinMPRT(double r, double theta, double Rsm);

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
	bool WithinMP(double x, double y, double z, double Rsm);
}
