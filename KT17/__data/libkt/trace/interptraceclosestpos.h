#ifndef __INTERPTRACECLOSESTPOS_H__
#define __INTERPTRACECLOSESTPOS_H__
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include "../spline/spline.h"
#include <math.h>
#endif

/***********************************************************************
 * NAME : void interptraceClosestPos(n,x,y,z,bx,by,bz,n0,x0,y0,z0,s0,
 * 								n1,x1,y1,z1,s1,xc0,yc0,zc0,xc1,yc1,zc1)
 * 
 * DESCRIPTION : Uses interpolation splines to work out the closest 
 * 			positions along adjacent field lines to help with working 
 * 			out h_alpha. The output positions have the same length as
 * 			the original field line (n,x,y,z) but are positions along
 * 			adjacent field lines 0 and 1.
 * 
 * INPUTS : 
 *		int n		Length of original field line.
 * 		double *x	x-coordinate of field line
 * 		double *y	y-coordinate of field line
 * 		double *z	z-coordinate of field line
 * 		double *bx	x-component of the magnetic field.
 * 		double *by	y-component of the magnetic field.
 * 		double *bz	z-component of the magnetic field.
 * 		int n0		length of adjacent field line 0.
 * 		double *x0	x-coordinate of field line 0
 * 		double *y0	y-coordinate of field line 0
 * 		double *z0	z-coordinate of field line 0
 * 		double *s0	distance along field line 0
 * 		int n1		length of adjacent field line 1.
 * 		double *x1	x-coordinate of field line 1
 * 		double *y1	y-coordinate of field line 1
 * 		double *z1	z-coordinate of field line 1
 * 		double *s1	distance along field line 1
 * 
 * OUPUTS :
 * 		double *xc0	x-coordinate along field line 0
 * 		double *yc0	y-coordinate along field line 0
 * 		double *zc0	z-coordinate along field line 0
 * 		double *xc1	x-coordinate along field line 1
 * 		double *yc1	y-coordinate along field line 1
 * 		double *zc1	z-coordinate along field line 1
 *		  
 * ********************************************************************/
void interptraceClosestPos(	int n, double *x, double *y, double *z,
						double *bx, double *by, double *bz,
						int n0, double *x0, double *y0, double *z0, double *s0,
						int n1, double *x1, double *y1, double *z1, double *s1,
						double *xc0, double *yc0, double *zc0,
						double *xc1, double *yc1, double *zc1 );
	
/***********************************************************************
 * NAME : double ClosestS(x,y,z,nt,xt,yt,zt,st)
 * 
 * DESCRIPTION : Find the distanc (s) along an adjacent field line 
 * 		closest to a point along the original field line
 * 
 * INPUTS : 
 * 		double x	x position
 * 		double y	y position
 * 		double z	z position
 * 		int nt		number of elements in adjacent field line.
 * 		double *xt	x-position along adjacent field line
 * 		double *yt	y-position along adjacent field line
 * 		double *zt	z-position along adjacent field line
 * 		double *st	distance along adjacent field line
 *	
 * RETURNS :
 * 		double stmin 	distance along adjacent field line closest to 
 * 						position defined by x,y,z. 
 *		  
 * ********************************************************************/					
double ClosestS(double x, double y, double z,
				int nt, double *xt, double *yt, double *zt,
				double *st);

/***********************************************************************
 * NAME : double AngleDiff(s,Sx,Sy,Sz,x,y,z,bx,by,bz)
 * 
 * DESCRIPTION : Calculate the difference in angle between the unit 
 * 		normal vector of the position on the original field (unit field 
 * 		vector + 90 degrees) and the unit vector from that position to
 * 		a position along the adjacent field line. Ideally this would be 
 * 		zero when we have the correct position.
 * 
 * INPUTS : 
 * 		double s		position along adjacent field line to test.
 * 		Spline Sx		x-position spline along adjacent field line.
 * 		Spline Sy		y-position spline along adjacent field line.
 * 		Spline Sz		z-position spline along adjacent field line.
 * 		double x		original x position
 * 		double y		original y position
 * 		double z		original z position
 * 		double bx		x-component of B
 * 		double by		y-component of B
 * 		double bz		z-component of B
 *	
 * RETURNS :
 *		double diff		difference in angle between the unit normal of 
 * 						the original field line and the unit vector from
 * 						the point on the original FL to the adjacent FL 
 * ********************************************************************/			
double AngleDiff( 	double s,								/* current position along the field line */
					Spline Sx, Spline Sy, Spline Sz,	/* Splines converting s to a  vector */
					double x, double y, double z,		/* this is the position along the original field line */
					double bx, double by, double bz);
					
/***********************************************************************
 * NAME : void OptimizePos(x,y,z,bx,by,bz,s0,Sx,Sy,Sz,xc,yc,zc)
 * 
 * DESCRIPTION : Uses Nelder-Mead algorithm to find the point along an
 * 		adjacent field line that is close to the original field line
 * 		normal.
 * 
 * INPUTS : 
 *		double x	x-position along original field line.
 *		double y	y-position along original field line.
 *		double z	z-position along original field line.
 * 		double bx		x-component of B
 * 		double by		y-component of B
 * 		double bz		z-component of B
 * 		double s0		starting position along adjacent field line.
 * 		Spline Sx		x-position spline along adjacent field line.
 * 		Spline Sy		y-position spline along adjacent field line.
 * 		Spline Sz		z-position spline along adjacent field line.
 * 
 * 
 * OUPUTS :
 * 		double *xc		closest x position along adjacent field line.
 * 		double *yc		closest y position along adjacent field line.
 * 		double *zc		closest z position along adjacent field line.
 *		  
 * ********************************************************************/					
void OptimizePos(	double x, double y, double z,
					double bx, double by, double bz,
					double s0, 
					Spline Sx, Spline Sy, Spline Sz,
					double *xc, double *yc, double *zc);
