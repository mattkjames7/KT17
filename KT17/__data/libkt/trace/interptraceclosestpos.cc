#include "interptraceclosestpos.h"

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
						double *xc1, double *yc1, double *zc1 ) {
							
	/* return an array of posisions along two nearby field lines which
	 * are the closest points to each element of the original field line */
	
	/* get a couple of splines (one for each field line) */
	Spline Sx0(n0,s0,x0);
	Spline Sy0(n0,s0,y0);
	Spline Sz0(n0,s0,z0);
	Spline Sx1(n1,s1,x1);
	Spline Sy1(n1,s1,y1);
	Spline Sz1(n1,s1,z1);
	
	

	/* find the closest position along the splines for each position */
	int i;
	double s0_0, s1_0;
	for (i=0;i<n;i++) {
		/* find the closest trace position to start the optimization at */
		s0_0 = ClosestS(x[i],y[i],z[i],n0,x0,y0,z0,s0);
		s1_0 = ClosestS(x[i],y[i],z[i],n1,x1,y1,z1,s1);
		
		/* Trace 0 */
		OptimizePos(x[i],y[i],z[i],bx[i],by[i],bz[i],s0_0,Sx0,Sy0,Sz0,&xc0[i],&yc0[i],&zc0[i]);
		
		/* Trace 1 */
		OptimizePos(x[i],y[i],z[i],bx[i],by[i],bz[i],s1_0,Sx1,Sy1,Sz1,&xc1[i],&yc1[i],&zc1[i]);
	}
							
}

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
				double *st) {

	int i, imin;
	double dx, dy, dz, d, dmin = INFINITY;
	for (i=0;i<nt;i++) {
		dx = x - xt[i];
		dy = y - yt[i];
		dz = z - zt[i];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		if (d < dmin) {
			imin = i;
			dmin = d;
		}
	}
	return st[imin];
					
}

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
					double bx, double by, double bz) {	/* B field unit vector */
	
	/* get the current position vector */
	double xc, yc, zc;
	Sx.Interpolate(1,&s,&xc);					
	Sy.Interpolate(1,&s,&yc);					
	Sz.Interpolate(1,&s,&zc);	
	
	/* get unit vector */
	double dx, dy, dz, d;
	dx = xc - x;
	dy = yc - y;
	dz = zc - z;
	d = sqrt(dx*dx + dy*dy + dz*dz);
	dx = dx/d;
	dy = dy/d;
	dz = dz/d;
	
	/* get the angle */
	double dot, angle;
	dot = dx*bx + dy*by + dz*bz;
	
	return fabs(M_PI/2 - acos(dot))*180.0/M_PI;
	
					
						
}

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
					double *xc, double *yc, double *zc) {
	
	/* Nelder-Mead settings */
	int MaxIter = 1000;
	double tola = 0.01;
	double tolf = 0.01;
	double alpha = 1.0;
	double gamma = 2.0;
	double rho = 0.5;
	double sigma = 0.5;
	
	/* initial/current positions */
	double s[] = {s0+0.01,s0-0.01};
	
	/* B unit vector */
	double B = sqrt(bx*bx + by*by + bz*bz);
	double bxu = bx/B;
	double byu = by/B;
	double bzu = bz/B;
	
	/* current difference between current angle and 90 degrees */
	double f[2];
	f[0] = AngleDiff(s[0],Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
	f[1] = AngleDiff(s[1],Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
		
	int best, wrst;
	bool cont = true, succ = false;
	bool shrink;
	double scnt, fcnt;
	double sr, se, sc;
	double fr, fe, fc;
	int n = 0;
	
	while (cont) {
		
		/* get best/worst indices */
		if (f[0] < f[1]) {
			best = 0;
			wrst = 1;
		} else {
			best = 1;
			wrst = 0;
		}
		shrink = false;
		
		/* centroid */
		fcnt = f[best];
		scnt = s[best];
		
		/* test reflection */
		sr = scnt + alpha*(scnt - s[wrst]);
		fr = AngleDiff(sr,Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
		
		if (fr < fcnt) {
			/* better than the best - try expanding */
			se = scnt + gamma*(sr - scnt);
			fe = AngleDiff(se,Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
			if (fe < fr) {
				/* accept expanded */
				s[wrst] = se;
				f[wrst] = fe;
			} else {
				/* accept reflected */
				s[wrst] = sr;
				f[wrst] = fr;
			}
		} else if (fr == fcnt) {
			/* accept the reflection */
			s[wrst] = sr;
			f[wrst] = fr;
		} else {
			/* worse than the worst or second worst - contract */
			if ((fr > fcnt) && (fr < f[wrst])) {
				/* outside contraction */
				sc = scnt + rho*(sr - scnt);
				fc = AngleDiff(sc,Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
				if (fc <= fr) {
					/* accept contraction */
					s[wrst] = sc;
					f[wrst] = fc;
				} else {
					shrink = true;
				}
			} else {
				/* inside contraction */
				sc = scnt + rho*(s[wrst] - scnt);
				fc = AngleDiff(sc,Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
				if (fc < f[wrst]) {
					/* accept contraction */
					s[wrst] = sc;
					f[wrst] = fc;
				} else {
					shrink = true;
				}				
			}
			
			/* shrink */
			if (shrink) {
				s[wrst] = s[wrst] + sigma*(scnt - s[wrst]);
				f[wrst] = AngleDiff(s[wrst],Sx,Sy,Sz,x,y,z,bxu,byu,bzu);
			}
		}
		
		/* check if we have more or less converged */
		if ((fabs(0.5*(f[1]+f[0])) <= tola) && (fabs(f[1]-f[0]) <= tolf)) {
			cont = false;
			succ = true;
		}
		
		/* or if we have ran out of iteration */
		if (n >= MaxIter) {
			cont = false;
		}
		n++;
	}

	/* use the average position */
	scnt = 0.5*(s[0] + s[1]);
	Sx.Interpolate(1,&scnt,xc);					
	Sy.Interpolate(1,&scnt,yc);					
	Sz.Interpolate(1,&scnt,zc);
	//double a = bxu*(x-xc[0]) + byu*(y-yc[0]) + bzu*(z-zc[0]);
	//double b = acos(a)*180.0/M_PI;
		
	//printf("s0: %f, scnt: %f, niter %d, d90: %f, da: %f, dot: %f, angle: %f\n",s0,scnt,n,fabs(0.5*(f[1]+f[0])),fabs(f[1]-f[0]),a,b);
}
