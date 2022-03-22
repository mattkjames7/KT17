#include "kt14.h"

const double R0 = 1.42;	//distance to fitted subsolar magnetopause

/***********************************************************************
 * NAME : 		void KT14(x,y,z,Rsm,t1,t2,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the total KT14 magnetic field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 * 			double Rsm	Subsolar standoff distance of the magnetopause
 * 						in Rm.
 * 			double t1	Scale factor for the disk field contribution to 
 * 						the model.
 * 			double t2	Scale factor for the quasi-harris sheet field 
 * 						contribution to the model.
 *
 * OUTPUTS :
 * 			double *Bx	x-component of the KT14 field (in nT).
 * 			double *By	y-component of the KT14 field (in nT).
 * 			double *Bz	z-component of the KT14 field (in nT).
 *
 * 
 * ********************************************************************/
void KT14(	double x, double y, double z,
			double Rsm, double t1, double t2, 
			double *Bx, double *By, double *Bz) {
				
	/* rescale the input coordinates based on the size of the MP */
	double k, k3, xk, yk, zk;
	k = R0/Rsm;
	xk = x*k;
	yk = y*k;
	zk = z*k;
	k3 = k*k*k;
	
	/* get each of the field components */
	double Bxi, Byi, Bzi, Bxd, Byd, Bzd, Bxq, Byq, Bzq;
	Dipole(xk,yk,zk,&Bxi,&Byi,&Bzi);
	Disk(xk,yk,zk,t1,&Bxd,&Byd,&Bzd);
	QHSheet(xk,yk,zk,t2,&Bxq,&Byq,&Bzq);

	/* combined them */
	Bx[0] = k3*Bxi + Bxd + Bxq;
	By[0] = k3*Byi + Byd + Byq;
	Bz[0] = k3*Bzi + Bzd + Bzq;

}
