#include "disk.h"

const double xshift = 0.3; // shift in x-direction
const double sc = 7.0;		//scaling factor
const double d0 = 0.09;	//half-width of current sheet in Z at inner edge of tail current [Rm]
const double deltax = 1.0;//expansion magnitudes of tail current sheet in x direction
const double deltay = 0.1;//expansion magnitudes of tail current sheet in y direction
const double SxD = 1.0;	//e-folding distance for the disk current sheet
const double SyD = 2.9;	//scale distance for the disk current sheet

/* Parameters relating to the vector potential, A, for the disk current */
const double f[5] = {59048.35734,  -135664.4246,  -913.4507339,   209989.1008,  -213142.9370};
const double b[5] = {19.69235037,  -18.16704312,   12.69175932,  -14.13692134,   14.13449724};
const double c[5] = {7.682813704,   9.663177797,  0.6465427021,   1.274059603,   1.280231032};


/* Shielding parameters for disk current field */
const double TailDiskA[36] = { -3.98467028e+02,  -1.14300168e+03,  -1.83630038e+03,  -7.39218042e+01,  -3.26398685e+02,  -2.99686811e+01,
							    -1.15703560e+03,  -6.04184603e+02,  -5.20487618e+01,  -2.03069124e+03,  -1.52912034e+03,  -6.38220995e+00,
							     2.58766603e+03,   2.13897918e+02,  -2.83022599e+01,   6.30130986e+02,   2.96855224e+03,   8.88632862e+02,
							     4.97386309e+02,   2.30425447e+03,   8.58417688e+02,   1.22695860e+03,   8.50168495e+02,  -2.09011094e+01,
							    -2.03918424e+02,  -7.92609902e+02,   1.11595569e+03,   5.27322683e+02,   2.24763404e+01,  -7.04405637e-02,
							    -1.40509314e+03,  -9.72040834e+01,   5.65673018e+00,  -1.38712910e+02,  -1.97975567e+03,   5.40760375e+00};
const double TailDiskP[6] = { 1.09108891,  0.67332998,  0.32667478,  0.95331615,  1.36276304,  0.00145152};


/***********************************************************************
 * NAME : 		double DiskThickness(x,y)
 * 
 * DESCRIPTION : 	Calculates the disk thickness in Rm.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 *
 * RETURNS :
 * 			double d	Current disk thickness (Rm).
 *
 * 
 * ********************************************************************/
double DiskThickness(double x, double y) {
	double y20 = y/20.0;
	return 7.0*d0 + 7.0*deltax*exp(x/7.0) + 7.0*deltay*y20*y20;
}

/***********************************************************************
 * NAME : 		void DiskField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the disk field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of magnetic field (in nT).
 * 			double *By	y-component of magnetic field (in nT).
 * 			double *Bz	z-component of magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void DiskField(	double x, double y, double z,
				double *Bx, double *By, double *Bz) {
	
	
	double xs, ys,zs;
	double d;
	double r;
	double dddx, dddy;
	double zeta, dzdx, dzdy, dzdz;
	double drdx, drdy, drdz;
	double S1, dS1drho, dS1dzeta, dS1dx, dS1dy, dS1dz;
	double S2, dS2drho, dS2dzeta, dS2dx, dS2dy, dS2dz;
	double As, dAsdS1, dAsdS2, dAsdx, dAsdy, dAsdz;
	double bpr, bmr, cpzeta, cpzeta2;
	double S1pS2, S1pS22, S1tS2;
	int i;
	
	/* init outputs */
	Bx[0] = 0.0;
	By[0] = 0.0;
	Bz[0] = 0.0;

	/* rescale inputs */
	xs = (x - xshift)*sc;
	ys = y*sc;
	zs = z*sc;

	/* current disk thickness */
	d = DiskThickness(xs,ys);

	dddy = deltay*7.0*ys*0.005;
	dddx = deltax*exp(xs/7.0);
	
	zeta = sqrt(zs*zs + d*d);
	dzdx = (d/zeta)*dddx;
	dzdy = (d/zeta)*dddy;
	dzdz = zs/zeta;
	

	r = sqrt(xs*xs + ys*ys);
	
	drdx = xs/r;
	drdy = ys/r;
	drdz = 0.0;
	
	if (isnan(drdx)) {
		drdx = 0.0;
		drdy = 0.0;
	}

	for (i=0;i<5;i++) {
		bpr = b[i] + r;
		bmr = b[i] - r;
		cpzeta = c[i] + zeta;
		cpzeta2 = cpzeta*cpzeta;
		
		S1 = sqrt(bpr*bpr + cpzeta2);
		S2 = sqrt(bmr*bmr + cpzeta2);
		
		dS1drho = bpr/S1;
		dS2drho = -bmr/S2;
	
		dS1dzeta = cpzeta/S1;
		dS2dzeta = cpzeta/S2;		


		dS1dx = dS1dzeta*dzdx + dS1drho*drdx;
		dS1dy = dS1dzeta*dzdy + dS1drho*drdy;
		dS1dz = dS1dzeta*dzdz + dS1drho*drdz;
	
		dS2dx = dS2dzeta*dzdx + dS2drho*drdx;
		dS2dy = dS2dzeta*dzdy + dS2drho*drdy;
		dS2dz = dS2dzeta*dzdz + dS2drho*drdz;	
		
		S1pS2 = S1 + S2;
		S1pS22 = S1pS2*S1pS2;
		S1tS2 = S1*S2;
		
		As = sqrt(S1pS22 - (4*b[i]*b[i]))/(S1tS2*S1pS22);
		
		dAsdS1 = (1.0/(S1tS2*S1pS2*sqrt(S1pS22 - 4*b[i]*b[i]))) - (As/S1pS22)/S1*(S2*S2 + S1*(3.0*S1 + 4.0*S2));
		dAsdS2 = (1.0/(S1tS2*S1pS2*sqrt(S1pS22 - 4*b[i]*b[i]))) - (As/S1pS22)/S2*(S1*S1 + S2*(3.0*S2 + 4.0*S1));		
		
		dAsdx = dAsdS1*dS1dx + dAsdS2*dS2dx;
		dAsdy = dAsdS1*dS1dy + dAsdS2*dS2dy;
		dAsdz = dAsdS1*dS1dz + dAsdS2*dS2dz;

		Bx[0] += -f[i]*xs*dAsdz;
		By[0] += -f[i]*ys*dAsdz;
		Bz[0] += f[i]*(2*As + ys*dAsdy + xs*dAsdx);

	}

}


/***********************************************************************
 * NAME : 		void DiskShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the magnetopause shielding field required
 * 					to contain the disk field.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of magnetopause field (in nT).
 * 			double *By	y-component of magnetopause field (in nT).
 * 			double *Bz	z-component of magnetopause field (in nT).
 *
 * 
 * ********************************************************************/
void DiskShield(	double x, double y, double z,
					double *Bx, double *By, double *Bz) {
	/* this calculates the shielding field using equation 12 of 
	 * Korth et al 2015 */
	ShieldField(x,y,z,6,TailDiskA,TailDiskP,Bx,By,Bz);
}


/***********************************************************************
 * NAME : 		void Disk(x,y,z,t1,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's disk field at a given 
 * 					postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 * 			double t1	Scale factor for the disk field contribution to 
 * 						the model.
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's disk field (in nT).
 * 			double *By	y-component of Mercury's disk field (in nT).
 * 			double *Bz	z-component of Mercury's disk field (in nT).
 *
 * 
 * ********************************************************************/
void Disk(	double x, double y, double z, double t1,
			double *Bx, double *By, double *Bz) {
	
	double Fx, Fy, Fz, Sx, Sy, Sz;
	
	/* get disk field and its corresponding shielding field, then add */
	DiskField(x,y,z,&Fx,&Fy,&Fz);
	DiskShield(x,y,z,&Sx,&Sy,&Sz);
	Bx[0] = t1*(Fx - Sx);
	By[0] = t1*(Fy - Sy);
	Bz[0] = t1*(Fz - Sz);
	
}

