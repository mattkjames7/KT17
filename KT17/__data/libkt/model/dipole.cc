#include "dipole.h"

/* dipole moment [nT Rm^3] */
const double mu = -190.0; 

/* Shielding parameters for the dipole field */
const double DipA[16] =  {    7.79240768,    74.37654983,     4.11964707,  -131.33086   ,
							  546.6006311 , -1077.694401  ,    52.46268495,  1057.273707  ,
							  -74.91550119,  -141.8047123 ,     3.87600489,   156.2250932 ,
							 -506.6470185 ,  1439.804381  ,   -64.55225925, -1443.754088  };
const double DipP[4] = { 0.14122971,  0.74398476,  1.04279834,  0.7057116 };


/***********************************************************************
 * NAME : 		void DipoleField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's dipole field at a given 
 * 					postion.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's magnetic field (in nT).
 * 			double *By	y-component of Mercury's magnetic field (in nT).
 * 			double *Bz	z-component of Mercury's magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void DipoleField(	double x, double y, double z,
					double *Bx, double *By, double *Bz) {
	
	/*calculate r and its unit vector rhat */
	double r, r2, r3;
	double rhatx, rhaty, rhatz;	
	r2 = x*x + y*y + z*z;
	r = sqrt(r2);
	r3 = r*r2;
	rhatx = x/r;
	rhaty = y/r;
	rhatz = z/r;
	
	/* this is the dot product of the magnetic moment vector with rhat 
	 * NOTE: dipole is aligned with z axis, so mhat = [0.0,0.0,1.0]
	 * therefore mhat.rhat = rhatz*/	
	double mdotr = rhatz;		
	
	/* times mdotr by 3 (this is done for each component in the old code,
	 * might as well save time by doing it here once)*/
	double mdotr3 = mdotr*3.0;
	
	/* another constant M/r**3 used for each component */
	double M_r3 = mu/r3;

	/* now calculate each component B = M*(3*mhat.rhat - mhat)/r**3 */
	Bx[0] = M_r3*mdotr3*rhatx;
	By[0] = M_r3*mdotr3*rhaty;
	Bz[0] = M_r3*(mdotr3*rhatz - 1.0);
}

/***********************************************************************
 * NAME : 		void DipoleShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the magnetopause shielding field required
 * 					to contain the dipole field.
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
void DipoleShield(	double x, double y, double z,
					double *Bx, double *By, double *Bz) {
	/* this calculates the shielding field using equation 12 of 
	 * Korth et al 2015 */
	ShieldField(x,y,z,4,DipA,DipP,Bx,By,Bz);
}


/***********************************************************************
 * NAME : 		void Dipole(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's dipole field at a given 
 * 					postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's magnetic field (in nT).
 * 			double *By	y-component of Mercury's magnetic field (in nT).
 * 			double *Bz	z-component of Mercury's magnetic field (in nT).
 *
 * 
 * ********************************************************************/
void Dipole(	double x, double y, double z, 
				double *Bx, double *By, double *Bz) {
	
	double Fx, Fy, Fz, Sx, Sy, Sz;
	
	/* get dipole field and its corresponding shielding field, then add */
	DipoleField(x,y,z,&Fx,&Fy,&Fz);
	DipoleShield(x,y,z,&Sx,&Sy,&Sz);
	Bx[0] = Fx - Sx;
	By[0] = Fy - Sy;
	Bz[0] = Fz - Sz;

	
}

