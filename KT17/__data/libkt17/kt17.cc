#include "kt17.h"

/* Mercury Magnetic Field Model (Korth et al. 2015)
 * 
 * Written By Matt James
 */

/* Some model parameters: */
const double mu = -190.0; //dipole moment [nT Rm^3]
const double dz = 0.196;	//dipole z offset
const double Rss = 1.41;	//distance from Mercury center to sub-solar magnetopause [Rm]
const double R0 = 1.42;	//distance to fitted subsolar magnetopause
const double t1 = 7.64;	//tail disk current magnitude
const double t2 = 2.06;	//harris sheet current magntiude
const double d0 = 0.09;	//half-width of current sheet in Z at inner edge of tail current [Rm]
const double deltax = 1.0;//expansion magnitudes of tail current sheet in x direction
const double deltay = 0.1;//expansion magnitudes of tail current sheet in y direction
const double SxD = 1.0;	//e-folding distance for the disk current sheet
const double SyD = 2.9;	//scale distance for the disk current sheet
const double SxQH = 1.5;	//e-folding distance for the sunward expansion of the harris sheet
const double SyQH = 9.0;	//scale distance for the flankward expansion of the harris sheet
const double zshift = 3.5;//location of image sheets for harris sheet
const double mptol = 0.001;

/* Parameters relating to the vector potential, A, for the disk current */
const double f[5] = {59048.35734,  -135664.4246,  -913.4507339,   209989.1008,  -213142.9370};
const double b[5] = {19.69235037,  -18.16704312,   12.69175932,  -14.13692134,   14.13449724};
const double c[5] = {7.682813704,   9.663177797,  0.6465427021,   1.274059603,   1.280231032};

/* Shielding parameters for disk current field */
const double TailDiskA[6][6] = {{ -3.98467028e+02,  -1.14300168e+03,  -1.83630038e+03,  -7.39218042e+01,  -3.26398685e+02,  -2.99686811e+01},
							   { -1.15703560e+03,  -6.04184603e+02,  -5.20487618e+01,  -2.03069124e+03,  -1.52912034e+03,  -6.38220995e+00},
							   {  2.58766603e+03,   2.13897918e+02,  -2.83022599e+01,   6.30130986e+02,   2.96855224e+03,   8.88632862e+02},
							   {  4.97386309e+02,   2.30425447e+03,   8.58417688e+02,   1.22695860e+03,   8.50168495e+02,  -2.09011094e+01},
							   { -2.03918424e+02,  -7.92609902e+02,   1.11595569e+03,   5.27322683e+02,   2.24763404e+01,  -7.04405637e-02},
							   { -1.40509314e+03,  -9.72040834e+01,   5.65673018e+00,  -1.38712910e+02,  -1.97975567e+03,   5.40760375e+00}};
const double TailDiskP[6] = { 1.09108891,  0.67332998,  0.32667478,  0.95331615,  1.36276304,  0.00145152};

/* Shielding parameters for the Quasi-Harris sheet */
const double TailQHA[6][6] = {{  -91.67686636,   -87.31240824,   251.8848107 ,    95.65629983,   -80.968107  ,   198.1447476 },
							{ -283.1968987 ,  -269.1514899 ,   504.632231  ,   166.027215  ,  -214.9025413 ,   623.7920115 },
							{  -35.99544615,  -322.864469  ,   345.710579  ,   928.8553184 ,   810.129509  ,    19.62627762},
							{  -12.70326428,   490.4662048 ,  -814.0985363 , -1781.184984  , -1371.261326  ,    60.3136479 },
							{  116.630551  ,  -178.3347065 ,   604.0308838 ,  1155.151174  ,   770.3896601 ,  -202.8545948 },
							{  298.6337705 ,   304.7964641 ,    33.70850254,   393.6080147 ,   308.1194271 ,  -660.1691658 }};
const double TailQHP[6] = { 1.67762971,  1.29222658,  0.31162534, -0.43926691,  0.75780748,  1.49777952};

/* Shielding parameters for the dipole field */
const double DipA[4][4] = {{    7.79240768,    74.37654983,     4.11964707,  -131.33086   },
							{  546.6006311 , -1077.694401  ,    52.46268495,  1057.273707  },
							{  -74.91550119,  -141.8047123 ,     3.87600489,   156.2250932 },
							{ -506.6470185 ,  1439.804381  ,   -64.55225925, -1443.754088  }};
const double DipP[4] = { 0.14122971,  0.74398476,  1.04279834,  0.7057116 };

bool WithinMP(double x, double y, double z, double Rsm) {
	double r, yz, theta, r_mp;
	r = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z,2.0));
	yz = sqrt(pow(y,2.0) + pow(z,2.0));
	theta = atan2(yz,x);
	r_mp = (Rsm+mptol)*sqrt(2/(1 + cos(theta)));
	
	if (r_mp >= r) {
		return true;
	} else {
		return false;
	} 
}

/*
void kt17DipoleField(double x, double y, double z, double *Bx, double *By, double *Bz) {
	double x2,y2,z2,R,q,v;
	x2 = x*x;
	y2 = y*y;
	z2 = z*z;
	R = sqrt(x2 + y2 + z2);
	v = 3.0*x*z;
	q = mu/pow(R,2.5);
	*Bx = -q*v;
	*By = -3.0*y*q*z;
	*Bz = q*(x2 + y2 - 2.0*z2);
}*/

void kt17DipoleField(double x, double y, double z, double *Bx, double *By, double *Bz) {
	double M = mu;
	double mhat[3] = {0.0/M, 0.0/M, mu/M};
	double r = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z,2.0));
	double rhat[3] = {x/r, y/r, z/r};
	double mdotr = (mhat[0]*rhat[0] + mhat[1]*rhat[1] + mhat[2]*rhat[2]);
	*Bx = M*((3*mdotr*rhat[0] - mhat[0])/pow(r,3));
	*By = M*((3*mdotr*rhat[1] - mhat[1])/pow(r,3));
	*Bz = M*((3*mdotr*rhat[2] - mhat[2])/pow(r,3));
	return;
}

void kt17DipoleShield(double x, double y, double z, double *Bx, double *By, double *Bz) {
	*Bx = 0.0;
	*By = 0.0;
	*Bz = 0.0;
	int i, k;
	double pik;
	for (i=0;i<4;i++) {
		for (k=0;k<4;k++) {
			pik = sqrt(pow(DipP[i],2.0) + pow(DipP[k],2.0));
			*Bx+= DipA[i][k]*pik*exp(pik*x)*cos(DipP[i]*y)*sin(DipP[k]*z);
			*By+=-DipA[i][k]*exp(pik*x)*DipP[i]*sin(DipP[i]*y)*sin(DipP[k]*z);
			*Bz+= DipA[i][k]*exp(pik*x)*cos(DipP[i]*y)*DipP[k]*cos(DipP[k]*z);
		}
	}
	return;
}

void kt17DipoleB(double x, double y, double z, double *Bx, double *By, double *Bz) {
	double Fx, Fy, Fz, Hx, Hy, Hz;
	kt17DipoleField(x,y,z,&Fx,&Fy,&Fz);
	kt17DipoleShield(x,y,z,&Hx,&Hy,&Hz);
	*Bx = Fx - Hx;
	*By = Fy - Hy;
	*Bz = Fz - Hz;

	
}

double kt17DiskThickness(double x, double y) {
	return 7.0*d0 + 7.0*deltax*exp(x/7.0) + 7.0*deltay*pow((y/20.0),2.0);
}

void kt17DiskField(double xin, double yin, double zin, double *Bx, double *By, double *Bz) {
	
	const double xshift = 0.3, sc = 7.0;
	/*Big-ass fudge here, not even defined all of the variables yet and
	 * there is something not mentioned in the paper! Shifting in the X 
	 * direction and then multiplying by 7...and it magically works!*/
	double x=(xin-xshift)*sc, y = yin*sc, z = zin*sc;
	double d = kt17DiskThickness(x,y);
	double r;
	double dddx, dddy;
	double zeta, dzdx, dzdy, dzdz;
	double drdx, drdy, drdz;
	double S1, dS1drho, dS1dzeta, dS1dx, dS1dy, dS1dz;
	double S2, dS2drho, dS2dzeta, dS2dx, dS2dy, dS2dz;
	double As, dAsdS1, dAsdS2, dAsdx, dAsdy, dAsdz;
	int i;
	
	*Bx = 0.0;
	*By = 0.0;
	*Bz = 0.0;
	
	//dddx = (deltax/SxD)*exp(x/SxD);
	//dddy = (2.0*deltay*y)/(pow(SyD,2.0));
	dddy=deltay*7.0*y*0.005;
	dddx=deltax*exp(x/7.0);
	
	zeta = sqrt(pow(z,2.0) + pow(d,2.0));
	dzdx = (d/zeta)*dddx;
	dzdy = (d/zeta)*dddy;
	dzdz = z/zeta;
	

	r = sqrt(pow(x,2.0) + pow(y,2.0));
	
	drdx = x/r;
	drdy = y/r;
	drdz = 0.0;
	
	if (isnan(drdx)) {
		drdx = 0.0;
		drdy = 0.0;
	}
	
	for (i=0;i<5;i++) {
		S1 = sqrt(pow((b[i] + r),2.0) + pow((c[i] + zeta),2.0));
		S2 = sqrt(pow((b[i] - r),2.0) + pow((c[i] + zeta),2.0));
		
		dS1drho = (b[i] + r)/S1;
		dS2drho = -(b[i] - r)/S2;
	
		dS1dzeta = (c[i] + zeta)/S1;
		dS2dzeta = (c[i] + zeta)/S2;		


		dS1dx = dS1dzeta*dzdx + dS1drho*drdx;
		dS1dy = dS1dzeta*dzdy + dS1drho*drdy;
		dS1dz = dS1dzeta*dzdz + dS1drho*drdz;
	
		dS2dx = dS2dzeta*dzdx + dS2drho*drdx;
		dS2dy = dS2dzeta*dzdy + dS2drho*drdy;
		dS2dz = dS2dzeta*dzdz + dS2drho*drdz;	
		
		As = sqrt(pow((S1 + S2),2.0) - pow((2*b[i]),2.0))/(S1*S2*pow((S1 + S2),2.0));
		
		dAsdS1 = (1.0/(S1*S2*(S1 + S2)*sqrt(pow((S1 + S2),2.0) - pow((2*b[i]),2.0)))) - (As/pow((S1 + S2),2.0))/S1*(pow(S2,2.0) + S1*(3.0*S1 + 4.0*S2));
		dAsdS2 = (1.0/(S1*S2*(S1 + S2)*sqrt(pow((S1 + S2),2.0) - pow((2*b[i]),2.0)))) - (As/pow((S1 + S2),2.0))/S2*(pow(S1,2.0) + S2*(3.0*S2 + 4.0*S1));		
		
		
		dAsdx = dAsdS1*dS1dx + dAsdS2*dS2dx;
		dAsdy = dAsdS1*dS1dy + dAsdS2*dS2dy;
		dAsdz = dAsdS1*dS1dz + dAsdS2*dS2dz;

		*Bx+= -f[i]*x*dAsdz;
		*By+= -f[i]*y*dAsdz;
		*Bz+= f[i]*(2*As + y*dAsdy + x*dAsdx);

	}
	
	return;
}


void kt17DiskShield(double x, double y, double z, double *Bx, double *By, double *Bz) {
	*Bx = 0.0;
	*By = 0.0;
	*Bz = 0.0;
	int i, k;
	double pik;
	for (i=0;i<6;i++) {
		for (k=0;k<6;k++) {
			pik = sqrt(pow(TailDiskP[i],2.0) + pow(TailDiskP[k],2.0));
			*Bx+= TailDiskA[i][k]*pik*exp(pik*x)*cos(TailDiskP[i]*y)*sin(TailDiskP[k]*z);
			*By+=-TailDiskA[i][k]*exp(pik*x)*TailDiskP[i]*sin(TailDiskP[i]*y)*sin(TailDiskP[k]*z);
			*Bz+= TailDiskA[i][k]*exp(pik*x)*cos(TailDiskP[i]*y)*TailDiskP[k]*cos(TailDiskP[k]*z);
		}
	}
	return;
}

void kt17DiskB(double x, double y, double z, double *Bx, double *By, double *Bz, double T1) {
	double Fx, Fy, Fz, Hx, Hy, Hz;
	kt17DiskField(x,y,z,&Fx,&Fy,&Fz);
	kt17DiskShield(x,y,z,&Hx,&Hy,&Hz);
	*Bx = T1*(Fx - Hx);
	*By = T1*(Fy - Hy);
	*Bz = T1*(Fz - Hz);
}


double kt17QuasiHarrisThickness(double x, double y) {
	return d0 + deltax*exp(x/SxQH) + deltay*pow((y/SyQH),2.0);
}

void kt17QuasiHarrisField(double x, double y, double z, double *Bx, double *By, double *Bz) {
	double d = kt17QuasiHarrisThickness(x,y);
	*Bx = (2.0/d)*tanh(z/d);
	*By = 0.0;
	*Bz = ((z*2.0)/pow(d,2.0))*tanh(z/d)*((deltax/SxQH)*exp(x/(SxQH)));
	return;
}

void kt17QuasiHarrisShield(double x, double y, double z, double *Bx, double *By, double *Bz) {
	*Bx = 0.0;
	*By = 0.0;
	*Bz = 0.0;
	int i, k;
	double pik;
	for (i=0;i<6;i++) {
		for (k=0;k<6;k++) {
			pik = sqrt(pow(TailQHP[i],2.0) + pow(TailQHP[k],2.0));
			*Bx+= TailQHA[i][k]*pik*exp(pik*x)*cos(TailQHP[i]*y)*sin(TailQHP[k]*z);
			*By+=-TailQHA[i][k]*exp(pik*x)*TailQHP[i]*sin(TailQHP[i]*y)*sin(TailQHP[k]*z);
			*Bz+= TailQHA[i][k]*exp(pik*x)*cos(TailQHP[i]*y)*TailQHP[k]*cos(TailQHP[k]*z);
		}
	}
	return;
}

void kt17QuasiHarrisB(double x, double y, double z, double *Bx, double *By, double *Bz, double T2) {
	double Fx, Fy, Fz, Hx, Hy, Hz;
	double F1x, F1y, F1z, F2x, F2y, F2z;
	kt17QuasiHarrisField(x,y,z,&Fx,&Fy,&Fz);
	kt17QuasiHarrisField(x,y,z+zshift,&F1x,&F1y,&F1z);
	kt17QuasiHarrisField(x,y,z-zshift,&F2x,&F2y,&F2z);
	kt17QuasiHarrisShield(x,y,z,&Hx,&Hy,&Hz);
	/* A fudge here, apparently the image current sheets must be 
	 * divided by 2 before being added to the rest. Then everything
	 * by 2.*/ 
	*Bx = 0.5*T2*(Fx - Hx + 0.5*(F1x + F2x));
	*By = 0.5*T2*(Fy - Hy + 0.5*(F1y + F2y));
	*Bz = 0.5*T2*(Fz - Hz - 0.5*(F1z + F2z));

}


void kt17B(double xin, double yin, double zin, double *Bx, double *By, double *Bz, double Rsm, double T1, double T2) {
	double k = R0/Rsm;
	double Bix, Biy, Biz, Bdx, Bdy, Bdz, Bqx, Bqy, Bqz;
	double x = xin*k, y = yin*k, z = zin*k;
	kt17DipoleB(x,y,z,&Bix,&Biy,&Biz);
	kt17DiskB(x,y,z,&Bdx,&Bdy,&Bdz,T1);
	kt17QuasiHarrisB(x,y,z,&Bqx,&Bqy,&Bqz,T2);
	
	*Bx = pow(k,3.0)*Bix + Bdx + Bqx;
	*By = pow(k,3.0)*Biy + Bdy + Bqy;
	*Bz = pow(k,3.0)*Biz + Bdz + Bqz;
	
}

/*void kt17Barray(int n, double *xin, double *yin, double *zin, double *Bx, double *By, double *Bz, double *Rsm, double *T1, double *T2) {
	int i;
	for (i=0;i<n;i++) {
		kt17B(xin[i],yin[i],zin[i],&Bx[i],&By[i],&Bz[i],Rsm[i],T1[i],T2[i]);
	}
}*/

void kt17Barray(int n, double *xin, double *yin, double *zin, double *Bx, double *By, double *Bz, int nP, int Plen, double *Params) {
	int i, pi;
	for (i=0;i<n;i++) {
		pi = (i % Plen/nP)*nP; //position in the Params array to start
		if (nP == 3) {
			kt17B(xin[i],yin[i],zin[i],&Bx[i],&By[i],&Bz[i],Params[pi],Params[pi+1],Params[pi+2]);
		} else {
			kt17(xin[i],yin[i],zin[i],&Bx[i],&By[i],&Bz[i],Params[pi],Params[pi+1]);
		}
	}
}

void kt17(double x, double y, double z, double *Bx, double *By, double *Bz, double Rsun, double DistIndex) {
	double f, Rsm, T1, T2;
	f = 2.06873 - 0.00279*DistIndex;
	Rsm = f*pow(Rsun,1.0/3.0);
	T1 = 6.495 + 0.0229*DistIndex;
	T2 = 1.6245 + 0.0088*DistIndex;
	kt17B(x,y,z,Bx,By,Bz,Rsm,T1,T2);
}
/*
void kt17array(int n, double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Rsun, double *DistIndex) {
	int i;
	for (i=0;i<n;i++) {
		kt17(x[i],y[i],z[i],&Bx[i],&By[i],&Bz[i],Rsun[i],DistIndex[i]);
	}
}*/
