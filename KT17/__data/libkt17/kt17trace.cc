#include "kt17trace.h"

void UnitVector(double ix, double iy, double iz, double *ox, double *oy, double *oz) {
	double I = sqrt(pow(ix, 2.0) + pow(iy,2.0) + pow(iz,2.0));
	*ox = ix/I;
	*oy = iy/I;
	*oz = iz/I;
	return;
}

void ReverseElements(double x[], int n) {
	double tmp[n];
	int i;
	for (i=0;i<n;i++) {
		tmp[i] = x[i];
	}
	for (i=0;i<n;i++) {
		x[i] = tmp[n-i-1];
	}
	return;
}

void GetBits(const int x, bool *bits) {
	int i;
	for (i=0;i<32;++i) {
		bits[i] = x & (1 << i);
	}
}

bool InsideTraceLims(double x, double y, double z, bool *bits, int dir, double Rsm) {
	double Rmsm = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z,2.0));
	double Rmso = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z+dz,2.0));

	bool 	WMP = true, 	//Within magnetopause
			OSP = true,		//Outside the planet
			InTail = true,	//Within 10 Rm
			InBox = true,	//Inside a box -6<x<2, -4<y<4, -4<z<4
			OSD = true,		//Outside the dipole-centered planet
			OSPcore = true, //Outside core of planet (2030km = 0.832Rm)
			OSDcore = true;	//Outside core of dipole
	
	/* check if within MP */
	if (bits[0]) {
		WMP = WithinMP(x,y,z,Rsm);
	}
	
	if (bits[1] && bits[2]) {
		/* check both planet and dipole */
		if (dir > 0) {
			if (Rmso > 1.0 || z < 0) {
				OSP = true;
				OSD = true;
			} else {
				OSP = false;
				OSD = false;
			}
		} else {
			if (Rmsm > 1.0 || z > 0) {
				OSP = true;
				OSD = true;
			} else {
				OSP = false;
				OSD = false;
			}		
		}
	} else if (bits[1] && !bits[2]) {
		/* Only check actual planet*/
		OSP = Rmso > 1.0;
	} else if (!bits[1] && bits[2]) {
		/* Only check if outside planets radius of the dipole */
		OSD = Rmsm > 1.0;
	}
	
	if (bits[3]) {
		/* Check if within 10 Rm of Mercury */
		InTail = Rmsm < 10.0;
	}
	
	if (bits[4] && bits[5]) {
		/* check both planet and dipole */
		if (dir > 0) {
			if (Rmso > 0.832 || z < 0) {
				OSPcore = true;
				OSDcore = true;
			} else {
				OSPcore = false;
				OSDcore = false;
			}
		} else {
			if (Rmsm > 0.832 || z > 0) {
				OSPcore = true;
				OSDcore = true;
			} else {
				OSPcore = false;
				OSDcore = false;
			}		
		}
	} else if (bits[4] && !bits[5]) {
		/* Only check actual planet*/
		OSPcore = Rmso > 0.832;
	} else if (!bits[4] && bits[5]) {
		/* Only check if outside planets radius of the dipole */
		OSDcore = Rmsm > 0.832;
	}
		
	/* check if within the box*/
	if (bits[6]) {
		if ((x >= -6) && (x <= 2) && (y >= -4) && (y <= 4) && (z >= -4) && (z <= 4)) {
			InBox = true;
		} else {
			InBox = false;
		}
	}	
	/* combine the results */
	return 	(WMP && InBox && InTail && OSP && OSD && OSPcore && OSDcore);
			
	
}

//bool InsideTraceLims(double x, double y, double z, int LimType, int dir, double Rsm) {
	//bool WMP = WithinMP(x,y,z,Rsm), OsP;
	//double Rmsm = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z,2.0));
	//double Rmso = sqrt(pow(x,2.0) + pow(y,2.0) + pow(z+dz,2.0));
	//if (dir > 0) {
		//if (Rmso > 1.0 || z < 0) {
			//OsP = true;
		//} else {
			//OsP = false;
		//}
	//} else {
		//if (Rmsm > 1.0 || z > 0) {
			//OsP = true;
		//} else {
			//OsP = false;
		//}		
	//}
	//switch (LimType) {
		//case 1:
			///* Confine to a box */
			//if ((x > -6 && x < 2) && (y > -4 && y < 4) && (z > -4 && z < 4)) {
				//return true;
			//} else { 
				//return false;
			//} 
			//break;
		//case 2:
			///* Confine to a box, also terminating at the planet */
			//if ((x > -6 && x < 2) && (y > -4 && y < 4) && (z > -4 && z < 4) && OsP == true) {
				//return true;
			//} else { 
				//return false;
			//} 
			//break;
		//case 3:
			///* Confine to just outside the planet */
			//if (OsP == true) {
				//return true;
			//} else {
				//return false;
			//}
			//break;
		//case 4:
			///* trace to MP and to planet, stop within 10 Rm */
			//if (OsP == true && WMP == true && Rmsm < 10.0) {
				//return true;
			//} else {
				//return false;
			//}
			//break;
		//default:
			///* Default action is to trace to MP and to planet */
			//if (OsP == true && WMP == true) {
				//return true;
			//} else {
				//return false;
			//}			
			//break;
	//}
//}

void NorthSouthFLs(double flx[],double fly[],double flz[], double Rmsm[], double Rmso[], int N, double **Nflx, double **Nfly, double **Nflz, double **NRmsm, double **NRmso, int *nn, double **Sflx, double **Sfly, double **Sflz, double **SRmsm, double **SRmso, int *ns) {
	int i,cn = 0, cs = 0;
	while (flz[cn] >= 0 && isfinite(flz[cn]) && cn < N) {
		cn++;
	}

	(*nn)=cn;
	if (cn > 0) {
		*Nflx = (double*) malloc(cn*sizeof(double));
		*Nfly = (double*) malloc(cn*sizeof(double));
		*Nflz = (double*) malloc(cn*sizeof(double));
		*NRmsm = (double*) malloc(cn*sizeof(double));
		*NRmso = (double*) malloc(cn*sizeof(double));
		for (i=0;i<cn;i++) {
			(*Nflx)[i]=flx[i];
			(*Nfly)[i]=fly[i];
			(*Nflz)[i]=flz[i];
			(*NRmso)[i]=Rmso[i];
			(*NRmsm)[i]=Rmsm[i];
		}
	} else { 
		*Nflx=NULL;
		*Nfly=NULL;
		*Nflz=NULL;
	}
	
	i=cn+1;
	while (flz[i] < 0 && isfinite(flz[i]) && i < N) {
		cs++;
		i++;
	}
	(*ns)=cs;
	if (cs > 0) {
		*Sflx = (double*) malloc(cs*sizeof(double));
		*Sfly = (double*) malloc(cs*sizeof(double));
		*Sflz = (double*) malloc(cs*sizeof(double));
		*SRmsm = (double*) malloc(cs*sizeof(double));
		*SRmso = (double*) malloc(cs*sizeof(double));
		for (i=0;i<cs;i++) {
			(*Sflx)[i]=flx[(cn+cs-1)-i];
			(*Sfly)[i]=fly[(cn+cs-1)-i];
			(*Sflz)[i]=flz[(cn+cs-1)-i];
			(*SRmso)[i]=Rmso[(cn+cs-1)-i];
			(*SRmsm)[i]=Rmsm[(cn+cs-1)-i];
		}
	} else { 
		*Sflx=NULL;
		*Sfly=NULL;
		*Sflz=NULL;
	}	
		
}

double linterp(double x0, double x1, double y0, double y1, const double xt) {
	double m;
	m = (y1 - y0)/(x1 - x0);
	return m*(xt - x0) + y0;
}

void EqFootprint(double *Nflx, double *Nfly, double *Nflz, int nN, double *Sflx, double *Sfly, double *Sflz, int nS, double *lshell, double *mlte) {
	//printf("nN:%d nS:%d \n",nN,nS);
	if (nN > 0 && nS > 0) {
	//	printf("nN:%d nS:%d \n",nN,nS);
		double fpx, fpy;
		fpx = linterp(Sflz[nS-1],Nflz[nN-1],Nflx[nN-1],Sflx[nS-1],0.0);
		fpy = linterp(Sflz[nS-1],Nflz[nN-1],Nfly[nN-1],Sfly[nS-1],0.0);
		*lshell = sqrt(pow(fpx,2.0) + pow(fpy,2.0));
		*mlte = fmod((M_PI + atan2(fpy,fpx)) , (2*M_PI)) * (180.0/M_PI)/15.0;
	//	printf("%f %f \n",*lshell,*mlte);
	} else { 
		*lshell = NAN;
		*mlte = NAN;
	}
	return;
}

double min(double r, double *a, int n) {
	int i;
	double mn = fabs(r-a[0]);
	double dr;
	for (i=1;i<n;i++) {
		dr = fabs(r-a[i]);
		if (dr < mn) {
			mn = dr;
		}
	}
	return mn;
}

void PlanetFootprints(double *flx, double *fly, double *flz, double *Rmsm, double *Rmso, int n,
					double *mlat, double *mlt, double *lat, double *lct, double *mlatc, double *mltc, double *latc, double *lctc) {
	double fpx, fpy, fpz, mnmsm, mnmso, mnmsmc, mnmsoc;
	int i;
	//These are the minimum distances to the surface and the core for the actual planet and the dipole-centered planet
	mnmsm = min(1.0,Rmsm,n);
	mnmso = min(1.0,Rmso,n);
	mnmsmc = min(0.832,Rmsm,n);
	mnmsoc = min(0.832,Rmso,n);
	
	for (i=0;i<n-1;i++) {
		if ((Rmsm[i] <= 0.832 && Rmsm[i+1] > 0.832) || (fabs(0.832-Rmsm[i]) == mnmsmc && fabs(0.832-Rmsm[i]) < 0.1)){
			fpx = linterp(Rmsm[i],Rmsm[i+1],flx[i],flx[i+1],0.832);
			fpy = linterp(Rmsm[i],Rmsm[i+1],fly[i],fly[i+1],0.832);
			fpz = linterp(Rmsm[i],Rmsm[i+1],flz[i],flz[i+1],0.832);
			*mlatc = fmod(atan2(fpz, sqrt(pow(fpx,2.0) + pow(fpy,2.0))),2*M_PI) * (180.0/M_PI); 
			*mltc = fmod((M_PI + atan2(fpy,fpx)) , (2*M_PI)) * (180.0/M_PI)/15.0;
		}

		if ((Rmsm[i] <= 1.0 && Rmsm[i+1] > 1.0) || (fabs(1-Rmsm[i]) == mnmsm && fabs(1-Rmsm[i]) < 0.1)){
			fpx = linterp(Rmsm[i],Rmsm[i+1],flx[i],flx[i+1],1.0);
			fpy = linterp(Rmsm[i],Rmsm[i+1],fly[i],fly[i+1],1.0);
			fpz = linterp(Rmsm[i],Rmsm[i+1],flz[i],flz[i+1],1.0);
			*mlat = fmod(atan2(fpz, sqrt(pow(fpx,2.0) + pow(fpy,2.0))),2*M_PI) * (180.0/M_PI); 
			*mlt = fmod((M_PI + atan2(fpy,fpx)) , (2*M_PI)) * (180.0/M_PI)/15.0;
			break;
		}
	}
	for (i=0;i<n-1;i++) {
		if ((Rmso[i] <= 0.832 && Rmso[i+1] > 0.832) || (fabs(0.832-Rmso[i]) == mnmsoc && fabs(0.832-Rmso[i]) < 0.1)){
			fpx = linterp(Rmso[i],Rmso[i+1],flx[i],flx[i+1],0.832);
			fpy = linterp(Rmso[i],Rmso[i+1],fly[i],fly[i+1],0.832);
			fpz = linterp(Rmso[i],Rmso[i+1],flz[i]+dz,flz[i+1]+dz,0.832);
			*latc = fmod(atan2(fpz, sqrt(pow(fpx,2.0) + pow(fpy,2.0))),2*M_PI) * (180.0/M_PI); 
			*lctc = fmod((M_PI + atan2(fpy,fpx)) , (2*M_PI)) * (180.0/M_PI)/15.0;
		}
		
		if ((Rmso[i] <= 1.0 && Rmso[i+1] > 1.0) || (fabs(1-Rmso[i]) == mnmso && fabs(1-Rmso[i]) < 0.2)) {
			fpx = linterp(Rmso[i],Rmso[i+1],flx[i],flx[i+1],1.0);
			fpy = linterp(Rmso[i],Rmso[i+1],fly[i],fly[i+1],1.0);
			fpz = linterp(Rmso[i],Rmso[i+1],flz[i]+dz,flz[i+1]+dz,1.0);
			*lat = fmod(atan2(fpz, sqrt(pow(fpx,2.0) + pow(fpy,2.0))),2*M_PI) * (180.0/M_PI); 
			*lct = fmod((M_PI + atan2(fpy,fpx)) , (2*M_PI)) * (180.0/M_PI)/15.0;
			break;
		}
		

	}
}

void FieldLineLength(double *x, double *y, double *z, int n, double *r, double *len, double *lenc) {
	len[0] = 0.0; //field line length outside the planet
	lenc[0] = 0.0; // field line length outside the core
	int i;
	for (i=0;i<n-1;i++) {
		if (r[i] >= 1.0 || r[i+1] >= 1.0) {
			len[0] += sqrt(pow(x[i+1]-x[i],2.0) + pow(y[i+1]-y[i],2.0) + pow(z[i+1]-z[i],2.0));
		} 
		if (r[i] >= 0.832 || r[i+1] >= 0.832) {
			lenc[0] += sqrt(pow(x[i+1]-x[i],2.0) + pow(y[i+1]-y[i],2.0) + pow(z[i+1]-z[i],2.0));
		} 
	}
	


}

void Rvecs(double x0, double y0, double z0, double *rx, double *ry, double *rz, double Rsm, double T1, double T2, double step3) {
	double bx, by, bz, bm;
	kt17B(x0,y0,z0,&bx,&by,&bz,Rsm,T1,T2);
	bm = step3/sqrt(bx*bx + by*by + bz*bz);
	(*rx) = bx*bm;
	(*ry) = by*bm;
	(*rz) = bz*bm;
	return;
}


void kt17StepRKM(double x0, double y0, double z0, double *step, double maxstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double Rsm, double T1, double T2) {
	double r11,r12,r13,r21,r22,r23,r31,r32,r33,r41,r42,r43,r51,r52,r53,step3=(*step)/3.0, Err;
	bool repeat=true;
	//int i=0;
	while (repeat) {
		Rvecs(x0,y0,z0,&r11,&r12,&r13,Rsm,T1,T2,step3);
		Rvecs(x0+r11,y0+r12,z0+r13,&r21,&r22,&r23,Rsm,T1,T2,step3);	
		Rvecs(x0+0.5*(r11+r21),y0+0.5*(r12+r22),z0+0.5*(r13+r23),&r31,&r32,&r33,Rsm,T1,T2,step3);
		Rvecs(x0+0.375*(r11+3*r31),y0+0.375*(r12+3*r32),z0+0.375*(r13+3*r33),&r41,&r42,&r43,Rsm,T1,T2,step3);
		Rvecs(x0+1.5*(r11-3*r31+4*r41),y0+1.5*(r12-3*r32+4*r42),z0+1.5*(r13-3*r33+4*r43),&r51,&r52,&r53,Rsm,T1,T2,step3);
	
		Err = fabs(r11 - 4.5*r31 + 4*r41 - 0.5*r51) + fabs(r12 - 4.5*r32 + 4*r42 - 0.5*r52) + fabs(r13 - 4.5*r33 + 4*r43 - 0.5*r53);
		
		if ((Err <= ErrMax) && (fabs(*step) <= maxstep)) {
			repeat = false;
		} else {
			if (Err > ErrMax) {
				if (*step > minstep) {
					(*step) = (*step)*0.5;
				}else{
					repeat = false;
				}
			}
			if (fabs(*step) > maxstep) {
				(*step) = maxstep;
			}
		}
		
		if ((Err < 0.04*ErrMax) && (fabs(*step) < (maxstep/1.5))) {
			(*step) = 1.5*(*step);
		}

	}
	
	(*x) = x0 + 0.5*(r11 + 4*r41 + r51);
	(*y) = y0 + 0.5*(r12 + 4*r42 + r52);
	(*z) = z0 + 0.5*(r13 + 4*r43 + r53);
	
	kt17B(*x,*y,*z,bx,by,bz,Rsm,T1,T2);
	
	
	return;
}

void kt17Trace(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int LimType, int nParams, double *Params, //inputs
				int *nstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double *Rmsm, double *Rmso, double *FP) {//outputs

	/* initialize variables */
	int i;
	bool LimBits[32];
	GetBits(LimType,LimBits);
	*nstep = 1;
	for (i=0;i<maxlen;i++) {
		x[i] = NAN;
		y[i] = NAN;
		z[i] = NAN;
		bx[i] = NAN;
		by[i] = NAN;
		bz[i] = NAN;
		Rmsm[i] = NAN;
		Rmso[i] = NAN;
	}		
	for (i=0;i<20;i++) {
		// These are the footprints and field line lengths
		FP[i] = NAN;
	}
	double step = initstep;

	/* convert Parameters */
	double Rsm, T1, T2;
	if (nParams == 3) {
		/* here the parameters are Rsm, T1, T2 for the KT14 model*/
		Rsm = Params[0];
		T1 = Params[1];
		T2 = Params[2];
	} else {
		/* there should be just two parameters for the KT17 model: Rsun and DistIndex*/
		Rsm = (2.06873 - 0.00279*Params[1])*pow(Params[0],1.0/3.0);
		T1 = 6.495 + 0.0229*Params[1];
		T2 = 1.6245 + 0.0088*Params[1];			
	}


	x[0] = x0;
	y[0] = y0;
	z[0] = z0;
	kt17B(x[0],y[0],z[0],&bx[0],&by[0],&bz[0],Rsm,T1,T2);
	Rmsm[0] = sqrt(pow(x[0],2.0) + pow(y[0],2.0) + pow(z[0],2.0));
	Rmso[0] = sqrt(pow(x[0],2.0) + pow(y[0],2.0) + pow(z[0]+dz,2.0));
	
	/* Start by tracing along the direction of the field */
	bool InLimits = InsideTraceLims(x[*nstep-1],y[*nstep-1],z[*nstep-1],LimBits,1,Rsm);
	while (InLimits && *nstep < maxlen/2-2) {
		kt17StepRKM(x[*nstep-1],y[*nstep-1],z[*nstep-1],&step,maxstep,&x[*nstep],&y[*nstep],&z[*nstep],&bx[*nstep],&by[*nstep],&bz[*nstep],Rsm,T1,T2);
		InLimits = InsideTraceLims(x[*nstep],y[*nstep],z[*nstep],LimBits,1,Rsm);
/*		if (InLimits == false) {
			x[*nstep] = NAN;
			y[*nstep] = NAN;
			z[*nstep] = NAN;				
			bx[*nstep] = NAN;
			by[*nstep] = NAN;
			bz[*nstep] = NAN;
			break;
		}*/
		Rmsm[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep],2.0));
		Rmso[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep]+dz,2.0));
		(*nstep)++;		
	}
	
	/* reverse order results so far */
	ReverseElements(x,*nstep);
	ReverseElements(y,*nstep);
	ReverseElements(z,*nstep);
	ReverseElements(bx,*nstep);
	ReverseElements(by,*nstep);
	ReverseElements(bz,*nstep);
	ReverseElements(Rmsm,*nstep);
	ReverseElements(Rmso,*nstep);
	
	step=-step;

	/* Start tracing in opposite direction from original location */
	InLimits = InsideTraceLims(x[*nstep-1],y[*nstep-1],z[*nstep-1],LimBits,-1,Rsm);
	while (InLimits == true && *nstep < maxlen-2) {
		kt17StepRKM(x[*nstep-1],y[*nstep-1],z[*nstep-1],&step,maxstep,&x[*nstep],&y[*nstep],&z[*nstep],&bx[*nstep],&by[*nstep],&bz[*nstep],Rsm,T1,T2);
		InLimits = InsideTraceLims(x[*nstep],y[*nstep],z[*nstep],LimBits,-1,Rsm);
	/*	if (InLimits == false) {
			x[*nstep] = NAN;
			y[*nstep] = NAN;
			z[*nstep] = NAN;				
			bx[*nstep] = NAN;
			by[*nstep] = NAN;
			bz[*nstep] = NAN;
			break;
		}*/
		Rmsm[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep],2.0));
		Rmso[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep]+dz,2.0));
		(*nstep)++;		
	}

	/* First of all, determine which parts of the field lines are north and south of the equatorial plane, if any */
	double *Nflx,*Nfly,*Nflz,*NRmsm,*NRmso,*Sflx,*Sfly,*Sflz,*SRmsm,*SRmso;
	int nN, nS;
	
	NorthSouthFLs(x,y,z,Rmsm,Rmso,(*nstep),&Nflx,&Nfly,&Nflz,&NRmsm,&NRmso,&nN,&Sflx,&Sfly,&Sflz,&SRmsm,&SRmso,&nS);

	
	/* Find equatorial footprint */
	EqFootprint(Nflx, Nfly, Nflz, nN, Sflx, Sfly, Sflz, nS, &FP[16], &FP[17]);
	
	/* Find northern hemisphere footprint */
	if (nN > 1) {
		PlanetFootprints(Nflx, Nfly, Nflz, NRmsm, NRmso, nN, &FP[0], &FP[4], &FP[2], &FP[6], &FP[8], &FP[12], &FP[10], &FP[14]);
	}
		
	/* Find southern hemisphere footprint */
	if (nS > 1) {
		PlanetFootprints(Sflx, Sfly, Sflz, SRmsm, SRmso, nS,  &FP[1], &FP[5], &FP[3], &FP[7], &FP[9], &FP[13], &FP[11], &FP[15]);
	}	
	
	/* Find field line length */
	if (!isnan(FP[2]) && !isnan(FP[3])) {
		FieldLineLength(x,y,z,*nstep,Rmso,&FP[18],&FP[19]);
	}
	
	if (nN > 0) {
		free(Nflx);
		free(Nfly);
		free(Nflz);
		free(NRmso);
		free(NRmsm);
	}
	if (nS > 0) {
		free(Sflx);
		free(Sfly);
		free(Sflz);
		free(SRmso);
		free(SRmsm);
	}
	return;

}

//void kt17TraceOld(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				//double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				//double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double Rsm, double T1, double T2){


	///* initialise some variables */
	//*nstep = 1;
	//*mlatn = NAN;
	//*mlats = NAN;
	//*latn = NAN;
	//*lats = NAN;
	//*mltn = NAN;
	//*mlts = NAN;
	//*lctn = NAN;
	//*lcts = NAN;
	//*lshell = NAN;
	//*mlte = NAN;
	//*fl_len = NAN;
	//int i;
	//for (i=1;i<maxlen;i++) {
		//x[i] = NAN;
		//y[i] = NAN;
		//z[i] = NAN;
		//bx[i] = NAN;
		//by[i] = NAN;
		//bz[i] = NAN;
		//Rmsm[i] = NAN;
		//Rmso[i] = NAN;
	//}
	//double step = initstep;
	
	//x[0] = x0;
	//y[0] = y0;
	//z[0] = z0;
	//kt17B(x[0],y[0],z[0],&bx[0],&by[0],&bz[0],Rsm,T1,T2);
	//Rmsm[0] = sqrt(pow(x[0],2.0) + pow(y[0],2.0) + pow(z[0],2.0));
	//Rmso[0] = sqrt(pow(x[0],2.0) + pow(y[0],2.0) + pow(z[0]+dz,2.0));
	
	///* Start by tracing along the direction of the field */
	//bool InLimits = InsideTraceLims(x[*nstep-1],y[*nstep-1],z[*nstep-1],LimType,1,Rsm);
	//while (InLimits && *nstep < maxlen/2-2) {
		//kt17StepRKM(x[*nstep-1],y[*nstep-1],z[*nstep-1],&step,maxstep,&x[*nstep],&y[*nstep],&z[*nstep],&bx[*nstep],&by[*nstep],&bz[*nstep],Rsm,T1,T2);
		//InLimits = InsideTraceLims(x[*nstep],y[*nstep],z[*nstep],LimType,1,Rsm);
		//if (InLimits == false) {
			//x[*nstep] = NAN;
			//y[*nstep] = NAN;
			//z[*nstep] = NAN;				
			//bx[*nstep] = NAN;
			//by[*nstep] = NAN;
			//bz[*nstep] = NAN;
			//break;
		//}
		//Rmsm[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep],2.0));
		//Rmso[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep]+dz,2.0));
		//(*nstep)++;		
	//}

	///* reverse order results so far */
	//ReverseElements(x,*nstep);
	//ReverseElements(y,*nstep);
	//ReverseElements(z,*nstep);
	//ReverseElements(bx,*nstep);
	//ReverseElements(by,*nstep);
	//ReverseElements(bz,*nstep);
	//ReverseElements(Rmsm,*nstep);
	//ReverseElements(Rmso,*nstep);
	
	//step=-step;

	///* Start tracing in opposite direction from original location */
	//InLimits = InsideTraceLims(x[*nstep-1],y[*nstep-1],z[*nstep-1],LimType,-1,Rsm);
	//while (InLimits == true && *nstep < maxlen-2) {
		//kt17StepRKM(x[*nstep-1],y[*nstep-1],z[*nstep-1],&step,maxstep,&x[*nstep],&y[*nstep],&z[*nstep],&bx[*nstep],&by[*nstep],&bz[*nstep],Rsm,T1,T2);
		//InLimits = InsideTraceLims(x[*nstep],y[*nstep],z[*nstep],LimType,-1,Rsm);
		//if (InLimits == false) {
			//x[*nstep] = NAN;
			//y[*nstep] = NAN;
			//z[*nstep] = NAN;				
			//bx[*nstep] = NAN;
			//by[*nstep] = NAN;
			//bz[*nstep] = NAN;
			//break;
		//}
		//Rmsm[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep],2.0));
		//Rmso[*nstep] = sqrt(pow(x[*nstep],2.0) + pow(y[*nstep],2.0) + pow(z[*nstep]+dz,2.0));
		//(*nstep)++;		
	//}
	
	///* more stuff will be added here to interpolate and find footprints */
	
	///* First of all, determine which parts of the field lines are north and south of the equatorial plane, if any */
	//double *Nflx,*Nfly,*Nflz,*NRmsm,*NRmso,*Sflx,*Sfly,*Sflz,*SRmsm,*SRmso;
	//int nN, nS;
	
	//NorthSouthFLs(x,y,z,Rmsm,Rmso,(*nstep),&Nflx,&Nfly,&Nflz,&NRmsm,&NRmso,&nN,&Sflx,&Sfly,&Sflz,&SRmsm,&SRmso,&nS);

	
	///* Find equatorial footprint */
	//EqFootprint(Nflx, Nfly, Nflz, nN, Sflx, Sfly, Sflz, nS, lshell, mlte);
	
	///* Find northern hemisphere footprint */
	//if (nN > 1) {
		//PlanetFootprints(Nflx, Nfly, Nflz, NRmsm, NRmso, nN, mlatn, mltn, latn, lctn);
	//}
		
	///* Find southern hemisphere footprint */
	//if (nS > 1) {
		//PlanetFootprints(Sflx, Sfly, Sflz, SRmsm, SRmso, nS, mlats, mlts, lats, lcts);
	//}	
	
	///* Find field line length */
	//if (!isnan(*latn) && !isnan(*lats)) {
		//*fl_len = FieldLineLength(x,y,z,*nstep,Rmso);
	//}
	
	//if (nN > 0) {
		//free(Nflx);
		//free(Nfly);
		//free(Nflz);
		//free(NRmso);
		//free(NRmsm);
	//}
	//if (nS > 0) {
		//free(Sflx);
		//free(Sfly);
		//free(Sflz);
		//free(SRmso);
		//free(SRmsm);
	//}
	//return;
//}
/*
void kt17TraceScaled(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double Rsun, double DistIndex) {
	
	double f, Rsm, T1, T2;
	f = 2.06873 - 0.00279*DistIndex;
	Rsm = f*pow(Rsun,1.0/3.0);
	T1 = 6.495 + 0.0229*DistIndex;
	T2 = 1.6245 + 0.0088*DistIndex;	
					
	kt17Trace(x0,y0,z0,maxlen,initstep,maxstep,nstep,x,y,z,bx,by,bz,Rmsm,Rmso,mlatn,mlats,latn,lats,mltn,mlts,lctn,lcts,lshell,mlte,fl_len,LimType,Rsm,T1,T2);				
}
*/
void kt17MultiTrace(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int LimType, int nParams, double *Params, //inputs
				int *nstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double *Rmsm, double *Rmso, double *FP) {
	
	int i;
	double *Pptr;
	int nP;
	//This little extra should give us the option to use one set of parameters for all traces
	// or to use different parameters for each trace
	if ((nParams == 2) || (nParams == 3)) {
		nP = nParams;
	} else {
		nP = nParams/n;
	}
	for (i=0;i<n;i++) {
		if ((nParams == 2) || (nParams == 3)) {
			Pptr = &Params[0];
		} else {
			Pptr = &Params[i*nP];
		}
		kt17Trace(x0[i],y0[i],z0[i],maxlen,initstep,maxstep,LimType,nP,Pptr,
				&nstep[i],&x[i*maxlen],&y[i*maxlen],&z[i*maxlen],&bx[i*maxlen],&by[i*maxlen],&bz[i*maxlen],&Rmsm[i*maxlen],&Rmso[i*maxlen],&FP[i*20]); 
	}				
}

/*
void kt17MultiTrace(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double *Rsm, double *T1, double *T2){
	int i;
	for (i=0;i<n;i++) {
		kt17Trace(x0[i],y0[i],z0[i],maxlen,initstep,maxstep,&nstep[i],&x[i*maxlen],&y[i*maxlen],&z[i*maxlen],&bx[i*maxlen],&by[i*maxlen],&bz[i*maxlen],&Rmsm[i*maxlen],
			&Rmso[i*maxlen],&mlatn[i],&mlats[i],&latn[i],&lats[i],&mltn[i],&mlts[i],&lctn[i],&lcts[i],&lshell[i],&mlte[i],&fl_len[i],LimType,Rsm[i],T1[i],T2[i]);
	}			
					
	return;
}
*/
/*
void kt17MultiTraceScaled(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double *Rsun, double *DistIndex){		
	int i;
	for (i=0;i<n;i++) {
		kt17TraceScaled(x0[i],y0[i],z0[i],maxlen,initstep,maxstep,&nstep[i],&x[i*maxlen],&y[i*maxlen],&z[i*maxlen],&bx[i*maxlen],&by[i*maxlen],&bz[i*maxlen],&Rmsm[i*maxlen],
			&Rmso[i*maxlen],&mlatn[i],&mlats[i],&latn[i],&lats[i],&mltn[i],&mlts[i],&lctn[i],&lcts[i],&lshell[i],&mlte[i],&fl_len[i],LimType,Rsun[i],DistIndex[i]);
	}			
					
	return;
}
*/
