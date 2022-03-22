#include "libkt.h"


/***********************************************************************
 * NAME : void ModelField(n,x,y,z,lP,nP,Params,Bx,By,Bz)
 * 
 * DESCRIPTION : Calculate the KT14/17 magnetic field.
 * 
 * INPUTS : 
 * 		int n			Number field vectors to find.
 * 		double *x		Position array (Rm).
 * 		double *y		Position array (Rm).
 * 		double *z		Position array (Rm).
 * 		int lP			Number or sets of parameters.
 * 		int nP 			Number of parameters (2 or 3)
 * 		double **Params	2D array of parameters, shape (lP,nP)
 * 
 * OUTPUTS :
 * 		double *Bx		Output magnetic field (nT).
 * 		double *By		Output magnetic field (nT).
 * 		double *Bz		Output magnetic field (nT).
 * 	  
 * ********************************************************************/
void ModelField(int n, double *x, double *y, double *z,
				int lP, int nP, double **Params,
				double *Bx, double *By, double *Bz) {
	
	/* create a mdoel object */
	KTModel kt;	
	
	int i, ip;
	bool inMP;
	if (nP == 2) {
		/* this would imply that we are using the KT17 parameters */
		for (i=0;i<n;i++) {
			ip = i % lP;
			kt.SetParams(Params[ip][0],Params[ip][1]);
			/* check that we are within the MP */
			inMP = WithinMP(x[i],y[i],z[i],kt.Rsm_);
			if (inMP) {
				kt.Field(x[i],y[i],z[i],&Bx[i],&By[i],&Bz[i]);
			} else {
				Bx[i] = NAN;
				By[i] = NAN;
				Bz[i] = NAN;
			}
		}
	} else {
		/* nP should equal 3 for the KT14 parameters */
		for (i=0;i<n;i++) {
			ip = i % lP;
			kt.SetParams(Params[ip][0],Params[ip][1],Params[ip][2]);
			/* check that we are within the MP */
			inMP = WithinMP(x[i],y[i],z[i],kt.Rsm_);
			if (inMP) {
				kt.Field(x[i],y[i],z[i],&Bx[i],&By[i],&Bz[i]);
			} else {
				Bx[i] = NAN;
				By[i] = NAN;
				Bz[i] = NAN;
			}
		}				
	}
	
}




/***********************************************************************
 * NAME : void TraceField(n,x0,y0,z0,nP,P0,P1,P2,BoundMP,BoundTail,
 * 						BoundSurface,MaxLen,MaxStep,InitStep,
 * 						MinStep,ErrMax,Delta,Verbose,TraceDir,nstep,
 * 						x,y,z,Bx,By,Bz,Rmsm,Rmso,S,Rnorm,FP,nalpha,
 * 						alpha,halpha)
 * 
 * DESCRIPTION : Trace field lines using the KT14/17 model.
 * 
 * INPUTS : 
 * 		int n		number of traces.
 * 		double *x	Trace start position MSM (Rm)
 * 		double *y	Trace start position MSM (Rm)
 * 		double *z	Trace start position MSM (Rm)
 * 		int nP		Number of parameters 2(KT17) or 3(KT14).
 * 		double *P0	Rsun in AU (KT17) or Rsm in Rm (KT14)
 * 		double *P1	Disturbance index from 0-97 (KT17) or t1 disk 
 * 					magnitude (KT14).
 * 		double *P2	Unused (KT17) or t2 quasi harris sheet magnitude 
 * 					(KT14).
 * 		bool BoundMP	If true, then the trace will stop at the 
 * 					magnetopause.
 * 		double BoundTail	X limit to stop trace in magnetotail,
 * 					ignored if positive.
 * 		int	BoundSurface	Surface on which to end the field trace:
 * 					1 - Stop at the planetary surface
 * 					2 - Stop at the planetary core (0.832 Rm).
 * 					3 - Stop at dipole at 1 Rm
 * 					4 - Stop at dipole at 0.832 Rm (core radius)
 * 					5 - Stop at northern surface, southern dipole at 
 * 						1 Rm (virtual surface).
 * 					6 - Stop at northern core and southern dipole at
 * 						0.832 Rm.
 * 		int MaxLen		Maximum length of each trace.
 * 		double MaxStep	Maximum step size (Rm)
 * 		double InitStep Initial step size (Rm)
 * 		double MinStep	Minimum step size (Rm)
 * 		double ErrMax	Maximum error.
 * 		double Delta	Distance between adjacent field lines (Rm) to be
 * 						used int he calculation of h_alpha.
 * 		bool Verbose	Display tracing progress.
 * 		int TraceDir	Direction in which to trace:
 * 						0 - trace in both directions
 * 						1 - along the field (towards north)
 * 						-1 - opposite to the field (towards south)		
 * 		int nalpha		The number of alpha values.
 * 		double *alpha	Array of alphas (degrees) 0 = toroidal.
 * 
 * OUTPUTS :
 * 		int *nstep		Number of trace steps
 * 		double **x		Trace positions
 * 		double **y		Trace positions
 * 		double **z		Trace positions
 * 		double **Bx		Trace field
 * 		double **By		Trace field
 * 		double **Bz		Trace field
 * 		double **Rmsm	Trace radial coordinates (MSM)
 * 		double **Rmso	Trace radial coordinates (MSO)
 * 		double **S		Trace distance.
 *		double **Rnorm	Array of Rnorms.
 * 		double **FP		Output footprint coords, shape (n,18), where n
 * 						is the number of traces and the elements in the
 * 						2nd dimension correspond to the following 
 * 						footprints:
 * 						0:	North planetary latitude
 * 						1:	North planetary local time
 * 						2:	South planetary latitude
 * 						3:	South planetary local time
 * 						4:	North core latitude
 * 						5:	North core local time
 * 						6:	South core latitude
 * 						7:	South core local time
 * 						8:	North dipole R=1 latitude
 * 						9:	North dipole R=1 local time
 * 						10:	South dipole R=1 latitude
 * 						11:	South dipole R=1 local time
 * 						12:	North dipole R=0.832 latitude
 * 						13:	North dipole R=0.832 local time
 * 						14:	South dipole R=0.832 latitude
 * 						15:	South dipole R=0.832 local time
 * 						16: L-shell
 * 						17:	Equatorial footprint magnetic local time.
 * 						
 * 						NOTE: core is assumed to be a sphere at 
 * 						Rmso = 0.832 Rm and the dipole footprints
 * 						are footprints on a sphere centered on the 
 * 						dipole rather than the planet itself.
 * 		double *halpha	1-D array containing all trace h_alphas,
 * 						reshape in Python to (n,nalpha,MaxLen)
 * 	  
 * ********************************************************************/
void TraceField(int n, double *x0, double *y0, double *z0,
				int nP, double *P0, double *P1, double *P2,
				bool BoundMP, double BoundTail, int BoundSurface,
				int MaxLen, double MaxStep, double InitStep,
				double MinStep, double ErrMax, double Delta,
				bool Verbose, int TraceDir,
				int *nstep,
				double **x, double **y, double **z,
				double **Bx, double **By, double **Bz,
				double **Rmsm, double **Rmso,
				double **S, double **Rnorm, double **FP,
				int nalpha, double *alpha, double *halpha) {

	/* create the trace object */
	Trace T;
	
	/* input the start positions */
	T.InputPos(n,x0,y0,z0);
	
	/* input the parameters for the model */
	if (nP == 2) {
		/* KT17 parameters */
		T.SetModelParams(P0,P1);
	} else {
		/* KT14 Parameters */
		T.SetModelParams(P0,P1,P2);
	}
	
	/* configure model boundary conditions */
	T.SetTraceBounds(BoundMP,BoundTail,BoundSurface);
	
	/* configure trace */
	T.SetTraceCFG(MaxLen,MaxStep,InitStep,MinStep,ErrMax,Delta,Verbose,TraceDir);
	
	/* run the trace */
	T.TraceField(nstep,x,y,z,Rmsm,Rmso,Bx,By,Bz);
	
	/* calculate the trace distance */
	T.CalculateTraceDist(S);
	
	/* calculate the footprints */
	T.CalculateTraceFP(FP);
	
	/* calculate the Rnorm */
	T.CalculateTraceRnorm(Rnorm);
	
	/* alpha if appropriate */
	if (nalpha > 0) {
		T.SetAlpha(nalpha,alpha);
		
		T.CalculateHalpha(halpha);
	}
	
}
