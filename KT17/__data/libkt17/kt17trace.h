#ifndef __kt17trace_h__
#define __kt17trace_h__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kt17.h"
using namespace std;

extern "C" {
	void UnitVector(double ix, double iy, double iz, double *ox, double *oy, double *oz);
	void ReverseElements(double x[], int n);
	bool InsideTraceLims(double x, double y, double z, bool* bits, int dir, double Rsm);
	void NorthSouthFLs(double flx[],double fly[],double flz[], double Rmsm[], double Rmso[], int N, double **Nflx, double **Nfly, double **Nflz, double **NRmsm, double **NRmso, int *nn, double **Sflx, double **Sfly, double **Sflz, double **SRmsm, double **SRmso, int *ns);
	double linterp(double x0, double x1, double y0, double y1, const double xt);
	void EqFootprint(double *Nflx, double *Nfly, double *Nflz, int nN, double *Sflx, double *Sfly, double *Sflz, int nS, double *lshell, double *mlte);
	double min(double r, double *a, int n);
	void PlanetFootprints(double *flx, double *fly, double *flz, double *Rmsm, double *Rmso, int n, double *mlat, double *mlt, double *lat, double *lct, double *mlatc, double *mltc, double *latc, double *lctc);
	void FieldLineLength(double *x, double *y, double *z, int n, double *r, double *len, double *lenc);
	void Rvecs(double x0, double y0, double z0, double *rx, double *ry, double *rz, double Rsm, double T1, double T2, double step3);
	void kt17StepRKM(double x0, double y0, double z0, double *step, double maxstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double Rsm, double T1, double T2);
/*	void kt17Trace(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double Rsm, double T1, double T2);
	void kt17TraceScaled(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double Rsun, double DistIndex);
	void kt17MultiTrace(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double *Rsm, double *T1, double *T2);
	void kt17MultiTraceScaled(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int *nstep, double x[], double y[], double z[],
				double bx[], double by[], double bz[], double Rmsm[], double Rmso[], double *mlatn, double *mlats, double *latn, double *lats, double *mltn, double *mlts,
				double *lctn, double *lcts, double *lshell, double *mlte, double *fl_len, int LimType, double *Rsun, double *DistIndex);	*/
	void kt17Trace(double x0, double y0, double z0, int maxlen, double initstep, double maxstep, int LimType, int nParams, double *Params, int *nstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double *Rmsm, double *Rmso, double *FP);
	void kt17MultiTrace(double *x0, double *y0, double *z0, int n, int maxlen, double initstep, double maxstep, int LimType, int nParams, double *Params, int *nstep, double *x, double *y, double *z, double *bx, double *by, double *bz, double *Rmsm, double *Rmso, double *FP);								
}

const double dz = 0.196;
#define ErrMax 0.0001
#define minstep 0.001
#endif
