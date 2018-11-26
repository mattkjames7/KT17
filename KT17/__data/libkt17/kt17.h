#ifndef __kt17_h__
#define __kt17_h__
#include <stdio.h>
#include <math.h>
using namespace std;

extern "C" {
	bool WithinMP(double x, double y, double z, double Rsm);
	void kt17DipoleField(double x, double y, double z, double *Bx, double *By, double *Bz);
	void kt17DipoleShield(double x, double y, double z, double *Bx, double *By, double *Bz);
	void kt17DipoleB(double x, double y, double z, double *Bx, double *By, double *Bz);
	double kt17DiskThickness(double x, double y);
	void kt17DiskField(double xin, double yin, double zin, double *Bx, double *By, double *Bz);
	void kt17DiskShield(double x, double y, double z, double *Bx, double *By, double *Bz);
	void kt17DiskB(double x, double y, double z, double *Bx, double *By, double *Bz, double T1);
	double kt17QuasiHarrisThickness(double x, double y);
	void kt17QuasiHarrisField(double x, double y, double z, double *Bx, double *By, double *Bz);
	void kt17QuasiHarrisShield(double x, double y, double z, double *Bx, double *By, double *Bz);
	void kt17QuasiHarrisB(double x, double y, double z, double *Bx, double *By, double *Bz, double T2);
	void kt17B(double xin, double yin, double zin, double *Bx, double *By, double *Bz, double Rsm, double T1, double T2);
	void kt17Barray(int n, double *xin, double *yin, double *zin, double *Bx, double *By, double *Bz, int nP, int Plen, double *Params);
	void kt17(double x, double y, double z, double *Bx, double *By, double *Bz, double Rsun, double DistIndex);
	//void kt17array(int n, double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Rsun, double *DistIndex);
}

#endif
