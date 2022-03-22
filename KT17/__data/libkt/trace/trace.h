#ifndef __TRACE_H__
#define __TRACE_H__
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../model/kt17.h"
#include "../model/ktmodel.h"
#include "../magnetopause/withinmp.h"
#include "interptraceclosestpos.h"
#include "conttrace.h"
#include "latlt.h"
#include "fieldlinernorm.h"

/***********************************************************************
 * This object will store a bunch of field traces within it.
 * 
 * It will have the ability to either allocate and store field vectors
 * and positions, or to accept pointers which can be created externally
 * (e.g. inside Python)
 * 
 * There will be optional member functions which obtain things like 
 * footprints and h_alphas.
 * 
 * The basic trace will be in GSM/GSW.
 * 
 * Other coordinate systems will be calculated as needed.
 * 
 * ********************************************************************/
class Trace {
	
	public:
		/* initialize the object */
		Trace();
		
		/* delete the object */
		~Trace();
		
		/* copy constructor */
	//	Trace(const Trace &);
		
		/* this will take in the input positions where the traces start*/
		void InputPos(int,double*,double*,double*);
		
		/* set model parameters */
		void SetModelParams(double*, double*, double *);
		void SetModelParams(double*, double*);
		
		/* set trace boundary conditions */
		void SetTraceBounds(bool,double,int);
		void SetTraceBounds();
		
		/* set the trace configuration */
		void SetTraceCFG(int, double,double,double,double,double,bool,int);
		void SetTraceCFG();
		
		/* polarization stuff */
		void SetAlpha(int,double*);

			
		/* trace functions */
		void TraceField(int*,double**,double**,double**,
						double**,double**,
						double**,double**,double**);
		void TraceField(int*,double**,double**,double**,
						double**,double**,double**);
		void TraceField(int*);
		void TraceField();
		void StepVector(double,double,double,double,double*,double*,double*);
		bool ContinueTrace(double,double,double,double*,double*,double*);
		void Step(double,double,double,double*,double*,double*,double*,
					double*,double*,double*,double*,double*,double*);
		void ReverseElements(int, double*);
		void RKMTrace(	double,double,double,int*,
						double*,double*,double*,
						double*,double*,double*,
						double*,double*,double*);		
	
		/* calculate trace distance,R,Rnorm */
		void CalculateTraceDist(double**);
		void CalculateTraceDist();
		void _CalculateTraceDist();
		void CalculateTraceRnorm(double**);
		void CalculateTraceRnorm();
		void _CalculateTraceRnorm();
	
		/* Calculate footprints */
		void CalculateTraceFP(double**);
		void CalculateTraceFP();
		void _CalculateTraceFP();
		
		/* calculate halpha */
		void CalculateHalpha();
		void CalculateHalpha(double*);
		void CalculateHalpha(double***);
		void CalculateHalpha(double*,double***);
	
		/* return things*/
		void GetTraceNstep(int*);
		void GetTrace(double**,double**,double**);
		void GetTrace(double**,double**,double**,double**,double**,double**);
		void GetTraceDist(double**);
		void GetTraceRmsm(double**);
		void GetTraceRmso(double**);
		void GetTraceRnorm(double**);
		void GetTraceFootprints(double**);
		void GetTraceHalpha(double*);	/* python will use this */
		void GetTraceHalpha(double***); /* no idea how to link this to python*/
		
		Trace TracePosition(int,double,double,double);
	

		/* input coords */
		int n_;
		bool *inMP_;
		double *x0_, *y0_, *z0_;  

		/* trace params */
		int MaxLen_;
		double MaxStep_, MinStep_, InitStep_;
		bool Verbose_;
		int TraceDir_;
		double ErrMax_;
		double MaxR_;
		
		/* model params */
		double *Rsm_, *t1_, *t2_;

		/* trace coords */
		int *nstep_;
		double **x_, **y_, **z_, **zmso_;
	
		/* trace fields */
		double **bx_, **by_, **bz_;

		/* trace end points */
		/* lots of variables here:
		 * x,y,z obviously refers to the coordinate
		 * n/s - north/south footprint
		 * v - on a virtual sphere
		 * c - core 
		 * e - magnetic equator */
		double *xfn_, *yfn_, *zfn_;
		double *xfs_, *yfs_, *zfs_;
		double *xfnc_, *yfnc_, *zfnc_;
		double *xfsc_, *yfsc_, *zfsc_;
		double *xfnv_, *yfnv_, *zfnv_;
		double *xfsv_, *yfsv_, *zfsv_;
		double *xfnvc_, *yfnvc_, *zfnvc_;
		double *xfsvc_, *yfsvc_, *zfsvc_;
		double *xfe_, *yfe_, *zfe_;

	private:
		/* booleans to tell the object what has been allocated */
		bool allocInputPos_;
		bool allocModelParams_;
		bool allocAlpha_;
		bool allocHalpha_;
		bool allocHalpha3D_;
		bool allocTrace_;
		bool allocNstep_;
		bool allocMP_;
		bool allocDist_;
		bool allocRnorm_;
		bool allocR_;
		bool allocFootprints_;
		bool allocEqFP_;
		bool allocEndpoints_;
		bool allocZmso_;
			
		
		/* booleans which tell the object what has been done */
		bool setModelParams_;
		bool setBounds_;
		bool setTraceConfig_;
		bool setAlpha_;
		bool setHalpha_;
		bool setTrace_;
		bool setDist_;
		bool setRnorm_;
		bool setFootprints_;
	
		/* boundary stuff */
		bool BoundMP_;
		double BoundTail_;
		int BoundSurface_;
		CTFuncMP ctfMP_;
		CTFuncTX ctfTX_;
		CTFuncSC ctfSC_;
	
		/* field length, R, Rnorm, Halpha, Footprints */
		int nalpha_;
		double *alpha0_, *alpha1_;
		double Delta_;
		double **S_;
		double **Rmsm_, **Rmso_;
		double **Rnorm_;
		double *Halpha_;
		double ***Halpha3D_;
		double **FP_;
		
		/* model */
		KTModel ktmodel;
		
		/* FP Functions */
		void _FPCoords();
		void _InterpPos(double*, double*, double*,
						double*, double*, double,
						double*, double*, double*);
		void _SingleTraceFP(int);
		
	
		/* hidden trace functions */
		void _TraceField();
		void _AllocTrace();
		void _AllocTraceR();
		void _AllocZmso();
		
		/* halpha functions */
		bool _CheckHalpha();
		void _CalculateHalpha();
		void _CalculateTraceHalpha(int,int,double*);
		void _CalculateHalphaStartPoints(int i, int j,
							double *xe0, double *ye0, double *ze0,
							double *xe1, double *ye1, double *ze1);
};



#endif

