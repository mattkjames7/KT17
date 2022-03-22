#include "trace.h"

/***********************************************************************
 * NAME : Trace()
 * 
 * DESCRIPTION : Constructor for the Trace object.
 * 
 * ********************************************************************/
Trace::Trace() {
	/* initialize all of the boolean parameters */
	allocInputPos_ = false;
	setModelParams_ = false;
	allocModelParams_ = false;
	setBounds_ = false;
	setTraceConfig_ = false;
	allocAlpha_ = false;
	setAlpha_ = false;
	allocHalpha_ = false;
	allocHalpha3D_ = false;
	setHalpha_ = false;
	allocTrace_ = false;
	setTrace_ = false;
	allocNstep_ = false;
	allocMP_ = false;
	allocDist_ = false;
	setDist_ = false;
	allocRnorm_ = false;	
	setRnorm_ = false;	
	allocFootprints_ = false;
	allocEqFP_ = false;
	allocEndpoints_ = false;
	setFootprints_ = false;
	allocR_ = false;
	allocZmso_ = false;

	
	/* default trace parameters */
	SetTraceCFG();
	
}

/***********************************************************************
 * NAME : ~Trace()
 * 
 * DESCRIPTION : Destructor for the Trace object.
 * 
 * ********************************************************************/
Trace::~Trace() {

	/* check for each allocated variable and delete it*/
	int i, j, k = 0;
	
	/* starting positions */
	if (allocInputPos_) {
		delete[] x0_;
		delete[] y0_;
		delete[] z0_;
	}

	/* parameters */
	if (allocModelParams_) {
		delete[] Rsm_;
		delete[] t1_;
		delete[] t2_;
	}

	/* traces */
	if (allocMP_) {
		delete[] inMP_;
	}
	if (allocTrace_) {
		for (i=0;i<n_;i++) {
			delete[] x_[i];
			delete[] y_[i];
			delete[] z_[i];
			delete[] bx_[i];
			delete[] by_[i];
			delete[] bz_[i];
		}
		delete[] x_;
		delete[] y_;
		delete[] z_;
		delete[] bx_;
		delete[] by_;
		delete[] bz_;
	}
	if (allocZmso_) {
		for (i=0;i<n_;i++) {
			delete[] zmso_[i];
		}
		delete[] zmso_;
	}
	
	/* nstep */
	if (allocNstep_) {
		delete[] nstep_;
	}

	/* field line footprints */
	if (allocFootprints_) {
		for (i=0;i<n_;i++) {
			delete[] FP_[i];
		}
		delete[] FP_;
	}

	/* field line distance */
	if (allocDist_) {
		for (i=0;i<n_;i++) {
			delete[] S_[i];
		}
		delete[] S_;
	}

	/* radial distance */
	if (allocR_) {
		for (i=0;i<n_;i++) {
			delete[] Rmsm_[i];
			delete[] Rmso_[i];
		}
		delete[] Rmsm_;
		delete[] Rmso_;
	}

	/* r norm distance */
	if (allocRnorm_) {
		for (i=0;i<n_;i++) {
			delete[] Rnorm_[i];
		}
		delete[] Rnorm_;
	}

	/* h alpha*/
	if (allocAlpha_) {
		delete[] alpha0_;
		delete[] alpha1_;
	}
	if (allocHalpha_) {
		delete[] Halpha_;
	}
	if (allocHalpha3D_) {
		for (i=0;i<n_;i++) {
			for (j=0;j<nalpha_;j++) {
				delete[] Halpha3D_[i][j];
			}
			delete[] Halpha3D_[i];
		}
		delete[] Halpha3D_;
	}

	/* footprint/endpoints */
	if (allocEndpoints_) {
		delete[] xfn_;
		delete[] yfn_;
		delete[] zfn_;
		delete[] xfs_;
		delete[] yfs_;
		delete[] zfs_;
		delete[] xfnc_;
		delete[] yfnc_;
		delete[] zfnc_;
		delete[] xfsc_;
		delete[] yfsc_;
		delete[] zfsc_;
		delete[] xfnv_;
		delete[] yfnv_;
		delete[] zfnv_;
		delete[] xfsv_;
		delete[] yfsv_;
		delete[] zfsv_;
		delete[] xfnvc_;
		delete[] yfnvc_;
		delete[] zfnvc_;
		delete[] xfsvc_;
		delete[] yfsvc_;
		delete[] zfsvc_;
	}	
	if (allocEqFP_) {
		delete[] xfe_;
		delete[] yfe_;
		delete[] zfe_;
	}
}

/***********************************************************************
 * NAME : void Trace::InputPos(n,x,y,z)
 * 
 * DESCRIPTION : Input the starting positions of each trace.
 * 
 * INPUTS : 
 * 		int n		Number of traces.
 * 		double x	x-MSM coordinate of starting points.
 * 		double y	y-MSM coordinate of starting points.
 * 		double z	z-MSM coordinate of starting points.
 *	
 * ********************************************************************/
void Trace::InputPos(	int n, double *x, double *y, double *z) {

	/* check that we have not already done this */
	if (allocInputPos_) {
		printf("Input positions already set, ignoring...\n");
		return;
	}

	/* allocate the memory to store the input coords */
	n_ = n;
	x0_ = new double[n];
	y0_ = new double[n];
	z0_ = new double[n];
	
	
	/* convert the coordinates to GSM */
	int i;	
	for (i=0;i<n_;i++) {
		x0_[i] = x[i];
		y0_[i] = y[i];
		z0_[i] = z[i];
	}	
	
	/* set the flag so we know to delete  again once the object is deleted */
	allocInputPos_ = true;
						
}						

/***********************************************************************
 * NAME : void Trace::SetModelParams(Rsm,t1,t2)
 * 
 * DESCRIPTION : Set the model parameters (KT14) for each trace.
 * 
 * INPUTS : 
 * 		double *Rsm		Array of subsolar standof distances (Rm).
 * 		double *t1		Tail disk magnitudes.
 * 		double *t2		Quasi-harris sheet magnitudes.
 *		
 * ********************************************************************/
void Trace::SetModelParams(double *Rsm, double *t1, double *t2) {
	
	/* set the pointers to parameters provided externally */
	Rsm_ = Rsm;
	t1_ = t1;
	t2_ = t2;

	/* set the flags to say that we have the parameters, but they don't
	* need to be deleted */
	setModelParams_ = true;
	allocModelParams_ = false;
}

/***********************************************************************
 * NAME : void Trace::SetModelParams(Rsun,DistIndex)
 * 
 * DESCRIPTION : Set model parameters (KT17) for each trace.
 * 
 * INPUTS : 
 * 		double *Rsun		Distance from the Sun (AU).
 * 		double *DistIndex	Anderson et al 2013 disturbance index (0-97).
 *	
 * ********************************************************************/
void Trace::SetModelParams(double *Rsun, double *DistIndex) {
	
	/* allocate parameters*/
	Rsm_ = new double[n_];
	t1_ = new double[n_];
	t2_ = new double[n_];
	
	/* convert from KT17 to KT14 */
	int i;
	for (i=0;i<n_;i++) {
		KT14Params(Rsun[i],DistIndex[i],&Rsm_[i],&t1_[i],&t2_[i]);
	}

	/* set the flags to say that we have the parameters, and that they
	* need to be deleted */
	setModelParams_ = true;
	allocModelParams_ = true;
}

/***********************************************************************
 * NAME : Trace::SetTraceBounds(MP,TailX,EndSurface)
 * 
 * DESCRIPTION : Tell trace routine where to stop tracing.
 * 
 * INPUTS : 
 * 		bool MP		If true, then the trace will stop if the MP is
 * 					encountered.
 *		double TailX	x-coordinate at which to stop tracing in the tail,
 * 					positive values will be ignored and tracing will
 * 					continue until MaxLen_.
 *		int	EndSurface	Surface on which to end the field trace:
 * 					1 - Stop at the planetary surface
 * 					2 - Stop at the planetary core (0.832 Rm).
 * 					3 - Stop at dipole at 1 Rm
 * 					4 - Stop at dipole at 0.832 Rm (core radius)
 * 					5 - Stop at northern surface, southern dipole at 
 * 						1 Rm (virtual surface).
 * 					6 - Stop at northern core and southern dipole at
 * 						0.832 Rm.
 *
 * 
 * ********************************************************************/
void Trace::SetTraceBounds(bool MP, double TailX, int EndSurface) {

	BoundMP_ = true;
	BoundTail_ = TailX;
	BoundSurface_ = EndSurface;
	
	/* set the function pointers */
	if (MP) {
		ctfMP_ = &ContTraceMP;
	} else { 
		ctfMP_ = &ContTraceMPDummy;
	}
	if (TailX < 0) {
		ctfTX_ = &ContTraceTail;
	} else { 
		ctfTX_ = &ContTraceTailDummy;
	}
	if (EndSurface == 1) {
		ctfSC_ = &ContTraceSurface1;
	} else if (EndSurface == 2) {
		ctfSC_ = &ContTraceSurface2;
	} else if (EndSurface == 3) {
		ctfSC_ = &ContTraceSurface3;
	} else if (EndSurface == 4) {
		ctfSC_ = &ContTraceSurface4;
	} else if (EndSurface == 5) {
		ctfSC_ = &ContTraceSurface5;
	} else {
		ctfSC_ = &ContTraceSurface6;
	}
	setBounds_ = true;
}

/***********************************************************************
 * NAME : Trace::SetTraceBounds()
 * 
 * DESCRIPTION : Tell trace routine where to stop tracing. (defaults)
 *
 * ********************************************************************/
void Trace::SetTraceBounds() {

	SetTraceBounds(true,-10.0,6);
}


/***********************************************************************
 * NAME : void Trace::SetTraceConfig(MaxLen,MaxStep,InitStep,MinStep,
 * 									ErrMax,Delta,Verbose,TraceDir)
 * 
 * DESCRIPTION : Sets a bunch of parameters for the trace.
 * 
 * INPUTS : 
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
 *	
 * ********************************************************************/
void Trace::SetTraceCFG(int MaxLen, double MaxStep, double InitStep,
						double MinStep, double ErrMax, double Delta,
						bool Verbose, int TraceDir) {
	
	/* set all of the params */					
	MaxLen_ = MaxLen;
	MaxStep_ = MaxStep;
	MinStep_ = MinStep;
	InitStep_ = InitStep;
	Verbose_ = Verbose;
	TraceDir_ = TraceDir;
	ErrMax_ = ErrMax;
	Delta_ = Delta;

	setTraceConfig_ = true;
}

/***********************************************************************
 * NAME : void Trace::SetTraceConfig()
 * 
 * DESCRIPTION : Sets a bunch of parameters for the trace.
 *		  
 * ********************************************************************/
void Trace::SetTraceCFG() {
	
	/* set default params */					
	MaxLen_ = 1000;
	MaxStep_ = 0.05;
	InitStep_ = 0.01;
	MinStep_ = 0.001;
	Verbose_ = false;
	TraceDir_ = 0;
	ErrMax_ = 0.0001;
	Delta_ = 0.05;

	setTraceConfig_ = true;
}


/***********************************************************************
 * NAME : void Trace::SetAlpha(nalpha,alpha)
 * 
 * DESCRIPTION : Set the alpha (polarization angles) to calculate 
 * 		h_alpha at.
 * 
 * INPUTS : 
 * 		int nalpha		The number of alpha values.
 * 		double *alpha	Array of alphas (degrees) 0 = toroidal.
 *	
 * OUPUTS :
 *		  
 * ********************************************************************/
void Trace::SetAlpha(int nalpha, double *alpha) {
	
	/*NOTE: for each alpha, there will be two traces - one for the 
	 * supplied value and one for alpha + 180 */
	/* set the alpha pointer */
	nalpha_ = nalpha;

	if (nalpha > 0) {
		alpha0_ = new double[nalpha_];
		alpha1_ = new double[nalpha_];
		allocAlpha_ = true;
		double dtor = M_PI/180.0;
		int i;
		for (i=0;i<nalpha;i++) {
			alpha0_[i] = alpha[i]*dtor;
			alpha1_[i] = fmod(alpha[i]*dtor + M_PI,2*M_PI);
		}
	}
	setAlpha_ = true;
}

/***********************************************************************
 * NAME : bool Trace::ContinueTrace(x,y,z,zmso,rmsm,rmso)
 * 
 * DESCRIPTION : Tells tracing routines whether to continue or not.
 * 	
 * INPUTS : 
 * 		double x		current x-coordinate.
 * 		double y		current y-coordinate.
 * 		double z		current z-coordinate.
 *	
 * OUPUTS :
 * 		double *zmso	z + 0.196
 *		double *rmsm	radial coordinate in MSM
 *		double *rmso	radial coordinate in MSO
 * 
 * RETURNS : 
 * 		bool cont		true if it wants to continue tracing.
 *		  
 * ********************************************************************/
bool Trace::ContinueTrace(	double x, double y, double z, 
							double *zmso, double *rmsm, double *rmso) {
	
	/* calculate zsmo, rmsm and rmso */
	zmso[0] = z + 0.196;
	double rho2 = x*x + y*y;
	rmsm[0] = sqrt(rho2 + z*z);
	rmso[0] = sqrt(rho2 + zmso[0]*zmso[0]);
	
	/* check if we're within MP, Tail and outside surface */
	bool goodMP, goodTX, goodSC;
	goodMP = ctfMP_(x,y,z,ktmodel.Rsm_);
	goodTX = ctfTX_(x,BoundTail_);
	goodSC = ctfSC_(z,rmsm[0],rmso[0]);
	
	return (goodMP && goodTX && goodSC);
}

/***********************************************************************
 * NAME : void Trace::Step(x0,y0,z0,step,x,y,z,zmso,rmsm,rmso,Bx,By,Bz)
 * 
 * DESCRIPTION : Take a single step using the Runge-Kutta-Merson algorithm.
 * 
 * INPUTS : 
 * 		double x0	Starting x-position
 * 		double y0	Starting y-position
 * 		double z0	Starting z-position
 *	
 * OUPUTS :
 * 		double *step	Step size.
 * 		double *x		New x-position
 * 		double *y		New y-position
 * 		double *z		New z-position
 * 		double *zmso	New z-position MSO
 * 		double *rmsm	Radial coordinate (MSM)
 * 		double *rmso	Radial coordinate (MSO)
 * 		double *Bx		x-component of magnetic field.
 * 		double *By		y-component of magnetic field.
 * 		double *Bz		z-component of magnetic field.
 *		  
 * ********************************************************************/
void Trace::Step(	double x0, double y0, double z0,
					double *step, 
					double *x, double *y, double *z,
					double *zmso, double *rmsm, double *rmso,
					double *Bx, double *By, double *Bz) {
	
	/* based on the STEP_08 function from GEOPACK */
	double rx1,ry1,rz1;
	double rx2,ry2,rz2;
	double rx3,ry3,rz3;
	double rx4,ry4,rz4;
	double rx5,ry5,rz5;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double x4,y4,z4;
	double step3;
	double Err;
	double xn, yn, zn, zmson, Rn, Rmson;
	bool cont;
	bool adjstep = false;
	bool repeat = true;	
	
	while (repeat) {
		/* this bit repeats until we get a desired step size */
		step3 = step[0]/3.0;
		StepVector(x0,y0,z0,step3,&rx1,&ry1,&rz1);
		x1 = x0 + rx1;
		y1 = y0 + ry1;
		z1 = z0 + rz1;
		StepVector(x1,y1,z1,step3,&rx2,&ry2,&rz2);
		x2 = x0 + 0.5*(rx1 + rx2);
		y2 = y0 + 0.5*(ry1 + ry2);
		z2 = z0 + 0.5*(rz1 + rz2);
		StepVector(x2,y2,z2,step3,&rx3,&ry3,&rz3);
		x3 = x0 + 0.375*(rx1 + 3*rx3);
		y3 = y0 + 0.375*(ry1 + 3*ry3);
		z3 = z0 + 0.375*(rz1 + 3*rz3);
		StepVector(x3,y3,z3,step3,&rx4,&ry4,&rz4);
		x4 = x0 + 1.5*(rx1 - 3*rx3 + 4*rx4);
		y4 = y0 + 1.5*(ry1 - 3*ry3 + 4*ry4);
		z4 = z0 + 1.5*(rz1 - 3*rz3 + 4*rz4);
		StepVector(x4,y4,z4,step3,&rx5,&ry5,&rz5);
		
		Err  = fabs(rx1 - 4.5*rx3 + 4*rx4 - 0.5*rx5);
		Err += fabs(ry1 - 4.5*ry3 + 4*ry4 - 0.5*ry5);
		Err += fabs(rz1 - 4.5*rz3 + 4*rz4 - 0.5*rz5);

		/* check that the next position is good */
		xn = x0 + 0.5*(rx1 + 4*rx4 + rx5);
		yn = y0 + 0.5*(ry1 + 4*ry4 + ry5);
		zn = z0 + 0.5*(rz1 + 4*rz4 + rz5);	
		cont = ContinueTrace(xn,yn,zn,&zmson,&Rn,&Rmson);
		if ((!cont) && (fabs(step[0]) > MinStep_)) {
			step[0] = 0.5*step[0];
			if (fabs(step[0]) < MinStep_) {
				step[0] = (step[0]/fabs(step[0]))*MinStep_;
			}
			adjstep = true;
		} else if (!cont) {
			step[0] = (step[0]/fabs(step[0]))*MinStep_;
			repeat = false;
		}
		
		if (cont) {
			if ((Err <= ErrMax_) && (fabs(step[0]) <= MaxStep_)) {
				repeat = false;
			} else {
				if (Err > ErrMax_) {
					if (fabs(step[0]) > MinStep_) {
						step[0] = step[0]*0.5;
					} else {
						repeat = false;
					}
				}
				if (fabs(step[0]) > MaxStep_) {
					step[0] = (step[0]/fabs(step[0]))*MaxStep_;
				}
			}
			
			if ((Err < 0.04*ErrMax_) && (fabs(step[0]) < (MaxStep_/1.5))) {
				step[0] = 1.5*step[0];
			}		
		}
	}
	
	x[0] = xn;
	y[0] = yn;
	z[0] = zn;	
	zmso[0] = zmson;
	rmsm[0] = Rn;
	rmso[0] = Rmson;
	
	ktmodel.Field(x[0],y[0],z[0],Bx,By,Bz);					
}

/***********************************************************************
 * NAME : void Trace::ReverseElements(n,x)
 * 
 * DESCRIPTION : Reverse the elements of an array.
 * 
 * INPUTS : 
 * 		int n		number of elements to be reversed.
 *		double *x	array to be reversed.
 *		  
 * ********************************************************************/
void Trace::ReverseElements(int n, double *x) {
	int i;
	double tmp;
	for (i=0;i<(n/2);i++) {
		tmp = x[i];
		x[i] = x[n-i-1];
		x[n-i-1] = tmp;
	}
}

/***********************************************************************
 * NAME : void Trace::StepVector(x,y,z,step3,rx,ry,rz)
 * 
 * DESCRIPTION : Calculate the vector of a given step.
 * 
 * INPUTS : 
 * 		double x		x-position
 * 		double y		y-position
 * 		double z		z-position
 * 		double step3	1/3 of the step size
 *	
 * OUPUTS :
 *		double *rx		output unit vector
 *		double *ry		output unit vector
 *		double *rz		output unit vector
 * 
 * ********************************************************************/
void Trace::StepVector(	double x, double y, double z, double step3,
						double *rx, double *ry, double *rz) {
	/* based on the RHAND_08 function from GEOPACK */
	
	double bx, by, bz, s3bm;
	ktmodel.Field(x,y,z,&bx,&by,&bz);
	
	s3bm = step3/sqrt(bx*bx + by*by + bz*bz); 
	/* this is a unit vector scaled by 1/3 of the step size */
	rx[0] = s3bm*bx;
	ry[0] = s3bm*by;
	rz[0] = s3bm*bz;
						
}

/***********************************************************************
 * NAME : void Trace::RKMTrace(x0,y0,z0,nstep,x,y,z,zmso,rmsm,rmso,
 * 								Bx,By,Bz)
 * 
 * DESCRIPTION : Perform the trace along the magnetic field line using
 * 				Runge-Kutta-Merson algorithm.
 * 
 * INPUTS : 
 * 		double x0		starting x-position
 * 		double y0		starting y-position
 * 		double z0		starting z-position
 *	
 * OUPUTS :
 *		int *nstep		Number of steps taken + 1 (actually the number 
 * 						of elements in each output trace.
 * 		double *x		Trace position
 * 		double *y		Trace position
 * 		double *z		Trace position
 * 		double *zmso	Trace z MSO coordinate.
 * 		double *Rmsm	MSM radial coordinate
 * 		double *Rmso	MSO radial coordinate
 * 		double *Bx		Trace field vector
 * 		double *By		Trace field vector
 * 		double *Bz		Trace field vector
 * 
 * ********************************************************************/
void Trace::RKMTrace(	double x0, double y0, double z0,
						int *nstep, 
						double *x, double *y, double *z,
						double *zmso, double *Rmsm, double *Rmso,
						double *Bx, double *By, double *Bz) {

	/* intialize the trace */
	nstep[0] = 1;
	x[0] = x0;
	y[0] = y0;
	z[0] = z0;
	ktmodel.Field(x0,y0,z0,&Bx[0],&By[0],&Bz[0]);
	double step;
	bool cont = ContinueTrace(x[0],y[0],z[0],&zmso[0],&Rmsm[0],&Rmso[0]);
	
	/* trace in one direction */
	if ((TraceDir_ == -1) || (TraceDir_ == 0)) {
		/* This will trace in the opposite of the field direction,
		 * towards the southern hemisphere*/ 
		step = -InitStep_;
		while ((cont) && (nstep[0] < (MaxLen_/2 - 1))) {
			Step(	x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],&step,
					&x[nstep[0]],&y[nstep[0]],&z[nstep[0]],
					&zmso[nstep[0]],&Rmsm[nstep[0]],&Rmso[nstep[0]],
					&Bx[nstep[0]],&By[nstep[0]],&Bz[nstep[0]]);
			cont = ContinueTrace(x[nstep[0]],y[nstep[0]],z[nstep[0]],
								&zmso[nstep[0]],&Rmsm[nstep[0]],&Rmso[nstep[0]]);
			nstep[0]++;
		}
	}
	
	/* reverse the elements of the trace */
	ReverseElements(nstep[0],x);
	ReverseElements(nstep[0],y);
	ReverseElements(nstep[0],z);
	ReverseElements(nstep[0],zmso);
	ReverseElements(nstep[0],Bx);
	ReverseElements(nstep[0],By);
	ReverseElements(nstep[0],Bz);
	ReverseElements(nstep[0],Rmsm);
	ReverseElements(nstep[0],Rmso);
	
	/* trace in the opposite direction */
	cont = ContinueTrace(x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],
						&zmso[nstep[0]-1],&Rmsm[nstep[0]-1],&Rmso[nstep[0]-1]);
	if ((TraceDir_ == 1) || (TraceDir_ == 0)) {
		/* hopefully this will go in the direction fo the field vectors
		 * towards the northern hemisphere */
		step = InitStep_;
		while ((cont) && (nstep[0] < (MaxLen_ - 1))) {
			Step(	x[nstep[0]-1],y[nstep[0]-1],z[nstep[0]-1],&step,
					&x[nstep[0]],&y[nstep[0]],&z[nstep[0]],
					&zmso[nstep[0]],&Rmsm[nstep[0]],&Rmso[nstep[0]],
					&Bx[nstep[0]],&By[nstep[0]],&Bz[nstep[0]]);
			cont = ContinueTrace(x[nstep[0]],y[nstep[0]],z[nstep[0]],
								&zmso[nstep[0]],&Rmsm[nstep[0]],&Rmso[nstep[0]]);
			nstep[0]++;
		}
	}
	
}

/***********************************************************************
 * NAME : Trace Trace::TracePostion(i,x,y,z)
 * 
 * DESCRIPTION : Return a trace object for  a position near to an 
 * 			existing trace.
 * 
 * INPUTS : 
 *		int i		Index of existing trace.
 * 		double x	x-coordinate near to existing trace
 * 		double y	y-coordinate near to existing trace
 * 		double z	z-coordinate near to existing trace
 * 		
 * RETURNS :
 * 		Trace T		Trace object instance.
 *		  
 * ********************************************************************/
Trace Trace::TracePosition(int i, double x, double y, double z) {
	/* return a new trace object at the supplied position using the
	 * parameters at time i */
	Trace T;
	
	/* input position and time - I am pretty certain that the midpoints
	 * of the field lines are stored in SM coords */
	T.InputPos(1,&x,&y,&z);
	
	/* set the model up */
	T.SetTraceBounds(BoundMP_,BoundTail_,BoundSurface_);
	T.SetModelParams(&Rsm_[i],&t1_[i],&t2_[i]);
	T.SetTraceCFG(MaxLen_,MaxStep_,InitStep_,MinStep_,ErrMax_,Delta_,false,0);
	
	/* run the trace */
	T.TraceField();

	/* calculate S*/
	T.CalculateTraceDist();
	
	return T;
	
}

/***********************************************************************
 * NAME : void Trace::_CalculateTraceHalpha(i,j,halpha)
 * 
 * DESCRIPTION : Calculates the h_alphas along a given trace and alpha.
 * 
 * INPUTS : 
 *		int i		Trace index
 * 		int j		Alpha index
 * 		
 * OUPUTS :
 *		double *halpha	h_alpha along this field line.
 *   
 * ********************************************************************/
void Trace::_CalculateTraceHalpha(	int i, int j, double *halpha) {

	/* some variables needed */
	double xe0,ye0,ze0,xe1,ye1,ze1;
	int k;
	
	/* get the trace starting points first */
	_CalculateHalphaStartPoints(i,j,&xe0,&ye0,&ze0,&xe1,&ye1,&ze1);

	/* do two traces */
	Trace T0 = TracePosition(i,xe0,ye0,ze0);
	Trace T1 = TracePosition(i,xe1,ye1,ze1);
	
	/* the traces above may only have 0 steps - in which case we can 
	 * just fill halpha with nans and leave the function */
	if ((T0.nstep_[0] == 0) | (T1.nstep_[0] == 0)) {
		for (k=0;k<nstep_[i];k++) {
			halpha[k] = NAN;
		}
		return;
	}
	
	/* get the closest points to each step of the original trace*/
	double *xc0 = new double[nstep_[i]];
	double *yc0 = new double[nstep_[i]];
	double *zc0 = new double[nstep_[i]];
	double *xc1 = new double[nstep_[i]];
	double *yc1 = new double[nstep_[i]];
	double *zc1 = new double[nstep_[i]];

	interptraceClosestPos(	nstep_[i],x_[i],y_[i],z_[i],
							bx_[i],by_[i],bz_[i],
							T0.nstep_[0],T0.x_[0],T0.y_[0],T0.z_[0],T0.S_[0],
							T1.nstep_[0],T1.x_[0],T1.y_[0],T1.z_[0],T1.S_[0],
							xc0,yc0,zc0,xc1,yc1,zc1);

	/* calculate distances and then halpha */
	double d, dx, dy, dz, h0, h1;
	
	for (k=0;k<nstep_[i];k++) {
		dx = x_[i][k] - xc0[k];
		dy = y_[i][k] - yc0[k];
		dz = z_[i][k] - zc0[k];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		h0 = d/Delta_;
		
		dx = x_[i][k] - xc1[k];
		dy = y_[i][k] - yc1[k];
		dz = z_[i][k] - zc1[k];
		d = sqrt(dx*dx + dy*dy + dz*dz);
		h1 = d/Delta_;
		
		halpha[k] = 0.5*(h0 + h1);
	}

	/* free up memory */
	delete[] xc0;
	delete[] yc0;
	delete[] zc0;
	delete[] xc1;
	delete[] yc1;
	delete[] zc1;
}

/***********************************************************************
 * NAME : void Trace::_CalculateHalpha()
 * 
 * DESCRIPTION : Loop through all traces and their alpha values to 
 * 			calculate h_alpha.
 *		  
 * ********************************************************************/
void Trace::_CalculateHalpha() {

	/* loop through each trace and alpha combination */
	int i, j, k, I, J;
	for (i=0;i<n_;i++) {
		I = i*(nalpha_*MaxLen_);
		if (isfinite(FP_[i][12])) {
			for (j=0;j<nalpha_;j++) {
				J = j*MaxLen_;
				_CalculateTraceHalpha(i,j,Halpha3D_[i][j]);
				for (k=0;k<MaxLen_;k++) {
					Halpha_[I + J + k] = Halpha3D_[i][j][k];
				}
			}
		}
	}
	setHalpha_ = true;
}

/***********************************************************************
 * NAME : bool Trace::_CheckHalpha()
 * 
 * DESCRIPTION : Checks that the object is set up correctly before it
 * 			attempts to calculate h_alpha.
 * 
 * RETURNS :
 * 		bool	true if it is ready.
 *		  
 * ********************************************************************/		
bool Trace::_CheckHalpha() {
	
	if ((!allocAlpha_) || (!setAlpha_)) {
		printf("Run the 'SetAlpha()' function prior to calculating h_alpha\n");
		return false;
	}
	
	if (nalpha_ <= 0) {
		printf("1 or more values of alpha must be provided to calculate h_alpha\n");
		return false;
	}

	if (setHalpha_) {
		printf("H alpha already calculated\n");
		return false;
	}
	
	return true;
	
}

/***********************************************************************
 * NAME : void Trace::CalculateHalpha()
 * 
 * DESCRIPTION : Calculate the h_alpha for all traces.
 * 
 *		  
 * ********************************************************************/		
void Trace::CalculateHalpha() {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate both 1D and 3D arrays */
	Halpha_ = new double[n_*nalpha_*MaxLen_];
	Halpha3D_ = new double**[n_];
	int i, j;
	for (i=0;i<n_;i++) {
		Halpha3D_[i] = new double*[nalpha_];
		for (j=0;j<nalpha_;j++) {
			Halpha3D_[i][j] = new double[MaxLen_];
		}
	}
	allocHalpha_ = true;
	allocHalpha3D_ = true;

	_CalculateHalpha();
}

/***********************************************************************
 * NAME : void Trace::CalculateHalpha()
 * 
 * DESCRIPTION : Calculate the h_alpha for all traces.
 *	
 * OUPUTS :
 *		 double *halpha	1-D array containing all trace h_alphas
 * 
 * ********************************************************************/
void Trace::CalculateHalpha(double *halpha) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate 3D array and use pointer for 1D */
	Halpha_ = halpha;
	Halpha3D_ = new double**[n_];
	int i, j;
	for (i=0;i<n_;i++) {
		Halpha3D_[i] = new double*[nalpha_];
		for (j=0;j<nalpha_;j++) {
			Halpha3D_[i][j] = new double[MaxLen_];
		}
	}
	allocHalpha3D_ = true;
	_CalculateHalpha();
}

/***********************************************************************
 * NAME : void Trace::CalculateHalpha()
 * 
 * DESCRIPTION : Calculate the h_alpha for all traces.
 * 
 * OUPUTS :
 *		double ***hapha3d	3D array containing all of the h_alphas  
 * 
 * ********************************************************************/
void Trace::CalculateHalpha(double ***halpha3d) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* allocate 1D and use pointer for 3D array */
	Halpha_ = new double[n_*nalpha_*MaxLen_];
	Halpha3D_ = halpha3d;

	allocHalpha_ = true;
	
	_CalculateHalpha();
}

/***********************************************************************
 * NAME : void Trace::CalculateHalpha()
 * 
 * DESCRIPTION : Calculate the h_alpha for all traces.
 *	
 * OUPUTS :
 * 		double *halpha	1-D array containing all trace h_alphas
 * 		double ***hapha3d	3D array containing all of the h_alphas  
 *		  
 * ********************************************************************/
void Trace::CalculateHalpha(double *halpha, double ***halpha3d) {

	if (!_CheckHalpha()) {
		return;
	}
	
	/* use pointer for both 1D and 3D arrays */
	Halpha_ = halpha;
	Halpha3D_ = halpha3d;
	
	_CalculateHalpha();
}

/***********************************************************************
 * NAME : void Trace::_CalculateHalphaStartPoints(i,j,xe0,ye0,ze0,xe1,ye1,ze1)
 * 
 * DESCRIPTION : Calcualte the starting points for two adjacent traces
 * 		to one original trace (one it 180 degrees around from the other)
 * 		in order to calculate h_alpha.
 * 
 * INPUTS : 
 * 		int i		Trace index
 * 		int j		Alpha index
 *	
 * OUPUTS :
 *		double *xe0		Starting position for adjacent field line 0 
 *		double *ye0		Starting position for adjacent field line 0 
 *		double *ze0		Starting position for adjacent field line 0 
 *		double *xe1		Starting position for adjacent field line 1 
 *		double *ye1		Starting position for adjacent field line 1 
 *		double *ze1		Starting position for adjacent field line 1 
 * 
 * ********************************************************************/
void Trace::_CalculateHalphaStartPoints(int i, int j,
							double *xe0, double *ye0, double *ze0,
							double *xe1, double *ye1, double *ze1) {
	
	/* calculate the tracing start points for each alpha */
	double dt, dp, beta, dx, dy;
	
	/* dt and dp are the toroidal and poloidal components of Delta */
	dt = Delta_*cos(alpha0_[j]); // alpha = 0.0 is toroidal
	dp = Delta_*sin(alpha0_[j]);
	
	/* rotate based on the local time */
	beta = atan2(-xfe_[i],-yfe_[i]);
	dy = dp*cos(beta) - dt*sin(beta);
	dx = dp*sin(beta) + dt*cos(beta);
	
	/* set the start points of the new field lines */
	xe0[0] = xfe_[i] + dx;
	ye0[0] = yfe_[i] + dy;
	ze0[0] = zfe_[i];
	xe1[0] = xfe_[i] - dx;
	ye1[0] = yfe_[i] - dy;
	ze1[0] = zfe_[i];

}

/***********************************************************************
 * NAME : void Trace::_AllocTrace()
 * 
 * DESCRIPTION : Allocates arrays to store trace positions and fields
 * 
 * ********************************************************************/
void Trace::_AllocTrace() {
	x_ = new double*[n_];					
	y_ = new double*[n_];					
	z_ = new double*[n_];					
	bx_ = new double*[n_];					
	by_ = new double*[n_];					
	bz_ = new double*[n_];
	int i;
	for (i=0;i<n_;i++) {
		x_[i] = new double[MaxLen_];					
		y_[i] = new double[MaxLen_];					
		z_[i] = new double[MaxLen_];					
		bx_[i] = new double[MaxLen_];					
		by_[i] = new double[MaxLen_];					
		bz_[i] = new double[MaxLen_];		
	}		
	allocTrace_ = true;
}

/***********************************************************************
 * NAME : void Trace::_AllocTraceR()
 * 
 * DESCRIPTION : Allocates Rmsm and Rmso arrays.
 *		  
 * ********************************************************************/
void Trace::_AllocTraceR() {

	Rmsm_ = new double*[n_];
	Rmso_ = new double*[n_];
	int i;
	for (i=0;i<n_;i++) {
		Rmsm_[i] = new double[MaxLen_];	
		Rmso_[i] = new double[MaxLen_];
	}	
	allocR_ = true;
}

/***********************************************************************
 * NAME : void Trace::_AllocZmso()
 * 1
 * DESCRIPTION : Allocates array to store zmso.
 *		  
 * ********************************************************************/
void Trace::_AllocZmso() {

	zmso_ = new double*[n_];
	int i;
	for (i=0;i<n_;i++) {
		zmso_[i] = new double[MaxLen_];		
	}
	allocZmso_ = true;
}

/***********************************************************************
 * NAME : void Trace::TraceField(nstep,x,y,z,Rmsm,Rmso,bx,by,bz)
 * 
 * DESCRIPTION : Trace the field lines.
 *	
 * OUTPUTS :
 * 		int *nstep		Number of trace steps
 * 		double **x		Trace positions
 * 		double **y		Trace positions
 * 		double **z		Trace positions
 * 		double **Rmsm	Trace radial coordinates (MSM)
 * 		double **Rmso	Trace radial coordinates (MSO)
 * 		double **Bx		Trace field
 * 		double **By		Trace field
 * 		double **Bz		Trace field
 *		  
 * ********************************************************************/
void Trace::TraceField(	int *nstep,
						double **x, double **y, double **z,
						double **Rmsm, double **Rmso,
						double **bx, double **by, double **bz) {
	
	/* link the pointers within the object to those supplied by this 
	 * function					*/
	nstep_ = nstep;
	x_ = x;					
	y_ = y;					
	z_ = z;					
	bx_ = bx;					
	by_ = by;					
	bz_ = bz;		
	Rmsm_ = Rmsm;
	Rmso_ = Rmso;
	_AllocZmso();
	
	/* call the tracing code */
	_TraceField();
}

/***********************************************************************
 * NAME : void Trace::TraceField(nstep,x,y,z,bx,by,bz)
 * 
 * DESCRIPTION : Trace the field lines.
 *	
 * OUPUTS :
 * 		int *nstep		Number of trace steps
 * 		double **x		Trace positions
 * 		double **y		Trace positions
 * 		double **z		Trace positions
 * 		double **Bx		Trace field
 * 		double **By		Trace field
 * 		double **Bz		Trace field
 *		  
 * ********************************************************************/
void Trace::TraceField(	int *nstep,
						double **x, double **y, double **z,
						double **bx, double **by, double **bz) {
	
	/* link the pointers within the object to those supplied by this 
	 * function					*/
	nstep_ = nstep;
	x_ = x;					
	y_ = y;					
	z_ = z;					
	bx_ = bx;					
	by_ = by;					
	bz_ = bz;		
	_AllocTraceR();
	_AllocZmso();
	
	/* call the tracing code */
	_TraceField();
}

/***********************************************************************
 * NAME : void Trace::TraceField(nstep)
 * 
 * DESCRIPTION : Trace the field lines.
 *	
 * OUPUTS :
 * 		int *nstep		Number of trace steps
 *		  
 * ********************************************************************/
void Trace::TraceField(	int *nstep) {
	
	/* link the pointers within the object to those supplied by this 
	 * function					*/
	nstep_ = nstep;
	_AllocTrace();	
	_AllocTraceR();	
	_AllocZmso();
	
	/* call the tracing code */
	_TraceField();	

}

/***********************************************************************
 * NAME : void Trace::TraceField()
 * 
 * DESCRIPTION : Trace the field lines.
 *		  
 * ********************************************************************/
void Trace::TraceField() {
	
	/* no pointers provided: allocate them*/
	if (!allocNstep_) {
		nstep_ = new int[n_];
		allocNstep_ = true;
	}
	_AllocTrace();
	_AllocTraceR();
	_AllocZmso();
	
	/* call the tracing code */
	_TraceField();
	
}

/***********************************************************************
 * NAME : void Trace::_TraceField()
 * 
 * DESCRIPTION : Run the field traces.
 *		  
 * ********************************************************************/
void Trace::_TraceField() {
	
	/* this function actually calls the tracing routines */
	
	/* check this hasn't already been done */
	if (setTrace_) {
		printf("Attempted to trace twice? not happening mate...\n");
		return;
	}
	
	/* check we have input positions */
	if (!allocInputPos_) {
		printf("Need InputPos() before trace\n");
		return;
	}
	
	
	/* check that we have model parameters */
	if (!setModelParams_) {
		printf("Run SetModelParams() before tracing\n");
		return;
	}
	

	
	/* check if all of the starting points are within the MP */
	inMP_ = new bool[n_];
	allocMP_ = true;
	int i;
	for (i=0;i<n_;i++) {
		inMP_[i] = WithinMP(x0_[i],y0_[i],z0_[i],Rsm_[i]);
	}
	
	
	for (i=0;i<n_;i++) {
		if (Verbose_) {
			printf("\rTracing field line %d of %d (%6.2f)%%",i+1,n_,((float) (i+1)*100.0)/n_);
		}

		
		if (inMP_[i]) {
			/* set current trace parameters */
			ktmodel.SetParams(Rsm_[i],t1_[i],t2_[i]);
		
			/* perform trace */
			RKMTrace(	x0_[i],y0_[i],z0_[i],&nstep_[i],
						x_[i],y_[i],z_[i],
						zmso_[i],Rmsm_[i],Rmso_[i],
						bx_[i],by_[i],bz_[i]);

		} else {
			/*fill with NaN*/
			nstep_[i] = 0;
		}
							
	}	
	if (Verbose_) { 
		printf("\n");
	}
	setTrace_ = true;
}


/***********************************************************************
 * NAME : void Trace::CalculateTraceDist()
 * 
 * DESCRIPTION : Calculate the distances along each field line.
 *		  
 * ********************************************************************/
void Trace::CalculateTraceDist() {
	int i;
	S_ = new double*[n_];
	for (i=0;i<n_;i++) {
		S_[i] = new double[MaxLen_];
	}
	allocDist_ = true;
	
	_CalculateTraceDist();
}

/***********************************************************************
 * NAME : void Trace::CalculateTraceDist(S)
 * 
 * DESCRIPTION : Calculate the distances along each field line.
 *	
 * OUPUTS :
 * 		double **S		Distance along each field line.
 *		  
 * ********************************************************************/
void Trace::CalculateTraceDist(double **S) {

	S_ = S;
	_CalculateTraceDist();
}

/***********************************************************************
 * NAME : void Trace_CalcualteTraceDist()
 * 
 * DESCRIPTION : Calcualtes the distance along all field lines.
 *		  
 * ********************************************************************/
void Trace::_CalculateTraceDist() {
	int i, j;
	double dx, dy, dz;
	for (i=0;i<n_;i++) {
		S_[i][0] = 0.0;
		for (j=1;j<nstep_[i];j++) {
			dx = x_[i][j] - x_[i][j-1];
			dy = y_[i][j] - y_[i][j-1];
			dz = z_[i][j] - z_[i][j-1];
			S_[i][j] = S_[i][j-1] + sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
	setDist_ = true;
}


/***********************************************************************
 * NAME : void Trace::CalculateTraceRnorm()
 * 
 * DESCRIPTION : Calcualtes the Rnorm value along each trace.
 *		  
 * ********************************************************************/
void Trace::CalculateTraceRnorm() {
	int i;
	Rnorm_ = new double*[n_];
	for (i=0;i<n_;i++) {
		Rnorm_[i] = new double[MaxLen_];
	}
	allocRnorm_ = true;
	
	_CalculateTraceRnorm();
}

/***********************************************************************
 * NAME : void Trace::CalculateTraceRnorm(Rnorm)
 * 
 * DESCRIPTION : Calcualtes the Rnorm value along each trace.
 * 
 * OUPUTS :
 *		double **Rnorm		Rnorm along each trace.  
 * 
 * ********************************************************************/
void Trace::CalculateTraceRnorm(double **Rnorm) {
	
	Rnorm_ = Rnorm;
	
	_CalculateTraceRnorm();
}

/***********************************************************************
 * NAME : void Trace::_CalculateTraceRnorm()
 * 
 * DESCRIPTION : Calcualtes Rnorm values for all traces.
 *		  
 * ********************************************************************/
void Trace::_CalculateTraceRnorm() {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			/* need footprints and R done first */
			FieldLineRnorm(nstep_[i],Rmsm_[i],FP_[i][16],Rnorm_[i]);
		}
	}
	setRnorm_ = true;
}

/***********************************************************************
 * NAME : void Trace::CalculateTraceFP()
 * 
 * DESCRIPTION : Works out where all of the footprints are for each 
 * 				trace.
 * 
 *		  
 * ********************************************************************/
void Trace::CalculateTraceFP() {
	int i;
	FP_ = new double*[n_];
	for (i=0;i<n_;i++) {
		FP_[i] = new double[18];
	}
	allocFootprints_ = true;
	
	_CalculateTraceFP();
}

/***********************************************************************
 * NAME : void Trace::CalculateTraceFP()
 * 
 * DESCRIPTION : Works out where all of the footprints are for each 
 * 				trace.
 *	
 * OUPUTS :
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
 * 
 *		  
 * ********************************************************************/
void Trace::CalculateTraceFP(double **FP) {
	
	FP_ = FP;
	
	_CalculateTraceFP();
}

/***********************************************************************
 * NAME : void Trace::_CalculateTraceFP()
 * 
 * DESCRIPTION : Calculates all footprints for each trace.
 *		  
 * ********************************************************************/
void Trace::_CalculateTraceFP() {

	/* before running this we should check that the following functions 
	 * have been called:
	 * 1. TraceField()
	 * 2. CalcualteTraceDist() */
	if (!setTrace_) {
		printf("Call TraceField() before calculating footprints\n");
		return;
	}
	if (!setDist_) {
		printf("Call CalcualteTraceDist() before calculating footprints\n");
		return;
	}
	
	/* allocate the endpoints */
	xfn_ = new double[n_];
	yfn_ = new double[n_];
	zfn_ = new double[n_];
	xfs_ = new double[n_];
	yfs_ = new double[n_];
	zfs_ = new double[n_];
	xfnc_ = new double[n_];
	yfnc_ = new double[n_];
	zfnc_ = new double[n_];
	xfsc_ = new double[n_];
	yfsc_ = new double[n_];
	zfsc_ = new double[n_];
	xfnv_ = new double[n_];
	yfnv_ = new double[n_];
	zfnv_ = new double[n_];
	xfsv_ = new double[n_];
	yfsv_ = new double[n_];
	zfsv_ = new double[n_];
	xfnvc_ = new double[n_];
	yfnvc_ = new double[n_];
	zfnvc_ = new double[n_];
	xfsvc_ = new double[n_];
	yfsvc_ = new double[n_];
	zfsvc_ = new double[n_];
	allocEndpoints_ = true;
	
	xfe_ = new double[n_];
	yfe_ = new double[n_];
	zfe_ = new double[n_];
	allocEqFP_ = true;
	
	int i, j;
	for (i=0;i<n_;i++) {
		_SingleTraceFP(i);
							
	}
	
	/* calculate the lats/ locals times for each footprint */
	_FPCoords();
	
	setFootprints_ = true;
}

/***********************************************************************
 * NAME : void Trace::FPCoords()
 * 
 * DESCRIPTION : Converts Cartesian footprint positions to latitude,
 * 			local times and L shell.
 *		  
 * ********************************************************************/
void Trace::_FPCoords() {
	
	int i;
	for (i=0;i<n_;i++) {
		/* surface */
		LatLT(xfn_[i],yfn_[i],zfn_[i],&FP_[i][0],&FP_[i][1]);
		LatLT(xfs_[i],yfs_[i],zfs_[i],&FP_[i][2],&FP_[i][3]);
		/* core */
		LatLT(xfnc_[i],yfnc_[i],zfnc_[i],&FP_[i][4],&FP_[i][5]);
		LatLT(xfsc_[i],yfsc_[i],zfsc_[i],&FP_[i][6],&FP_[i][7]);
		/* dipole R = 1 */
		LatLT(xfnv_[i],yfnv_[i],zfnv_[i],&FP_[i][8],&FP_[i][9]);
		LatLT(xfsv_[i],yfsv_[i],zfsv_[i],&FP_[i][10],&FP_[i][11]);
		/* dipole R = 0.832 */
		LatLT(xfnvc_[i],yfnvc_[i],zfnvc_[i],&FP_[i][12],&FP_[i][13]);
		LatLT(xfsvc_[i],yfsvc_[i],zfsvc_[i],&FP_[i][14],&FP_[i][15]);
		/* equatorial */
		LshellMLT(xfe_[i],yfe_[i],zfe_[i],&FP_[i][16],&FP_[i][17]);
	}
	
}

/***********************************************************************
 * NAME : void Trace::_InterpPos(xi,yi,zi,s,t,target,xo,yo,zo)
 * 
 * DESCRIPTION : Attempt to work out a position along a trace which
 * 		corresponds to some target value along a target array (t) using
 * 		linear interpolation.
 * 		
 * 		e.g. For the northern footprint on the planet, we provide
 * 		the two positions surrounding Rmso = 1, their corresponding 
 * 		trace distances (s), target arrays (R) and target value (R=1):
 * 		_InterpPos(xi,yi,zi,s,R,1.0,xo,yo,zo)
 * 		
 * 		It works out what trace distance corresponds to t == target 
 * 		first, then uses interpolation to work out what x,y,z 
 * 		corresponds to that same distance.
 * 
 * 
 * INPUTS : 
 * 		double *xi		Position along trace around target.
 * 		double *yi		Position along trace around target.
 * 		double *zi		Position along trace around target.
 * 		double *s		distance along field line around target.
 * 		double *t		target array.
 * 		double *target 	target values (~between t[0] and t[1])
 *	
 * OUPUTS :
 * 		double *x0		Output position
 * 		double *y0		Output position
 * 		double *z0		Output position
 *		  
 * ********************************************************************/
void Trace::_InterpPos(	double *xi, double *yi, double *zi,
						double *s, double *t, double target,
						double *xo, double *yo, double *zo) {
	
	/* calculate the gradients wrt to t and s */
	double s_targ, ds;
	double m, mx, my, mz;
	double c, cx, cy, cz;
	ds = (s[1] - s[0]);
	m = ds/(t[1] - t[0]);
	c = s[0] - m*t[0];
	mx = (xi[1] - xi[0])/ds;
	cx = xi[0] - mx*s[0];
	my = (yi[1] - yi[0])/ds;
	cy = yi[0] - my*s[0];
	mz = (zi[1] - zi[0])/ds;
	cz = zi[0] - mz*s[0];
		
	/* this is the distance along the field line where our target
	 * array (t) crosses the target value (target) */
	s_targ = c + target*m;
	
	/* calculate target crossing footprints */
	xo[0] = cx + mx*s_targ;
	yo[0] = cy + my*s_targ;
	zo[0] = cz + mz*s_targ;

}

/***********************************************************************
 * NAME : void Trace::_SingleTraceFP(I)
 * 
 * DESCRIPTION : Calculate the footprints for a trace.
 * 
 * INPUTS : 
 * 		int I		Trace index.
 *		  
 * ******************************************!**************************/
void Trace::_SingleTraceFP(	int I) {

	int i;
	int eq_ind = -1;
	int n_ind = -1;
	int s_ind = -1;
	int nv_ind = -1;
	int sv_ind = -1;
	int nc_ind = -1;
	int sc_ind = -1;
	int nvc_ind = -1;
	int svc_ind = -1;

	/* check inMP */
	if ((inMP_[I]) && (nstep_[I] > 1)) {

		
		/* check if there's an equator crossing */
		for (i=0;i<nstep_[I]-1;i++) {
			if ((z_[I][i+1] >= 0) && (z_[I][i] < 0)) {
				eq_ind = i;
				break;
			}
		}
		
	
		/* southern footprints */
		for (i=0;i<nstep_[I]-1;i++) {
			/* check that we're still in the south */
			if ((Rmsm_[I][i+1] < Rmsm_[I][i]) || (Rmsm_[I][i] >= 2.0)) {
				break;
			}
			
			/* find where we cross Rmso == 1.0 */
			if ((Rmso_[I][i] <= 1.0) && (Rmso_[I][i+1] > 1.0)) {
				s_ind = i;
			}
			
			/* find where we cross Rmso == 0.832 */
			if ((Rmso_[I][i] <= 0.832) && (Rmso_[I][i+1] > 0.832)) {
				sc_ind = i;
			}
		
			/* find where we cross Rmsm == 1.0 */
			if ((Rmsm_[I][i] <= 1.0) && (Rmsm_[I][i+1] > 1.0)) {
				sv_ind = i;
			}
			
			/* find where we cross Rmsm == 0.832 */
			if ((Rmsm_[I][i] <= 0.832) && (Rmsm_[I][i+1] > 0.832)) {
				svc_ind = i;
			}
		}
		
		/* northern footprints */
		for (i=nstep_[I]-2;i>=0;i--) {
			/* check that we're still in the south */
			if ((Rmsm_[I][i+1] > Rmsm_[I][i]) || (Rmsm_[I][i] >= 2.0)) {
				break;
			}
			
			/* find where we cross Rmso == 1.0 */
			if ((Rmso_[I][i] > 1.0) && (Rmso_[I][i+1] <= 1.0)) {
				n_ind = i;
			}
			
			/* find where we cross Rmso == 0.832 */
			if ((Rmso_[I][i] > 0.832) && (Rmso_[I][i+1] <= 0.832)) {
				nc_ind = i;
			}
		
			/* find where we cross Rmsm == 1.0 */
			if ((Rmsm_[I][i] > 1.0) && (Rmsm_[I][i+1] <= 1.0)) {
				nv_ind = i;
			}
			
			/* find where we cross Rmsm == 0.832 */
			if ((Rmsm_[I][i] > 0.832) && (Rmsm_[I][i+1] <= 0.832)) {
				nvc_ind = i;
			}
		}
		

	}
	
	/* set equatorial footprint */
	if (eq_ind >= 0) {
		/* interpolate for the point where z = 0*/
		_InterpPos(	&x_[I][eq_ind],&y_[I][eq_ind],&z_[I][eq_ind],
					&S_[I][eq_ind],&z_[I][eq_ind],0.0,
					&xfe_[I],&yfe_[I],&zfe_[I]);		
	} else {
		xfe_[I] = NAN;
		yfe_[I] = NAN;
		zfe_[I] = NAN;
	}
	
	/* northern surface */
	if (n_ind >= 0) {
		_InterpPos(	&x_[I][n_ind],&y_[I][n_ind],&zmso_[I][n_ind],
					&S_[I][n_ind],&Rmso_[I][n_ind],1.0,
					&xfn_[I],&yfn_[I],&zfn_[I]);	
	} else {
		xfn_[I] = NAN;
		yfn_[I] = NAN;
		zfn_[I] = NAN;
	}

	
	/* southern surface */
	if (s_ind >= 0) {
		_InterpPos(	&x_[I][s_ind],&y_[I][s_ind],&zmso_[I][s_ind],
					&S_[I][s_ind],&Rmso_[I][s_ind],1.0,
					&xfs_[I],&yfs_[I],&zfs_[I]);	
	} else {
		xfs_[I] = NAN;
		yfs_[I] = NAN;
		zfs_[I] = NAN;
	}
	
	/* northern core */
	if (nc_ind >= 0) {
		_InterpPos(	&x_[I][nc_ind],&y_[I][nc_ind],&zmso_[I][nc_ind],
					&S_[I][nc_ind],&Rmso_[I][nc_ind],0.832,
					&xfnc_[I],&yfnc_[I],&zfnc_[I]);	
	} else {
		xfnc_[I] = NAN;
		yfnc_[I] = NAN;
		zfnc_[I] = NAN;
	}

	
	/* southern core */
	if (sc_ind >= 0) {
		_InterpPos(	&x_[I][sc_ind],&y_[I][sc_ind],&zmso_[I][sc_ind],
					&S_[I][sc_ind],&Rmso_[I][sc_ind],0.832,
					&xfsc_[I],&yfsc_[I],&zfsc_[I]);	
	} else {
		xfsc_[I] = NAN;
		yfsc_[I] = NAN;
		zfsc_[I] = NAN;
	}
	
	/* northern dipole surface */
	if (nv_ind >= 0) {
		_InterpPos(	&x_[I][nv_ind],&y_[I][nv_ind],&z_[I][nv_ind],
					&S_[I][nv_ind],&Rmsm_[I][nv_ind],1.0,
					&xfnv_[I],&yfnv_[I],&zfnv_[I]);	
	} else {
		xfnv_[I] = NAN;
		yfnv_[I] = NAN;
		zfnv_[I] = NAN;
	}

	
	/* southern dipole surface */
	if (sv_ind >= 0) {
		_InterpPos(	&x_[I][sv_ind],&y_[I][sv_ind],&z_[I][sv_ind],
					&S_[I][sv_ind],&Rmsm_[I][sv_ind],1.0,
					&xfsv_[I],&yfsv_[I],&zfsv_[I]);	
	} else {
		xfsv_[I] = NAN;
		yfsv_[I] = NAN;
		zfsv_[I] = NAN;
	}
	
	/* northern dipole core */
	if (nvc_ind >= 0) {
		_InterpPos(	&x_[I][nvc_ind],&y_[I][nvc_ind],&z_[I][nvc_ind],
					&S_[I][nvc_ind],&Rmsm_[I][nvc_ind],0.832,
					&xfnvc_[I],&yfnvc_[I],&zfnvc_[I]);	
	} else {
		xfnvc_[I] = NAN;
		yfnvc_[I] = NAN;
		zfnvc_[I] = NAN;
	}

	
	/* southern dipole core */
	if (svc_ind >= 0) {
		_InterpPos(	&x_[I][svc_ind],&y_[I][svc_ind],&z_[I][svc_ind],
					&S_[I][svc_ind],&Rmsm_[I][svc_ind],0.832,
					&xfsvc_[I],&yfsvc_[I],&zfsvc_[I]);	
	} else {
		xfsvc_[I] = NAN;
		yfsvc_[I] = NAN;
		zfsvc_[I] = NAN;
	}
	
}

/***********************************************************************
 * NAME : void Trace::GetTrace(x,y,z)
 * 
 * DESCRIPTION : Get trace positions.
 *	
 * OUPUTS :
 * 		double **x		Positions along traces.
 * 		double **y		Positions along traces.
 * 		double **z		Positions along traces.
 *		  
 * ********************************************************************/
void Trace::GetTrace(double **x,double **y, double **z) {
	/* copy GSM position into output arrays*/
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			x[i][j] = x_[i][j];
			y[i][j] = y_[i][j];
			z[i][j] = z_[i][j];
		}
	}
}

/***********************************************************************
 * NAME : void Trace::GetTrace(x,y,z)
 * 
 * DESCRIPTION : Get trace positions and field vectors.
 *	
 * OUPUTS :
 * 		double **x		Positions along traces.
 * 		double **y		Positions along traces.
 * 		double **z		Positions along traces.
 * 		double **Bx		Field along traces
 * 		double **By		Field along traces
 * 		double **Bz		Field along traces
 *		  
 * ********************************************************************/
void Trace::GetTrace(	double **x,double **y, double **z,
						double **Bx,double **By, double **Bz) {
	/* copy GSM field into output arrays*/
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Bx[i][j] = bx_[i][j];
			By[i][j] = by_[i][j];
			Bz[i][j] = bz_[i][j];
		}
	}
	
	/* get the position */
	GetTrace(x,y,z);
}

/***********************************************************************
 * NAME : void Trace::GetTraceDist(S)
 * 
 * DESCRIPTION : Get the distance along the field traces.
 * 
 * OUPUTS :
 * 		double **S		Trace distance.
 *		  
 * ********************************************************************/
void Trace::GetTraceDist(double **S) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			S[i][j] = S_[i][j];
		}
	}	
}

/***********************************************************************
 * NAME : void Trace::GetTraceRmsm(Rmsm)
 * 
 * DESCRIPTION : Get the radial coodinate along each trace (MSM)
 *	
 * OUPUTS : 
 * 		double **Rmsm		radial coordinates (MSM)
 * 		  
 * ********************************************************************/
void Trace::GetTraceRmsm(double **Rmsm) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Rmsm[i][j] = Rmsm_[i][j];
		}
	}	
}

/***********************************************************************
 * NAME : void Trace::GetTraceRmso(Rmso)
 * 
 * DESCRIPTION : Get the radial coodinate along each trace (MSO)
 *	
 * OUPUTS : 
 * 		double **Rmso		radial coordinates (MSO)
 *		  
 * ********************************************************************/
void Trace::GetTraceRmso(double **Rmso) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Rmso[i][j] = Rmso_[i][j];
		}
	}	
}

/***********************************************************************
 * NAME : void Trace::GetTraceRnorm(Rnorm)
 * 
 * DESCRIPTION : Get the Rnorms for each trace.
 *	
 * OUPUTS :
 *		  double **Rnorm		Array of Rnorms.
 * 
 * ********************************************************************/
void Trace::GetTraceRnorm(double **Rnorm) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<nstep_[i];j++) {
			Rnorm[i][j] = Rnorm_[i][j];
		}
	}	
}

/***********************************************************************
 * NAME : void Trace::GetTraceFootprints(FP)
 * 
 * DESCRIPTION : Get the trace footprints for all traces.
 * 
 * OUPUTS :
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
 *		  
 * ********************************************************************/
void Trace::GetTraceFootprints(double **FP) {
	int i, j;
	for (i=0;i<n_;i++) {
		for (j=0;j<15;j++) {
			FP[i][j] = FP_[i][j];
		}
	}	
}

/***********************************************************************
 * NAME : void Trace::GetTraceNstep(nstep)
 * 
 * DESCRIPTION : Get the number of steps in each trace.
 *	
 * OUPUTS :
 * 		int *nstep		Number of steps in each trace.
 *		  
 * ********************************************************************/
void Trace::GetTraceNstep(int *nstep) {
	int i;
	for (i=0;i<n_;i++) {
		nstep[i] = nstep_[i];
	}
}
