#include "ktmodel.h"

/***********************************************************************
 * NAME : KTModel()
 * 
 * DESCRIPTION : Constructor for the KTModel object.
 * 
 * ********************************************************************/
KTModel::KTModel() {
	
	/* initialize the model parameters using default KT14 values */
	SetParams(1.42,7.64,2.06);
	
}

/***********************************************************************
 * NAME : ~KTModel()	
 * 
 * DESCRIPTION : Destructor for the KTModel object.
 *
 * ********************************************************************/
KTModel::~KTModel() {
	/* nothing really to do here */
}

/***********************************************************************
 * NAME : KTModel(obj)
 * 
 * DESCRIPTION : Copy constructor for the KTModel object.
 * 
 * INPUTS : 
 * 		const KTModel &obj	Object to be copied.
 *
 * ********************************************************************/
KTModel::KTModel(const KTModel &obj) {
	/* copy parameters */
	Rsm_ = obj.Rsm_;
	t1_ = obj.t1_;
	t2_ = obj.t2_;
	Rsun_ = obj.Rsun_;
	DistIndex_ = obj.DistIndex_;
	
}

/***********************************************************************
 * NAME : KTModel::SetParams(Rsm,t1,t2)	
 * 
 * DESCRIPTION : Set the model parameters using KT14 parameters.
 * 
 * INPUTS : 
 * 		double Rsm		Subsolar magnetopause standoff distance (Rm).
 * 		double t1		Magnitude of the disk current.
 * 		double t2		Magnitude of the quasi-harris sheet.
 * 
 * ********************************************************************/
void KTModel::SetParams(double Rsm, double t1, double t2) {
	
	/* set the KT14 parameters */
	Rsm_ = Rsm;
	t1_ = t1;
	t2_ = t2;
	
	/* then the KT17 parameters */
	KT17Params(Rsm,t1,t2,&Rsun_,&DistIndex_);
	
}

/***********************************************************************
 * NAME : KTModel::SetParams(Rsun,DistIndex)
 * 
 * DESCRIPTION : Set the model parameters up using the KT17 parameters,
 * 		which will be converted to KT14 ones.
 * 
 * INPUTS : 
 * 		double Rsun			Distance from Mercury to the Sun (AU).
 * 		double DistIndex	Anderson et al 2013 disturbance index (0-97)
 * 
 * ********************************************************************/
void KTModel::SetParams(double Rsun, double DistIndex) {

	/* set the KT17 Parameters */
	Rsun_ = Rsun;
	DistIndex_ = DistIndex;
	
	/* calculate the KT14 parameters */
	KT14Params(Rsun,DistIndex,&Rsm_,&t1_,&t2_);
	
}

/***********************************************************************
 * NAME : KTModel::GetParams(Rsm,t1,t2)
 * 
 * DESCRIPTION : Get the KT14 model parameters.
 *
 * OUTPUTS :
 * 		double *Rsm		Subsolar magnetopause standoff distance (Rm).
 * 		double *t1		Magnitude of the disk current.
 * 		double *t2		Magnitude of the quasi-harris sheet.
 * 
 * ********************************************************************/
void KTModel::GetParams(double *Rsm, double *t1, double *t2) {
	
	Rsm[0] = Rsm_;
	t1[0] = t1_;
	t2[0] = t2_;

}

/***********************************************************************
 * NAME : KTModel::GetParams(Rsun,DistIndex)	
 * 
 * DESCRIPTION : Get the KT17 model parameters.
 * 
 * OUTPUTS :
 * 		double Rsun			Distance from Mercury to the Sun (AU).
 * 		double DistIndex	Anderson et al 2013 disturbance index (0-97)
 * 
 * ********************************************************************/
void KTModel::GetParams(double *Rsun, double *DistIndex) {
	
	Rsun[0] = Rsun_;
	DistIndex[0] = DistIndex_;

}

/***********************************************************************
 * NAME : KTModel::Field(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : Get the model field vector at a position in MSM
 * 		coordinates.
 * 
 * INPUTS : 
 * 		double x	x-MSM coordinate.
 * 		double y	y-MSM coordinate.
 * 		double z	z-MSM coordinate.
 *
 * OUTPUTS :
 * 		double *Bx	x-component of the magnetic field.
 * 		double *By	y-component of the magnetic field.
 * 		double *Bz	z-component of the magnetic field.
 * 
 * ********************************************************************/
void KTModel::Field(double x, double y, double z,
					double *Bx, double *By, double *Bz) {
						
	/* call the KT14 version of the code */
	KT14(x,y,z,Rsm_,t1_,t2_,Bx,By,Bz);
	
}
