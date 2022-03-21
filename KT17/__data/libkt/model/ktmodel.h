#ifndef __KTMODEL_H__
#define __KTMODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include "kt14.h"
#include "kt17.h"



class KTModel {
	public:
		KTModel();
		~KTModel();
		KTModel(const KTModel& obj);
		
		/* functions to alter model parameters */
		void SetParams(double Rsm, double t1, double t2);
		void SetParams(double Rsun, double DistIndex);
		void GetParams(double *Rsm, double *t1, double *t2);
		void GetParams(double *Rsun, double *DistIndex);
		
		/* this function will call the model and provide a field vector*/
		void Field(	double x, double y, double z, 
					double *Bx, double *By, double *Bz);

		/* some local model parameters */
		double Rsm_, t1_, t2_, Rsun_, DistIndex_;
		
};
#endif
