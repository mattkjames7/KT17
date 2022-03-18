#ifndef __KTMODEL_H__
#define __KTMODEL_H__
#include <stdio.h>
#include <stdlib.h>
#include "kt14.h"
#include "kt17.h"
#endif


class KTModel {
	public:
		KTModel();
		~KTModel();
		KTModel(const KTModel& obj);
		
		void SetParams(double Rsm, double t1, double t2);
		void SetParams(double Rsun, double DistIndex);
		void GetParams(double *Rsm, double *t1, double *t2);
		void GetParams(double *Rsun, double *DistIndex);
	private:
	
}
