#include "qhsheet.h"

const double d0 = 0.09;	//half-width of current sheet in Z at inner edge of tail current [Rm]
const double deltax = 1.0;//expansion magnitudes of tail current sheet in x direction
const double deltay = 0.1;//expansion magnitudes of tail current sheet in y direction
const double SxQH = 1.5;	//e-folding distance for the sunward expansion of the harris sheet
const double SyQH = 9.0;	//scale distance for the flankward expansion of the harris sheet
const double zshift = 3.5;//location of image sheets for harris sheet

/* Shielding parameters for the Quasi-Harris sheet */
const double TailQHA[36] = {  -91.67686636,   -87.31240824,   251.8848107 ,    95.65629983,   -80.968107  ,   198.1447476 ,
							 -283.1968987 ,  -269.1514899 ,   504.632231  ,   166.027215  ,  -214.9025413 ,   623.7920115 ,
							  -35.99544615,  -322.864469  ,   345.710579  ,   928.8553184 ,   810.129509  ,    19.62627762,
							  -12.70326428,   490.4662048 ,  -814.0985363 , -1781.184984  , -1371.261326  ,    60.3136479 ,
							  116.630551  ,  -178.3347065 ,   604.0308838 ,  1155.151174  ,   770.3896601 ,  -202.8545948 ,
							  298.6337705 ,   304.7964641 ,    33.70850254,   393.6080147 ,   308.1194271 ,  -660.1691658 };
const double TailQHP[6] = { 1.67762971,  1.29222658,  0.31162534, -0.43926691,  0.75780748,  1.49777952};


/***********************************************************************
 * NAME : 		double QHThickness(x,y)
 * 
 * DESCRIPTION : 	Calculates the thickness of the quasi-harris sheet.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 *
 * RETURNS :
 * 			double d	Thickness of the current sheet in Rm.
 *
 * 
 * ********************************************************************/
double QHThickness(double x, double y) {
	double y_SyQH = y/SyQH;
	return d0 + deltax*exp(x/SxQH) + deltay*y_SyQH*y_SyQH;
}

/***********************************************************************
 * NAME : 		void QHSheetField(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's quasi-harris sheet field at a 
 * 					given postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet field (in nT).
 * 			double *By	y-component of Mercury's QH sheet field (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet field (in nT).
 *
 * 
 * ********************************************************************/
void QHSheetField(	double x, double y, double z, 
					double *Bx, double *By, double *Bz) {
	double d = QHThickness(x,y);
	double dddx = deltax/SxQH*exp(x/SxQH);
	double zpzi = z+zshift;
    double zmzi = z-zshift;
    double tanhzd = tanh(z/d);
    double tanhzmd = tanh(zmzi/d);
    double tanhzpd = tanh(zpzi/d);
	Bx[0] = (tanhzd-0.5*(tanhzmd + tanhzpd))/d;
	By[0] = 0.0;
	Bz[0] = (z*tanhzd-0.5*(zmzi*tanhzmd + zpzi*tanhzpd))*dddx/(d*d);

}

/***********************************************************************
 * NAME : 		void QHSheetShield(x,y,z,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates the shielding field to contain the 
 * 					quasi-harris sheet field within the magnetopause.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet shield (in nT).
 * 			double *By	y-component of Mercury's QH sheet shield (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet shield (in nT).
 *
 * 
 * ********************************************************************/
void QHSheetShield(	double x, double y, double z, 
					double *Bx, double *By, double *Bz) {
	
	/* this calculates the shielding field using equation 12 of 
	 * Korth et al 2015 */
	ShieldField(x,y,z,6,TailQHA,TailQHP,Bx,By,Bz);
}	


/***********************************************************************
 * NAME : 		void QHSheet(x,y,z,t2,Bx,By,Bz)
 * 
 * DESCRIPTION : 	Calculates Mercury's quasi-harris sheet field at a 
 * 					given postion including the magnetopause shielding
 * 					contribution.
 * 
 * INPUTS : 
 * 			double x	x-MSM coordinate (in Rm).
 * 			double y	y-MSM coordinate (in Rm).
 * 			double z	z-MSM coordinate (in Rm).
 * 			double t2	Scale factor for the field contribution to the
 * 						model.
 *
 * OUTPUTS :
 * 			double *Bx	x-component of Mercury's QH sheet field (in nT).
 * 			double *By	y-component of Mercury's QH sheet field (in nT).
 * 			double *Bz	z-component of Mercury's QH sheet field (in nT).
 *
 * 
 * ********************************************************************/
void QHSheet(	double x, double y, double z, double t2,
				double *Bx, double *By, double *Bz) {
	
	double Fx, Fy, Fz, Sx, Sy, Sz;
	
	/* get disk field and its corresponding shielding field, then add */
	QHSheetField(x,y,z,&Fx,&Fy,&Fz);
	QHSheetShield(x,y,z,&Sx,&Sy,&Sz);
	Bx[0] = t2*(Fx - Sx);
	By[0] = t2*(Fy - Sy);
	Bz[0] = t2*(Fz - Sz);
	
}
