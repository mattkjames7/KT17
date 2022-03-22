import numpy as np
import matplotlib.pyplot as plt
from .TraceField import TraceField
from .PlotMP import PlotKT17MercuryMP,PlotMercuryMP
from .PlotPlanet import PlotPlanetXZ

def TestTrace(degstep=5.0,Params=(0.39,50.0)):
	'''
	Produces a simple trace in X-Z plane.
	
	Inputs
	======
	degstep : float
		Latitudinal separaation between field traces in degrees.
	Params : list|tuple
		Either 2 or 3 element array containing the parameters which 
		control the model field:
			KT14: Three parameters [Rsm,t1,t2], where
				Rsm = subsolar distance of the MP in R_m
				t1,t2 = tail current scaling parameters
			KT17: Two parameters [Rsun,DistIndex], where
				Rsun = radial distance of mercury from the Sun in AU
				DistIndex = Anderson et al 2013 disturbance index.
				
	Returns
	=======
	matplotlib.pyplot.axes
	'''
	
	n = np.int32(360.0/degstep)
	a = np.arange(n)*np.pi*2.0/n
	
	x = 1.01*np.cos(a)
	y = np.zeros(n)
	z = 1.01*np.sin(a)
	
	kwargs = {}
	if len(Params) == 2:
		kwargs['Rsun'] = Params[0]
		kwargs['DistIndex'] = Params[1]
	else:
		kwargs['Rsm'] = Params[0]
		kwargs['t1'] = Params[1]
		kwargs['t2'] = Params[2]
	
	T = TraceField(x,y,z,**kwargs)
	
	fig = plt
	fig.figure()
	PlotPlanetXZ(fig,Center=[0.0,0.0,-0.196],NoonTop=False)
	
	if np.size(Params) == 2:
		PlotKT17MercuryMP(fig,Params[0],Params[1],NoonTop=False)
	else:	
		PlotMercuryMP(fig,Params[0],0.5,NoonTop=False)
	fig.plot(T.x.T,T.z.T,color=[0.0,0.0,0.0],lw=1.0)
	fig.axis([2.0,-10.0,-4.0,4.0])
	ax = fig.gca()
	ax.set_aspect(1.0)

	return ax
