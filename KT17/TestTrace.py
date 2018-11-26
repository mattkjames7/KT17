import numpy as np
import matplotlib.pyplot as plt
from .TraceField import TraceField
from .PlotMP import PlotKT17MercuryMP,PlotMercuryMP
from .PlotPlanet import PlotPlanetXZ

def TestTrace(degstep=5.0,Params=(0.39,50.0)):
	n = np.int32(360.0/degstep)
	a = np.arange(n)*np.pi*2.0/n
	
	x = 1.01*np.cos(a)
	y = np.zeros(n)
	z = 1.01*np.sin(a)
	
	T = TraceField(x,y,z,Params=Params)
	
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
