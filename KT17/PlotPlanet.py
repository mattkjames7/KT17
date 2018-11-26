import numpy as np

def PlotPlanetXZ(fig,R=1.0,Center=[0.0,0.0,0.0],zorder=10,NoBlack=False,NoonTop=True):
	

	a = 2*np.pi*np.arange(361,dtype='float32')/360
	x = R*np.sin(a) + Center[0]
	z = R*np.cos(a) + Center[2]
	
	if NoonTop:	
		fig.fill(z,x,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(z,x,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(z[180:360],x[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)
	else:
		fig.fill(x,z,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(x,z,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(x[180:360],z[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)


def PlotPlanetXY(fig,R=1.0,Center=[0.0,0.0,0.0],zorder=10,NoBlack=False,NoonTop=True):
	a = 2*np.pi*np.arange(361,dtype='float32')/360
	x = R*np.sin(a) + Center[0]
	y = R*np.cos(a) + Center[1]
	
	if NoonTop:	
		fig.fill(y,x,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(y,x,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(y[180:360],x[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)
	else:
		fig.fill(x,y,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(x,y,color=[0,0,0],zorder=zorder+1)
		if NoBlack == False:
			fig.fill(x[180:360],y[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)

def PlotPlanetYZ(fig,R=1.0,Center=[0.0,0.0,0.0],Side='day',zorder=10,NoFill=False,Color=[0.0,0.0,0.0],linestyle='-'):
	
	a = 2*np.pi*np.arange(361,dtype='float32')/360
	y = R*np.sin(a) + Center[1]
	z = R*np.cos(a) + Center[2]
	
	if NoFill == False:	
		if Side == 'day':
			fig.fill(y,z,color=[1.0,1.0,1.0],zorder=zorder)
		else:
			fig.fill(y,z,color=[0.0,0.0,0.0],zorder=zorder)
	fig.plot(y,z,color=Color,zorder=zorder+1,linestyle='-')

