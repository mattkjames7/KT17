import numpy as np

def PlotPlanetXZ(ax,R=1.0,Center=[0.0,0.0,0.0],zorder=10,NoShadow=False,NoonTop=True):
	'''
	Adds a planet to a plot in the X-Z plane.
	
	Inputs
	======
	ax : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	R : float
		Radius of the planet.
	Center : list
		Cartesean position of the center of the planet.
	zorder : float
		Height order of this plot object.
	NoShadow : bool 
		When False, the black bit (the shadow of the planet) is 
		plotted, when True, thent he whole planet appears white.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.
	
	'''
	a = 2*np.pi*np.arange(361,dtype='float32')/360
	x = R*np.sin(a) + Center[0]
	z = R*np.cos(a) + Center[2]
	if NoonTop:	
		ax.fill(z,x,color=[1.0,1.0,1.0],zorder=zorder)
		ax.plot(z,x,color=[0,0,0],zorder=zorder+1)
		if NoShadow == False:
			ax.fill(z[180:360],x[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)
	else:
		ax.fill(x,z,color=[1.0,1.0,1.0],zorder=zorder)
		ax.plot(x,z,color=[0,0,0],zorder=zorder+1)
		if NoShadow == False:
			ax.fill(x[180:360],z[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)


def PlotPlanetXY(fig,R=1.0,Center=[0.0,0.0,0.0],zorder=10,NoShadow=False,NoonTop=True):
	'''
	Adds a planet to a plot in the X-Y plane.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	R : float
		Radius of the planet.
	Center : list
		Cartesean position of the center of the planet.
	zorder : float
		Height order of this plot object.
	NoShadow : bool 
		When False, the black bit (the shadow of the planet) is 
		plotted, when True, thent he whole planet appears white.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	
	'''

	a = 2*np.pi*np.arange(361,dtype='float32')/360
	x = R*np.sin(a) + Center[0]
	y = R*np.cos(a) + Center[1]
	
	if NoonTop:	
		fig.fill(y,x,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(y,x,color=[0,0,0],zorder=zorder+1)
		if NoShadow == False:
			fig.fill(y[180:360],x[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)
	else:
		fig.fill(x,y,color=[1.0,1.0,1.0],zorder=zorder)
		fig.plot(x,y,color=[0,0,0],zorder=zorder+1)
		if NoShadow == False:
			fig.fill(x[180:360],y[180:360],color=[0.0,0.0,0.0],zorder=zorder+1)

def PlotPlanetYZ(fig,R=1.0,Center=[0.0,0.0,0.0],Side='day',zorder=10,NoFill=False,Color=[0.0,0.0,0.0],linestyle='-'):
	'''
	Adds a planet to a plot in the Y-Z plane.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	R : float
		Radius of the planet.
	Center : list
		Cartesean position of the center of the planet.
	zorder : float
		Height order of this plot object.
	Side : str
		'day' (planet is white) or 'night' (planet appears black).
	NoFill : bool
		If True then the centre of the planet will not be filled.
	Color : list
			RGB values for the color of the planet outline.
	linestyle : str
		Line style for the planet outline.
	'''	
	a = 2*np.pi*np.arange(361,dtype='float32')/360
	y = R*np.sin(a) + Center[1]
	z = R*np.cos(a) + Center[2]
	
	if NoFill == False:	
		if Side == 'day':
			fig.fill(y,z,color=[1.0,1.0,1.0],zorder=zorder)
		else:
			fig.fill(y,z,color=[0.0,0.0,0.0],zorder=zorder)
	fig.plot(y,z,color=Color,zorder=zorder+1,linestyle='-')

