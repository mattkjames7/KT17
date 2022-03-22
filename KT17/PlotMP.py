import numpy as np

def PlotShueMP(fig,Rsm,Alpha,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	'''
	Plots the Shue et al 1997 magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsm : float
		Subsolar distance of the magnetopause.
	Alpha : float
	 	Flaring parameter of the magnetopause.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	
	
	'''

	theta = np.pi*(np.arange(359)-179)/180.0
	r = Rsm*(2.0/(1+np.cos(theta)))**Alpha
	x = r*np.cos(theta)
	z = r*np.sin(theta)
	
	if NoonTop:
		fig.plot(z,x,color=color,lw=lw,zorder=zorder,label='Magnetopause')
	else:
		fig.plot(x,z,color=color,lw=lw,zorder=zorder,label='Magnetopause')



def PlotShueMPYZ(fig,Rsm,Alpha,color=[0.0,0.0,1.0],lw=2,zorder=10):
	'''
	Plots the Shue et al 1997 magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsm : float
		Subsolar distance of the magnetopause.
	Alpha : float
	 	Flaring parameter of the magnetopause.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	

	
	'''
	theta = np.pi*(np.arange(359)-179)/180.0
	r = Rsm*(2.0)**Alpha
	y = r*np.cos(theta)
	z = r*np.sin(theta)
	
	fig.plot(y,z,color=col,lw=lw,zorder=zorder,label='Magnetopause')


def PlotMercuryMP(fig,Rsm=1.42,Alpha=0.5,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	'''
	Plots Mercury's magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsm : float
		Subsolar distance of the magnetopause.
	Alpha : float
	 	Flaring parameter of the magnetopause.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.		
	
	'''
	PlotShueMP(fig,Rsm,Alpha,color,lw,zorder,NoonTop)

def PlotMercuryMPYZ(fig,Rsm=1.42,Alpha=0.5,color=[0.0,0.0,1.0],lw=2,zorder=10):
	'''
	Plots Mercury's magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsm : float
		Subsolar distance of the magnetopause.
	Alpha : float
	 	Flaring parameter of the magnetopause.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	
	'''
	PlotShueMPYZ(fig,Rsm,Alpha,color,lw,zorder)

def PlotKT17MercuryMP(fig,Rsun=0.39,DistInd=50.0,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	'''
	Plots Mercury's magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsun : float
		radial distance of mercury from the Sun in AU
	DistIndex : float
		Anderson et al 2013 disturbance index.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	
	
	'''

	f=2.06873-0.00279*DistInd

	Rsm=f*Rsun**(1.0/3.0) 	
	
	PlotShueMP(fig,Rsm,0.5,color,lw,zorder,NoonTop)

def PlotKT17MercuryMPYZ(fig,Rsun=0.39,DistInd=50.0,color=[0.0,0.0,1.0],lw=2,zorder=10):
	'''
	Plots Mercury's magnetopause.
	
	Inputs
	======
	fig : obj
		Either a matplotlib.pyplot or matplotlib.pyplot.axes instance
	Rsun : float
		radial distance of mercury from the Sun in AU
	DistIndex : float
		Anderson et al 2013 disturbance index.
	color : list
	 	RGB color to plot the magnetopause in.
	lw : float
		Width of the magnetopause line.
	zorder : float
		Height order of this plot object.
	NoonTop : bool
		When True, the plot is oriented such that noon is at the top
		otherwise, noon is to the left of the plot.	

	'''


	f=2.06873-0.00279*DistInd
	Rsm=f*Rsun**(1.0/3.0) 	
	PlotShueMPYZ(fig,Rsm,0.5,color,lw,zorder)
