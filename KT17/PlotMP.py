import numpy as np

def PlotShueMP(fig,Rsm,Alpha,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	

	theta = np.pi*(np.arange(359)-179)/180.0
	r = Rsm*(2.0/(1+np.cos(theta)))**Alpha
	x = r*np.cos(theta)
	z = r*np.sin(theta)
	
	if NoonTop:
		fig.plot(z,x,color=color,lw=lw,zorder=zorder,label='Magnetopause')
	else:
		fig.plot(x,z,color=color,lw=lw,zorder=zorder,label='Magnetopause')



def PlotShueMPYZ(fig,Rsm,Alpha,color=[0.0,0.0,1.0],lw=2,zorder=10):

	theta = np.pi*(np.arange(359)-179)/180.0
	r = Rsm*(2.0)**Alpha
	y = r*np.cos(theta)
	z = r*np.sin(theta)
	
	fig.plot(y,z,color=col,lw=lw,zorder=zorder,label='Magnetopause')


def PlotMercuryMP(fig,Rsm=1.42,Alpha=0.5,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	PlotShueMP(fig,Rsm,Alpha,color,lw,zorder,NoonTop)

def PlotMercuryMPYZ(fig,Rsm=1.42,Alpha=0.5,color=[0.0,0.0,1.0],lw=2,zorder=10):
	PlotShueMPYZ(fig,Rsm,Alpha,color,lw,zorder)

def PlotKT17MercuryMP(fig,Rsun=0.39,DistInd=50.0,color=[0.0,0.0,1.0],lw=2,zorder=10,NoonTop=True):
	f=2.06873-0.00279*DistInd
	Rsm=f*Rsun**(1.0/3.0) 	
	
	PlotShueMP(fig,Rsm,0.5,color,lw,zorder,NoonTop)

def PlotKT17MercuryMPYZ(fig,Rsun=0.39,DistInd=50.0,color=[0.0,0.0,1.0],lw=2,zorder=10):
	f=2.06873-0.00279*DistInd
	Rsm=f*Rsun**(1.0/3.0) 	
	PlotShueMPYZ(fig,Rsm,0.5,color,lw,zorder)
