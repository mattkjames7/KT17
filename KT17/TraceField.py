import numpy as np
from ._CFunctions import _CTraceField
import PyFileIO as pf
from .ct import ctBool,ctInt,ctIntPtr,ctDouble,ctDoublePtr,ctDoublePtrPtr
import matplotlib.pyplot as plt
from .PlotPlanet import PlotPlanetXY, PlotPlanetXZ, PlotPlanetYZ

class TraceField(object):
	def __init__(self,*args,**kwargs):
		'''
		Either create a new set of KT14/17 traces, or load from file.
		
		Inputs: *args
		===========
		Either a single argument:
			TraceField(fname)
		Or three:
			TraceField(x0,y0,z0)
		
		fname : str
			Name of a file which contains the contents of a TraceField
			object.
		x0 : float
			Trace starting position(s) in MSM (Rm)
		y0 : float
			Trace starting position(s) in MSM (Rm)
		z0 : float
			Trace starting position(s) in MSM (Rm)
		
		Keywords: **kwargs
		==================
		Rsm : float
			Subsolar magnetopause standoff distance (Rm).
		t1 : float
			Disk current magntitude
		t2 : float
			Quasi-harris sheet magnitude
		Rsun : float
			Distance of Mercury from the Sun (AU).
		DistIndex : float
			Anderson et al 2013 disturbance index (0.0-97.0).	
		MPStop : bool
			If True, tracing will stop at the magnetopause
		TailX : float
			Limit at which tracing will stop in the magnetotail (Rm).
		EndSurface : int
			1 - Stop at the planetary surface
			2 - Stop at the planetary core (0.832 Rm).
			3 - Stop at dipole at 1 Rm
			4 - Stop at dipole at 0.832 Rm (core radius)
			5 - Stop at northern surface, southern dipole at 
				1 Rm (virtual surface).
			6 - Stop at northern core and southern dipole at
				0.832 Rm.			
		MaxLen : int
			Maximum number of trace steps
		MaxStep : float
			Maximum step size (Rm).
		InitStep :  float
			Initial step size (Rm).
		MinStep : float
			Minimum step size (Rm).
		ErrMax : float
			Maximum trace error.
		Delta : float
			Distance between adjacent traces (Rm) used for calculating
			h_alpha
		Verbose : bool
			If True then tracing progress will be displayed
		TraceDir : int
			0 - Trace in both directions
			1 - Only trace in the direction of the field (towards north)
			-1 - Only trace in the opposite direction to the field
				(towards south)
		alpha : float
			Polarization angles to calculate h_alpha for (degrees), 
			0 degrees is toroidal, 90 is poloidal.

		Member Functions
		================
		TraceDict()
			Return a dictionary containing all of the traces.
		Save()
			Save the trace object contents to a file.
		GetTrace()
			Return a dictionary for a single trace.
		PlotXY()
			Plot Field line(s) in X-Y plane.
		PlotXZ()
			Plot Field line(s) in X-Z plane.
		PlotRhoZ()
			Plot Field line(s) in Rho-Z plane.
		PlotHalpha()
			Plot H_alpha along a field line.
		
		Attributes
		==========
		x : float
			Position along traces (Rm).
		y : float
			Position along traces (Rm).
		z : float
			Position along traces (Rm).
		Bx : float
			Field along traces (nT).
		By : float
			Field along traces (nT).
		Bz : float
			Field along traces (nT).
		nstep : int
			Number of elements in each trace.
		s : float
			Distance along each trace (Rm).
		Rmsm : float
			Radial coordinate MSM (Rm).
		Rmso : float
			Radial coordinate MSO (Rm).
		Rnorm : float
			Normalized radius.
		halpha : float
			Arrays of halphas for each trace and polarization.
		Lshell : float
			L-shell of each field line if it crosses the equator.,
		MLTe : float
			Local time at which each field line crosses the magnetic equator.
		
		Footprints:
			Each footprint has a latitude (Lat) and local time (LT), where
			Lat and LT are followed by N or S (for North or South). 
			Footprints preceded by "M" are those which sit on a dipole-
			centered surface, rather than the planet itself. Ones which end 
			with a "c" are ones which map to the core, or a surface the same
			size as the core.
			The list of footprint attributes:
			'LatN','LTN','LatS','LTS','LatNc','LTNc','LatSc','LTSc',
			'MLatN','MLTN','MLatS','MLTS','MLatNc','MLTNc','MLatSc','MLTSc'
			
		'''
		
		#check if we are loading from file, or creating new traces
		if len(args) == 1:
			#read from file or dict
			if isinstance(args[0],dict):
				#assume that the dictionary provided is a TraceField dict
				self.__dict__ = args[0]
			else:
				#load from file
				self._Load(*args)
		elif len(args) == 3:
			#new traces
			self._Trace(*args,**kwargs)
		else:
			#something's wrong
			print('TraceField was supplied with {:d} arguments...'.format(len(args)))
			print('Either use 1 string (file name), or')
			print('use 3 inputs (x,y,z)')
			return None
		
	def _Load(self,fname):
		'''
		Load the object data from a file containing the object __dict__
		'''
		self.__dict__ = pf.LoadObject(fname)
	

	def _Trace(self,*args,**kwargs):

		'''
		Run the tracing code, see __init__ docstring.
					
		'''
		

		#Convert input variables to appropriate numpy dtype:
		self.x0 = ctDoublePtr(args[0])
		self.y0 = ctDoublePtr(args[1])
		self.z0 = ctDoublePtr(args[2])
		self.n = ctInt(np.size(self.x0))

		#figure out what parameters we have and which to pass to the C code
		kt14 = np.array(['Rsm' in kwargs,'t1' in kwargs,'t2' in kwargs])
		kt17 = np.array(['Rsun' in kwargs,'DistIndex' in kwargs])
		
		if kt14.all():
			#use just kt14 parameters
			self.Params = np.zeros((self.n,3),dtype='float64')
			self.Params[:,0] = kwargs['Rsm']
			self.Params[:,1] = kwargs['t1']
			self.Params[:,2] = kwargs['t2']
		elif kt17.all():
			#use only kt17 parameters
			self.Params = np.zeros((self.n,2),dtype='float64')
			self.Params[:,0] = kwargs['Rsun']
			self.Params[:,1] = kwargs['DistIndex']
		elif kt14.any():
			#use some kt14 parameters + defaults
			self.Params = np.zeros((self.n,3),dtype='float64')
			self.Params[:,0] = kwargs.get('Rsm',1.42)
			self.Params[:,1] = kwargs.get('t1',7.37)
			self.Params[:,2] = kwargs.get('t2',2.16)		
		elif kt17.any():
			#use some kt17 parameters + defaults
			self.Params = np.zeros((self.n,2),dtype='float64')
			self.Params[:,0] = kwargs.get('Rsun',0.427)
			self.Params[:,1] = kwargs.get('DistIndex',50.0)		
		else:
			#use defaults
			self.Params = np.zeros((self.n,3),dtype='float64')
			self.Params[:,0] = 1.42
			self.Params[:,1] = 7.37
			self.Params[:,2] = 2.16		

		_,_nP = self.Params.shape
		
		#switch to three arrays
		_P0 = ctDoublePtr(self.Params[:,0])
		_P1 = ctDoublePtr(self.Params[:,1])
		if _nP == 2:
			_P2 = np.zeros(self.n,dtype='float64')
		else:
			_P2 = ctDoublePtr(self.Params[:,2])

		#some config options
		self.BoundMP = ctBool(kwargs.get('MPStop',True))
		self.BoundTail = ctDouble(kwargs.get('TailX',-10.0))
		self.BoundSurface = ctInt(kwargs.get('EndSurface',6))
		
		self.MaxLen = ctInt(kwargs.get('MaxLen',1000))
		self.MaxStep = ctDouble(kwargs.get('MaxStep',0.05))
		self.InitStep = ctDouble(kwargs.get('InitStep',0.01))
		self.MinStep = ctDouble(kwargs.get('MinStep',0.001))
		self.ErrMax = ctDouble(kwargs.get('ErrMax',0.0001))
		self.Delta = ctDouble(kwargs.get('Delta',0.05))
		self.Verbose = ctBool(kwargs.get('Verbose',False))
		self.TraceDir = ctInt(kwargs.get('TraceDir',0))
		
		#alpha
		self.alpha = ctDoublePtr(kwargs.get('alpha',[]))
		self.nalpha = np.int32(self.alpha.size)
		
		#some output arrays
		self.x = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.y = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.z = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Bx = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.By = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Bz = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan

		self.nstep = np.zeros(self.n,dtype="int32")

		self.s = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Rmsm = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Rmso = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.Rnorm = np.zeros((self.n,self.MaxLen),dtype="float64") + np.nan
		self.halpha = np.zeros((self.n*self.MaxLen*self.nalpha,),dtype="float64") + np.nan #hopefully this will be reshaped to (n,nalpha,MaxLen)
		self.FP = np.zeros((self.n,18),dtype="float64")

		_x = ctDoublePtrPtr(self.x)
		_y = ctDoublePtrPtr(self.y)
		_z = ctDoublePtrPtr(self.z)

		_Bx = ctDoublePtrPtr(self.Bx)
		_By = ctDoublePtrPtr(self.By)
		_Bz = ctDoublePtrPtr(self.Bz)
	
		
		_s = ctDoublePtrPtr(self.s)
		_Rmsm = ctDoublePtrPtr(self.Rmsm)
		_Rmso = ctDoublePtrPtr(self.Rmso)
		_Rnorm = ctDoublePtrPtr(self.Rnorm)		
		_FP = ctDoublePtrPtr(self.FP)
	
		#call the C++ function 
		_CTraceField(	self.n,self.x0,self.y0,self.z0,
						_nP,_P0,_P1,_P2,
						self.BoundMP,self.BoundTail,self.BoundSurface,
						self.MaxLen,self.MaxStep,self.InitStep,
						self.MinStep,self.ErrMax,self.Delta,
						self.Verbose,self.TraceDir,
						self.nstep,
						_x,_y,_z,
						_Bx,_By,_Bz,
						_Rmsm,_Rmso,_s,_Rnorm,_FP,
						self.nalpha,self.alpha,self.halpha)

		
		
		
		#add footprints as attributes to this object		
		fpnames = [	'LatN','LTN','LatS','LTS',
					'LatNc','LTNc','LatSc','LTSc',
					'MLatN','MLTN','MLatS','MLTS',
					'MLatNc','MLTNc','MLatSc','MLTSc',
					'Lshell','MLTe']

		for i in range(0,len(fpnames)):
			setattr(self,fpnames[i],self.FP[:,i])

		#reshape halpha
		if self.nalpha > 0 :
			self.halpha = self.halpha.reshape(self.n,self.nalpha,self.MaxLen)
	
	def TraceDict(self,RemoveNAN=True):
		'''
		Return a dictionary with all of the outputs of the field trace.
		
		Inputs
		======
		RemoveNAN : bool
			If True then arrays will be shortened by removing nans.
			
		Returns
		=======
		out : dict
			Contains the field traces coordinates, field components etc.
		
		'''
		#we could save a fair bit of space by removing NANs - this will
		#mean that simple 2D arrays will become arrays of objects
		if RemoveNAN:
			ptrs = ['x','y','z','Bx','By','Bz','s','Rmsm','Rmso','Rnorm']
			out = {}
			keys = list(self.__dict__.keys())
			for k in keys:
				if k in ptrs:
					#fix these
					tmp = np.zeros(self.n,dtype='object')
					for i in range(0,self.n):
						tmp[i] = self.__dict__[k][i,:self.nstep[i]]
					out[k] = tmp
				elif k == 'halpha':
					#3D
					tmp = np.zeros(self.halpha.shape[:2],dtype='object')
					for i in range(0,self.n):
						for j in range(0,self.nalpha):
							tmp[i,j] = self.halpha[i,j,:self.nstep[i]]
					out[k] = tmp
				else:
					out[k] = self.__dict__[k]
		else:
			out = self.__dict__
		return out


	def Save(self,fname,RemoveNAN=True):
		'''
		Save the data in this object to file.
		
		Inputs
		======
		fname : str
			Path to the file where this trace will be save on disk.
		RemoveNAN : bool
			If True then arrays will be shortened by removing nans.
			
		'''
		out = self.TraceDict(RemoveNAN)
		
		print('Saving file: {:s}'.format(fname))
		
		pf.SaveObject(out,fname)


	def GetTrace(self,i):
		'''
		Return a Trace
		
		Inputs
		======
		i : int
			Index of the trace to be returned.
		
		Returns
		=======
		out : dict 
			Dictionary containing the following fields:
			x : float
				x-coordinate (Rm)
			y : float
				y-coordinate (Rm)
			z : float
				z-coordinate (Rm)
			Bx : float
				x-component of the magnetic field (nT)
			By : float
				y-component of the magnetic field (nT)
			Bz : float
				z-component of the magnetic field (nT)
			Rmsm : float
				radial distance MSM (Rm)
			Rmso : float
				radial distance MSO (Rm)
			Rnorm : float
				Normalised radial distance (Rnorm = 1.0 at Rmax)
			s : float
				Distance along the field line trace (Rm)
			h : float
				H_alpha array.
		
		'''
		out = {}
		out['x'] = self.x[i][:self.nstep[i]]
		out['y'] = self.y[i][:self.nstep[i]]
		out['z'] = self.z[i][:self.nstep[i]]
		out['Bx'] = self.Bx[i][:self.nstep[i]]
		out['By'] = self.By[i][:self.nstep[i]]
		out['Bz'] = self.Bz[i][:self.nstep[i]]

		out['Rmsm'] = self.Rmsm[i][:self.nstep[i]]
		out['Rmso'] = self.Rmso[i][:self.nstep[i]]
		out['Rnorm'] = self.Rnorm[i][:self.nstep[i]]
		out['s'] = self.s[i][:self.nstep[i]]
		if self.nalpha > 0:
			out['h'] = self.halpha[i,:][:self.nstep[i]]
		else:
			out['h'] = None
			
		return out



	def PlotXZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the X-Z plane
		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		
		x = self.x[ind]
		z = self.z[ind]
	
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(x[i],z[i],color=color)
			mx = np.nanmax([mx,np.abs(x[i]).max(),np.abs(z[i]).max()])
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{MSM}$ (R$_M$)')
		ax.set_xlabel('$x_{MSM}$ (R$_M$)')

		mx = 1.1*mx	
		ax.set_xlim(-mx,mx)
		ax.set_ylim(-mx,mx)
		
		PlotPlanetXZ(ax,Center=[0.0,0.0,-0.196],NoonTop=False)
		ax.set_aspect(1.0)

		return ax
	
	def PlotXY(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the X-Y plane
		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines		
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		x = self.x[ind]
		y = self.y[ind]

			
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(y[i],x[i],color=color)
			mx = np.nanmax([mx,np.abs(x[i]).max(),np.abs(y[i]).max()])

		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		yl = ax.get_xlim()
		ax.set_xlim(yl[::-1])
		
		ax.set_xlabel('$y_{MSM}$ (R$_M$)')
		ax.set_ylabel('$x_{MSM}$ (R$_M$)')


		mx = 1.1*mx
		ax.set_xlim(mx,-mx)
		ax.set_ylim(-mx,mx)
		
		PlotPlanetXY(ax)
		ax.set_aspect(1.0)
		return ax
	
	def PlotRhoZ(self,ind='all',fig=None,maps=[1,1,0,0],label=None,color='black'):
		'''
		Plot field lines in the rho-Z plane

		
		Inputs
		======
		ind : int|str
			Index of trace to plot. Can be scalar or an array. If set 
			ind='all' then all traces will be plotted.
		fig : None|pyplot|pyplot.Axes instance
			None - new figure will be created
			pyplot - new subplot will be created on existing figure
			pyplot.Axes - existing subplot will be used
		maps : list
			4-element array-like to determine the subplot position,
			ignored when fig=pyplot.Axes.
			maps = [xmaps,ymaps,xmap,ymap]
			xmaps - number of subplots in x-direction
			ymaps - number of subplots in y-direction
			xmap - x position of this subplot
			ymap - y position of this subplot
		label : None|str
			Add label to traces.
		color : str|array-like
			Colour to plot the field lines		
		'''
		
		if ind == 'all':
			ind = np.arange(self.n)
		elif np.size(ind) == 1:
			ind = np.array([ind]).flatten()
		else:
			ind = np.array(ind)
			
		
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig
		
		x = self.x[ind]
		y = self.y[ind]
		z = self.z[ind]

		
		r = np.array([np.sqrt(x[i]**2 + y[i]**2) for i in range(0,x.shape[0])],dtype='object')
		mx = 1.5
		for i in range(0,ind.size):
			ln = ax.plot(r[i],z[i],color=color)
			mx = np.nanmax([mx,np.abs(r[i]).max(),np.abs(z[i]).max()])
		if not label is None:
			hs,ls = GetLegendHandLab(ax)
			hs.append(ln[0])
			ls.append(label)
			ax.legend(hs,ls)
		
		ax.set_ylabel('$z_{MSM}$ (R$_M$)')
		ax.set_xlabel(r'$\rho_{MSM}$ (R$_M$)')


		mx = 1.1*mx		
		ax.set_xlim(-mx,mx)
		ax.set_ylim(-mx,mx)
		
		PlotPlanetXZ(ax,NoShadow=True,Center=[0.0,0.0,-0.196],NoonTop=False)
		ax.set_aspect(1.0)
		return ax
	
	
	def PlotHalpha(self,TI='all',AI='all',fig=None,maps=[1,1,0,0]):
		'''
		Plot h_alpha (see Singer et al 1982) for a field line.
		
		Inputs
		======
		TI : int|str
			Index of trace to plot. TI='all' will plot for all traces.
		AI : int|str
			Index of alpha angle to plot for. AI will plot all alphas.
		fig : None|matplotlib.pyplot|matplotlib.pyplot.Axes
			None - a new figure will be created with new axes
			matplotlib.pyplot - existing figure, new axes
			matplotlib.pyplot.Axes - existing axes instance to be used
				(maps ignored in the case).
		maps : list|tuple|numpy.ndarray
			Four element array-like, denoting subplot position,
			e.g. [xmaps,ymaps,xmap,ymap]
				xmaps : number of subplots in x-direction
				ymaps : number of subplots in y-direction
				xmap : position index (0 is left)
				ymap : position index (0 is top)
		
		
		'''
		if AI == 'all':
			AI = np.arange(self.nalpha)
		
		if np.size(AI) == 1:
			AI = np.array([AI]).flatten()
			
		if TI == 'all':
			TI = np.arange(self.n)
		
		if np.size(TI) == 1:
			TI = np.array([TI]).flatten()
			
		if fig is None:
			fig = plt
			fig.figure()
		if hasattr(fig,'Axes'):	
			ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
		else:
			ax = fig		
			
		for t in TI:
			for a in AI:
				ax.plot(self.s[t],self.halpha[t,a],label=r'Trace {:d} $\alpha=${:5.1f}'.format(t,self.alpha[a]))

		ax.legend()
		ax.set_xlabel(r'$s$ (R$_M$)')
		ax.set_ylabel(r'$h_{\alpha}$')

		return ax
