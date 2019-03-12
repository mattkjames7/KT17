import numpy as np
from ._CFunctions import _Ckt17MultiTrace
from scipy.interpolate import InterpolatedUnivariateSpline

###### File created automatically using PopulateCtypes ######
'''
LimTypes include:
	Confining to a box: -6 < x < 2, -4 < y < 4, -4 < z < 4
	Confine to box and planet's surface
	Confine to just outside the planet
	Confine to within MP, outside planet and within 10Rm
	Default is trace to surface and to MP

New LimType options
	Each = 1 bit
	Default: 1111000 (Big Endian[I think]) = 15
	MP						+1
	Planet					+2
	Dipole of Planet		+4
	Tail Limit at 10Rm		+8
	Core					+16
	Core Dipole				+32
	Box						+64
	


'''


class TraceField(object):
	def __init__(self,x0,y0,z0,Params=[1.42,7.37,2.16],maxlen=1000,
				initstep=0.01,maxstep=0.05,LimType=15,
				FlattenSingleTraces=True,SmoothSurfaceTransition=False,
				AssumeOutsidePlanet=True):

		'''
		The TraceField object performs field line traces from initial
		positions provided using the KT17/KT14 model.
		
		Inputs:
			x0,y0,z0: Scalars or arrays containing the x-MSM (Mercury Solar
					Magnetic) coordinate of the	starting point for the trace.
			Params: 2 or 3 element array, list or tuple of model parameters.
					For the KT14 model - Params=[Rss,T1,T2] where Rss is 
					the	subsolar magnetopause distance in Rm, T1 is the 
					magnitude of the tail disk current and T2 is the 
					magnitude of the quasi-harris sheet. By default,
					Params = [1.42,7.37,2.16].
					For the KT17 model - Params=[Rsun,DistIndex], where
					Rsun is the distance of Mercury from the Sun in AU,
					DistIndex is the disturbance index (0-100) as calculated
					in Anderson et al., 2013.
			maxlen: Maximum number of steps used to trace the field.
			initstep: Initial step size for trace in Rm.
			maxstep: Maximum length of a trace step in Rm.
			LimType: This tells the C++ code where to terminate the trace.
					The different options can be enabled by setting the 
					appropriate bit of and 8-bit integer to 1, where the
					options are:
						MP						+1 	(Stop if magnetopause is reached)
						Planet					+2	(Stop at planetary surface)
						Dipole of Planet		+4	(Stop 1 Rm from centre of planetary dipole)
						Tail Limit at 10Rm		+8	(Stop if trace reaches x MSM < -10 (in the magnetotail))
						Core					+16	(Stop at the iron core of Mercury (R~2030km))
						Core Dipole				+32 (Stop at a distance from the centre of the dipole equivalent to the size of Mercury's core)
						Box						+64 (Stop within a box where -6 < x < 2, -4 < y < 4, -4 < z < 4)
					Default: 1111000 (Big Endian) = 15 (stop at MP, Planetary surface in the south, 1Rm from the dipole in the north and 10Rm down-tail)
			FlattenSingleTraces: Flattens the position and field magnitude 
					arrays created if only a single trace is performed.
			SmoothSurfaceTransition: Smooths the transition from within
					the planets mantle/crust to outside of the planet
					by inserting extra points along the trace within a small
					distance of the surface. The InPlanet attribute will 
					vary smoothly from 1.0 to 0.0 where the field trace 
					is near to the crust. Requires AssumeOutsidePlanet=False
					and traces to terminate at the planets core.
			AssumeOutsidePlanet: Assumes the the magnetosphere extends 
					all the way to the surface of the core - sets 
					InPlanet[:] = 0.0.
					
			
		Attributes:
			nstep: (n,) array with number of steps taken in each trace.
			x,y,z:	(n,maxlen) arrays containing positions in MSM of the
				field line along the trace, where n is the number of traces.
				For trace i, x[i,0:nstep[i]],y[i,0:nstep[i]],z[i,0:nstep[i]] 
				contain real values positions and x[i,nstep[i]:],
				y[i,nstep[i]:],z[i,nstep[i]:] = NAN. If n==1 and 
				FlattenSingleTraces=True, then the dimensions are reduced 
				to (maxlen,).
			Bx,By,Bz: Magnetic field along the trace in nT.
			Rmsm: Radial distance of trace points in Rm in the MSM
				coordinate system (centred upon the magnetic dipole).
			Rmso: Radial distance of each trace point in Rm in the MSO
				(Mercury Solar Orbital) coordinate system, where
				xMSO,yMSO = xMSM,yMSM and zMSO (centered on planet) = 
				zMSM+0.19.
			InPlanet: Floating point array for each trace, where InPlanet
				= 1.0 inside the planet and InPlanet = 0.0 outside. This
				is useful when modeling MHD waves at Mercury.
				
			MlatN,MlatS: Magnetic latitudes of the northern and southern fieldline footprints.
			MltN,MltS: Magnetic local times of the northern and southern fieldline footprints.
			GlatN,GlatS: Hermeographic (is that a word?) latitudes of the northern and southern fieldline footprints.
			GltN,GltS: Hermeographic local times of the northern and southern fieldline footprints.
			
			MlatNcore,MlatScore: Magnetic latitudes of the northern and southern core fieldline footprints.
			MltNcore,MltScore: Magnetic local times of the northern and southern core fieldline footprints.
			GlatNcore,GlatScore: Hermeographic (is that a word?) latitudes of the northern and southern core fieldline footprints.
			GltNcore,GltScore: Hermeographic local times of the northern and southern core fieldline footprints.
			
			MltE: Magnetic local time at the magnetic equatorial (zMSM=0) footprint.
			Lshell: Radial distance of magnetic equatorial footprint.
			FlLen: Total length of field line (NAN if not closed)
			FlLencore: Total length of field line traced to the core.
					
		'''
		

		#Convert input variables to appropriate numpy dtype:
		_x0 = np.array([x0]).astype("float64").flatten()
		_y0 = np.array([y0]).astype("float64").flatten()
		_z0 = np.array([z0]).astype("float64").flatten()
		_n = np.int32(np.size(x0))
		_maxlen = np.int32(maxlen)
		_initstep = np.float64(initstep)
		_maxstep = np.float64(maxstep)
		_nstep = np.zeros(_n,dtype='int32')
		_x = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_y = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_z = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_bx = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_by = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_bz = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_Rmsm = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_Rmso = np.zeros(_n*_maxlen,dtype='float64')+np.nan
		_FP = np.zeros(_n*20,dtype='float64')+np.nan
		_LimType = np.int32(LimType)
		_Params = np.float64(Params)
		_nParams = np.int32(_Params.size)
		#call the C++ function 
		_Ckt17MultiTrace(_x0, _y0, _z0, _n, _maxlen, _initstep, _maxstep,
						_LimType, _nParams, _Params, _nstep, _x, _y, _z,
						_bx, _by, _bz, _Rmsm, _Rmso, _FP)

		#reshape the footprints
		_FP = _FP.reshape((_n,20))
		
		fpnames = ['MlatN','MlatS','GlatN','GlatS','MltN','MltS','GltN',
					'GltS','MlatNcore','MlatScore','GlatNcore','GlatScore',
					'MltNcore','MltScore','GltNcore','GltScore','Lshell',
					'MltE','FlLen','FlLencore']

		self.n = _n
		if _n == 1 and FlattenSingleTraces:
			self.nstep = _nstep[0]
			self.x = _x
			self.y = _y
			self.z = _z
			self.Bx = _bx
			self.By = _by
			self.Bz = _bz
			self.Rmsm = _Rmsm
			self.Rmso = _Rmso
			for i in range(0,20):
				setattr(self,fpnames[i],_FP[0,i])

		else:		
			self.nstep = _nstep
			self.x = _x.reshape((_n,_maxlen))
			self.y = _y.reshape((_n,_maxlen))
			self.z = _z.reshape((_n,_maxlen))
			self.Bx = _bx.reshape((_n,_maxlen))
			self.By = _by.reshape((_n,_maxlen))
			self.Bz = _bz.reshape((_n,_maxlen))
			self.Rmsm = _Rmsm.reshape((_n,_maxlen))
			self.Rmso = _Rmso.reshape((_n,_maxlen))
			for i in range(0,20):
				setattr(self,fpnames[i],_FP[:,i])


		#make a copy of the original traces before tey are butchered by the options 
		self.nstep_full = np.copy(self.nstep)
		self.x_full = np.copy(self.x)
		self.y_full = np.copy(self.y)
		self.z_full = np.copy(self.z)
		self.Bx_full = np.copy(self.Bx)
		self.By_full = np.copy(self.By)
		self.Bz_full = np.copy(self.Bz)	
		self.Rmsm_full = np.copy(self.Rmsm)	
		self.Rmso_full = np.copy(self.Rmso)	
		
		
		bits = np.array(list('{:032b}'.format(LimType))[::-1]).astype('bool8')
		if bits[4]:
			Rcutoff = 0.832
			Ruse = self.Rmso_full
		elif bits[5] and not bits[4]:
			Rcutoff = 0.832
			Ruse = self.Rmsm_full
		elif bits[2] and not bits[1]:
			Rcutoff = 1.0
			Ruse = self.Rmsm_full
		else:
			Rcutoff = 1.0
			Ruse = self.Rmso_full
		
		if _n == 1 and FlattenSingleTraces:
			good = np.where(np.isfinite(self.x) & (Ruse >= Rcutoff))[0]
			if good.size > 0:
				if good[0] > 0:
					if np.isfinite(Ruse[good[0]-1]):
						good = np.append(good[0]-1,good)
				if good[-1] < Ruse.size-1:
					if np.isfinite(Ruse[good[-1]+1]):
						good = np.append(good,good[-1]+1)
				
				

				self.x = self.x_full[good]
				self.y = self.y_full[good]
				self.z = self.z_full[good]
				self.Bx = self.Bx_full[good]
				self.By = self.By_full[good]
				self.Bz = self.Bz_full[good]
				self.Rmsm = self.Rmsm_full[good]
				self.Rmso = self.Rmso_full[good]
				
				
				self.nstep = good.size		
				
				good = np.where(np.isfinite(self.x) & (self.Rmso >= 1.0))[0]
				if good[0] > 0:
					if np.isfinite(self.Rmso[good[0]-1]):
						good = np.append(good[0]-1,good)
				if good[-1] < self.Rmso.size-1:
					if np.isfinite(self.Rmso[good[-1]+1]):
						good = np.append(good,good[-1]+1)		
				self.InPlanet = np.ones(self.nstep,dtype='float32')
				self.InPlanet[good] = 0.0
			else:
				self.InPlanet = np.zeros(self.nstep,dtype='float32')
		else:
			for i in range(0,_n):
				good = np.where(np.isfinite(self.x[i]) & (Ruse[i] >= Rcutoff))[0]
				bad = np.where((np.isfinite(self.x[i]) == False) | (Ruse[i] < Rcutoff))[0]	
				gd = np.arange(good.size)
				bd = np.arange(bad.size) + good.size
				self.nstep[i] = good.size	
				self.x[i,gd] = self.x_full[i,good]
				self.x[i,bd] = np.nan
				self.y[i,gd] = self.y_full[i,good]
				self.y[i,bd] = np.nan
				self.z[i,gd] = self.z_full[i,good]
				self.z[i,bd] = np.nan				
				self.Bx[i,gd] = self.Bx_full[i,good]
				self.Bx[i,bd] = np.nan
				self.By[i,gd] = self.By_full[i,good]
				self.By[i,bd] = np.nan
				self.Bz[i,gd] = self.Bz_full[i,good]
				self.Bz[i,bd] = np.nan						
				self.Rmsm[i,gd] = self.Rmsm_full[i,good]
				self.Rmsm[i,bd] = np.nan						
				self.Rmso[i,gd] = self.Rmso_full[i,good]
				self.Rmso[i,bd] = np.nan						

			self.InPlanet = (self.Rmso < 1.0).astype('float32')
		
		
		if AssumeOutsidePlanet:
			self.InPlanet[:] = 0.0			
		

		if SmoothSurfaceTransition and _n == 1 and FlattenSingleTraces:
			# Here we need to add some extra points along the field line in order to  
			#smooth the transition between inside the planet and outside
			J = np.where(np.abs(self.InPlanet[1:]-self.InPlanet[:-1]) == 1)[0]
			J = J[np.where((J != 0) & (J != self.nstep-2))[0]]
			nJ = J.size
			if nJ > 0:
				#get distance along field line
				s0 = np.zeros(self.nstep,dtype='float32')
				for i in range(1,self.nstep):
					s0[i] = s0[i-1] + np.sqrt((self.x[i]-self.x[i-1])**2 + (self.y[i]-self.y[i-1])**2 + (self.z[i]-self.z[i-1])**2)
				
				nInsert = 19 # number of extra points to add between steps
				nstep_new = self.nstep + nJ*nInsert
				
				#get new distance array along field line
				s1 = np.zeros(nstep_new,dtype='float32')
				new_InPlanet = np.zeros(nstep_new,dtype='float32')
				p = 0
				for i in range(0,nJ+1):
					if i == 0:
						j0 = 0
						j1 = J[0]
					elif i == nJ:
						j0 = J[-1]+1
						j1 = self.nstep-1
					else:
						j0 = J[i-1]+1
						j1 = J[i]
					
					for j in range(j0,j1+1):
						s1[p] = s0[j]
						new_InPlanet[p] = self.InPlanet[j]
						p += 1
						
					if i < nJ:
						ds = s0[J[i]+1] - s0[J[i]]
						dip = self.InPlanet[J[i]+1] - self.InPlanet[J[i]]
						for j in range(0,nInsert):
							dj = (j+1)/(nInsert+1)
							s1[p] = s0[J[i]] + ds*dj
							new_InPlanet[p] = self.InPlanet[J[i]] + dip*(3*dj**2-2*dj**3)
							p+=1
				self.InPlanet = new_InPlanet
				self.s0 = s0
				self.s1 = s1
				Params = ['x','y','z','Bx','By','Bz','Rmsm','Rmso']
				for par in Params:
					tmp0 = np.copy(getattr(self,par))[0:self.nstep]
					
					f = InterpolatedUnivariateSpline(s0,tmp0)
					tmp1 = f(s1).astype('float32')
					
					setattr(self,par,tmp1)
				self.nstep = nstep_new
