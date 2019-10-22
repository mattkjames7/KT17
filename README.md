# KT17
KT17/KT14 magnetic field model for Mercury written in C++ with a Python wrapper. See Korth et al., 2015 (JGR) and Korth et al., 2017 (GRL) for more details.

## Installation

Best to install using ```pip3```:

```
pip3 install KT17 --user
```

The installation will require the following packages:

* numpy
* scipy
* matplotlib

## Usage 

### Model parameters

The `Params` keyword present in `ModelField` and `TraceField` is used to
determine the scaling factors of the KT14/KT17 magnetic field models. 

For the KT14 model, `Params` must contain three numbers: ```Params = [Rsm,t1,t2]```
where `Rsm` is the radial distance of the subsolar magnetopause and `t1` 
and `t2` are the two tail current scaling factors. By default 
```Params = [1.42,7.37,2.16]```.

For the KT17 model, `Params` must contain just two numbers: 
```Params = [Rsun,DistIndex]```, where `Rsun` is the radial distance of
Mercury from the Sun in AU and `DistIndex` is the Anderson et al., 2013 
disturbance index.

### Obtaining model field vectors

To get model field vectors for some position(s) in the magnetosphere:

```python
import KT17
Bx,By,Bz = KT17.ModelField(x,y,z,Params=[1.42,7.37,2.16])
```

where `x`, `y` and `z` are either scalars or arrays of position(s) in
the MSM coordinate system, and `Params` contains the model parameters as
described above. The returned `Bx`, `By` and `Bz` contain the magnetic
field vecotr(s) at each position.

### Tracing the magnetic field

To trace the magnetic field, use the `TraceField` object:

```python
T = KT17.TraceField(x0,y0,z0,Params=[1.42,7.37,2.16],maxlen=1000,
				initstep=0.01,maxstep=0.05,LimType=15,
				FlattenSingleTraces=True)
```

where `x0`, `y0` and `z0` are the starting position(s) of the traces in 
MSM coordinates. `maxlen` defines the maximum number of steps for the 
traces. `initstep` sets the initial step size in R<sub>M</sub>. `maxstep`
sets the maximum step size. `FlattenSingleTraces` flattens the arrays 
stored in the `TraceField` object if there is only a single trace.
`LimType` determines where the limits of the trace are (i.e. where it 
will stop tracing). The different options can be enabled by setting the 
appropriate bit of and 8-bit integer to 1, where the options are:

|   |   |  |
|:----|:---:|:---|
| MP					|	+1 	| Stop if magnetopause is reached |
| Planet				|	+2	| Stop at planetary surface |
| Dipole of Planet	|	+4	| Stop 1 Rm from centre of planetary dipole | 
| Tail Limit at 10Rm	|	+8	| Stop if trace reaches x MSM < -10 (in the magnetotail) |
| Core				|	+16	| Stop at the iron core of Mercury (R~2030km) |
| Core Dipole			|	+32 | Stop at a distance from the centre of the dipole equivalent to the size of Mercury's core |
| Box					|	+64 | Stop within a box where -6 < x < 2, -4 < y < 4, -4 < z < 4 |


Default: 1111000 (Big Endian) = 15 (stop at MP, Planetary surface in the 
south, 1Rm from the dipole in the north and 10Rm down-tail)

The `TraceField` object, `T` in the above code snippet, contains the 
following arrays:

| | |
|:--|:---|
| `x` |		x coordinate along the field trace(s) |
| `y` |			y coordinate along the field trace(s) |
| `z` |			z coordinate along the field trace(s) |
| `Bx` |			x component of the magnetic field along the trace(s) |
| `By` |			y component of the magnetic field along the trace(s) |
| `Bz` |			z component of the magnetic field along the trace(s) |	
| `nstep` |	 	number of steps along the trace(s) |
| `GlatN` |	 	Geographic latitude of the northern footprint(s) |
| `GlatS` |		Geographic latitude of the southern footprint(s) |
| `MlatN` |		Magnetic latitude of the northern footprint(s) |
| `MlatS` |		Magnetic latitude of the southern footprint(s) |
| `GlonN` |		Geographic longitude of the northern footprint(s) |
| `GlonS` |		Geographic longitude of the southern footprint(s) |
| `MlonN` |		Magnetic longitude of the northern footprint(s) |
| `MlonS` |		Magnetic longitude of the southern footprint(s) |
| `GltN` |			Geographic local time of the northern footprint(s) |
| `GltS` |			Geographic local time of the southern footprint(s) |
| `MltN` |			Magnetic local time of the northern footprint(s) |
| `MltS` |			Magnetic local time of the southern footprint(s) |
| `Lshell` |		L-shell of the field line(s) at the equator |
| `MltE` |			Magnetic local time of the equatorial footprint(s) |
| `FlLen` |		Field line length in planetary radii |
| `GlatNcore` |	 	Geographic latitude of the northern footprint(s) when traced to the outer surface of the core |
| `GlatScore` |		Geographic latitude of the southern footprint(s) when traced to the outer surface of the core |
| `MlatNcore` |		Magnetic latitude of the northern footprint(s) when traced to the outer surface of the core |
| `MlatScore` |		Magnetic latitude of the southern footprint(s) when traced to the outer surface of the core |
| `GlonNcore` |		Geographic longitude of the northern footprint(s) when traced to the outer surface of the core |
| `GlonScore` |		Geographic longitude of the southern footprint(s) when traced to the outer surface of the core |
| `MlonNcore` |		Magnetic longitude of the northern footprint(s) when traced to the outer surface of the core |
| `MlonScore` |		Magnetic longitude of the southern footprint(s) when traced to the outer surface of the core |
| `GltNcore` |			Geographic local time of the northern footprint(s) when traced to the outer surface of the core |
| `GltScore` |			Geographic local time of the southern footprint(s) when traced to the outer surface of the core |
| `MltNcore` |			Magnetic local time of the northern footprint(s) when traced to the outer surface of the core |
| `MltScore` |			Magnetic local time of the southern footprint(s) when traced to the outer surface of the core |
| `FlLencore` |		Field line length in planetary radii when traced to the outer surface of the core |
| `Rmso` |			`R = sqrt(x**2 + y**2 + (z+0.196)**2)`	 |
| `Rmsm` |			`R = sqrt(x**2 + y**2 + z**2)`	 |


### Other routines

For a quick plot of the magnetic field in the X-Z plane run:

```python
KT17.TestTrace()
```

To find out if a position in MSM coordinates is actually within the 
magnetopause or not, run the following:

```python
KT17.WithinMP(x,y,z,Rsm=1.42)
```

To get the latitude and local times of the northern and southern open-closed
field line boundaries, run:

```python
ocb = KT17.ReadOCB()
```

where `ocb` is a `numpy.recarray` object which contains the following fields:

| | |
|:--|:---|
|	`Mlt`  |	Magnetic local time of the boundary. |
|	`LctN` |	Local time in the northern hemisphere. |
|	`LctS` |	Local time in the southern hemisphere. |
|	`Mlat` |	Invariant latitude of the boundary. |
|	`LatN` |	Latitude of the northern boundary. |
|	`LatS` |	Latitude of the southern boundary. |
 
