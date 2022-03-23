# KT17

KT17/KT14 magnetic field model for Mercury written in C++ with a Python wrapper. See Korth et al., 2015 (JGR) and Korth et al., 2017 (GRL) for more details on the models themselves.

## Installation

Best to install using ```pip3```:

```
pip3 install KT17 --user
```

The installation will require the following packages:

* numpy
* scipy
* matplotlib

Alternatively clone this repo and build a wheel:

```bash
git clone https://github.com/mattkjames7/KT17.git
cd KT17
python3 setup.py bdist_wheel
pip3 install dist/KT17-1.0.0-py3-none-any.whl
```

NOTE: This module includes C++ code which does all of the model calculating and tracing. The package includes `libkt.so` for Linux and `libkt.dll` for Windows 10 - if either of these files can't be loaded, the code will attempt to recompile. If the compilation fails, then make sure that you have `g++`(version >= 5) and `make` installed on your system. For Windows users, install TDM-GCC.

## Usage

### Model parameters

Both `KT17.ModelField()` and `KT17.TraceField()` functions accept five different keyword arguments (`**kwargs`) which affect the model:

| Parameter   | Model | Description                                                |
| ----------- | ----- | ---------------------------------------------------------- |
| `Rsm`       | KT14  | Magnetopause standoff distance (*R<sub>M</sub>*).          |
| `t1`        | KT14  | Disk current field magnitude.                              |
| `t2`        | KT14  | Quasi-Harris sheet field magnitude.                        |
| `Rsun`      | KT17  | Distance of Mercury from the Sun (AU)                      |
| `DistIndex` | KT17  | Anderson et al., 2013 disturbance index, in the range 0-97 |

If the KT17 parameters are provided, they are used to calculate the appropriate KT14 set of parameters as defined in Korth et al., 2017.

### Obtaining model field vectors

To get model field vectors for some position(s) in the magnetosphere:

```python
import KT17

#default model
Bx,By,Bz = KT17.ModelField(x,y,z)

#Custom KT14 model
Bx,By,Bz = KT17.ModelField(x,y,z,Rsm=1.3,t1=7.0,t2=2.5)

#Custom KT17 model
Bx,By,Bz = KT17.ModelField(x,y,z,Rsun=0.4,DistIndex=25.0)

```

where `x`, `y` and `z` are either scalars or arrays of position(s) in
the MSM coordinate system. The returned `Bx`, `By` and `Bz` contain the magnetic
field vector(s) at each position in nT. For positions which are outside of the magnetopause, the magnetic field components are returned as NAN.

### Tracing the magnetic field

To trace the magnetic field, use the `TraceField` object:

```python
T = KT17.TraceField(x0,y0,z0,**kwargs)
```

where `x0`, `y0` and `z0` are the starting position(s) of the traces in 
MSM coordinates. The parameters of the trace can be provided using `**kwargs` - see the `KT17.Tracefield` dosctring for more information using `help(KT17.TraceField)` or `KT17.TraceField?` in ipython.

`TraceField` can be saved to a file and reloaded, e.g.:

```python
#save the Trace object
T.Save('Trace.bin')

#load the file
T = KT17.TraceField('Trace.bin')
```

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

|        |                                        |
|:------ |:-------------------------------------- |
| `Mlt`  | Magnetic local time of the boundary.   |
| `LctN` | Local time in the northern hemisphere. |
| `LctS` | Local time in the southern hemisphere. |
| `Mlat` | Invariant latitude of the boundary.    |
| `LatN` | Latitude of the northern boundary.     |
| `LatS` | Latitude of the southern boundary.     |

### References

Anderson, B. J., Johnson, C. L., and Korth, H. (2013), A magnetic disturbance index for Mercury's magnetic field derived from MESSENGER Magnetometer data, *Geochem. Geophys. Geosyst.*, 14, 3875– 3886, doi:[10.1002/ggge.20242](https://doi.org/10.1002/ggge.20242 "Link to external resource: 10.1002/ggge.20242")

Korth, H., Tsyganenko, N. A., Johnson, C. L., Philpott, L. C., Anderson, B. J., Al Asad, M. M., Solomon, S. C., and McNutt, R. L. (2015), Modular model for Mercury's magnetospheric magnetic field confined within the average observed magnetopause. *J. Geophys. Res. Space Physics*, 120, 4503– 4518. doi: [10.1002/2015JA021022](https://doi.org/10.1002/2015JA021022 "Link to external resource: 10.1002/2015JA021022").Korth, H., Tsyganenko, N. A., Johnson, C. L., Philpott, L. C., Anderson, B. J., Al Asad, M. M., Solomon, S. C., and McNutt, R. L. (2015), Modular model for Mercury's magnetospheric magnetic field confined within the average observed magnetopause. J. Geophys. Res. Space Physics, 120, 4503– 4518. doi: 10.1002/2015JA021022.

Korth, H., Johnson, C. L., Philpott, L., Tsyganenko, N. A., & Anderson, B. J. (2017). A dynamic model of Mercury's magnetospheric magnetic field. *Geophysical Research Letters*, 44, 10,147– 10,154. https://doi.org/10.1002/2017GL074699
