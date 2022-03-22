import setuptools
from setuptools.command.install import install
import os


with open("README.md", "r") as fh:
    long_description = fh.read()

def getversion():
	'''
	read the version string from __init__
	
	'''
	#get the init file path
	thispath = os.path.abspath(os.path.dirname(__file__))+'/'
	initfile = thispath + 'KT17/__init__.py'
	
	#read the file in
	f = open(initfile,'r')
	lines = f.readlines()
	f.close()
	
	#search for the version
	version = 'unknown'
	for l in lines:
		if '__version__' in l:
			s = l.split('=')
			version = s[-1].strip().strip('"').strip("'")
			break
	return version
	
version = getversion()

setuptools.setup(
    name="KT17",
    version=version,
    author="Matthew Knight James",
    author_email="mattkjames7@gmail.com",
    description="KT17/KT14 magnetic field model for Mercury written in C++ with a Python wrapper. See Korth et al., 2015 (JGR) and Korth et al., 2017 (GRL) for more details.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mattkjames7/KT17",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    install_requires=[
		'numpy',
		'scipy',
		'matplotlib',
		'PyFileIO',
	],
	include_package_data=True,
)



