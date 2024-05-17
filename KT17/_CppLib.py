import numpy as np
import os
import subprocess
import ctypes
import platform
import fnmatch

def getLibFilename(isShort=False):
    """
    Return library filename string
    
    Inputs
    ======
    isShort : bool 
        If False return filename with full path, if True return only filename
        default - False
    
    Returns
    =======
    libFilename    : str
        Filename of the source library

    """
    if(isShort):
        libFilename = "libkt."
    else:
        libFilename = os.path.dirname(__file__) + "/__data/libkt/lib/libkt."

    systype = platform.system()
    if systype == 'Linux':
        extension = "so"
    elif systype == 'Windows':
        extension = "dll"
    elif systype == 'Darwin':
        extension = 'dylib'
    else:
        raise Exception("The Operating System is not supported")
    
    return libFilename + extension


def checkLibExists():
    """Check if library file exist, and start compilation script if not."""
    if not os.path.isfile(getLibFilename()):
        print(getLibFilename(isShort=True)+" not found, try reinstalling")
        raise SystemError


def getWindowsSearchPaths():
    '''Scan the directories within PATH and look for std C++ libs'''
    paths = os.getenv('PATH')
    paths = paths.split(';')

    pattern = 'libstdc++*.dll'

    out = []
    for p in paths:
        if os.path.isdir(p):
            files = os.listdir(p)
            mch = any(fnmatch.fnmatch(f,pattern) for f in files)
            if mch:
                out.append(p)
    
    return out



def addWindowsSearchPaths():

    paths = getWindowsSearchPaths()
    for p in paths:
        if os.path.isdir(p):
            os.add_dll_directory(p)

    


