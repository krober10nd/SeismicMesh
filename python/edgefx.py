import numpy as np 
import segyio as io 

def edgefx(bbox,hmin,**kwargs):
    """
    edgefx: build a mesh size function for seismic problems

    Usage
    -------
    >>>> fd,fh = edgefx(bbox,hmin,**kwargs) 


    Parameters 
    -------
        bbox: Bounding box, (xmin, ymin, xmax, ymax)
        hmin: Initial edgelength populating domain, (meters)
        segy: Segy file containing velocity model, (assumes m/s) 
        wl:   number of nodes per wavelength for given freq, (num. of nodes per wl.)
        freq: maximum source frequency for which to estimate wl (hertz)


    Returns 
    -------
        fd: lambda function representing a signed distance function (SDF) describing domain 
        fh: lambda function representing isotropic mesh sizes in domain


    Example
    ------


    """
    for key, value in kwargs.items() :
        print("key is " , key , " and the value is ", value)
    
    # call the desired mesh size functions 
    #,hmax=1e3,segy=fname,wl=0,freq=0):
    return

