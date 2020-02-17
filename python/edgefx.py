import numpy as np
import segyio
import os


def edgefx(bbox, hmin, **kwargs):
    """
    edgefx: build a mesh size function for seismic problems

    Usage
    -------
    >>>> fd,fh = edgefx(bbox,hmin,**kwargs) 


    Parameters 
    -------
        bbox: Bounding box, (xmin, ymin, xmax, ymax)
        hmin: Initial edgelength populating domain, (meters)
        segy: Segy file containing velocity model, (assumes velocity is m/s) 
        wl:   number of nodes per wavelength for given freq, (num. of nodes per wl.)
        freq: maximum source frequency for which to estimate wl (hertz)


    Returns 
    -------
        fd: lambda function representing a signed distance function (SDF) describing domain 
        fh: lambda function representing isotropic mesh sizes in domain


    Example
    ------


    """
    # call the desired mesh size functions
    for key, value in kwargs.items():
        if(key == "wl"):
            print('INFO: wavelength sizing function is activated...')

            freq=2 # 2 hz by default
            for key2, value2 in kwargs.items(): 
                if(key2 == "freq"): freq=value2
            print('INFO: mesh will be built for a max. source freq. of '+str(freq)+" hz...")

            for ftype, fname in kwargs.items():
                if(ftype == "segy"):
                    with segyio.open(fname,ignore_geometry=True) as f:
                        print("Succesfully read in file "+str(f)+"...")
                        for x in f:
                            
                        
                        
                        # do something with file
                        # now visualize the velocity data


    return
