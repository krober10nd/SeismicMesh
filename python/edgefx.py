import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import matplotlib
import os,sys

import utils 


def edgefx(bbox, hmin, segy,  **kwargs):
    """
    edgefx: build a mesh size function for seismic problems

    Usage
    -------
    >>>> fh = edgefx(bbox,hmin,segy,**kwargs)


    Parameters
    -------
        bbox: Bounding box, (xmin, xmax, ymin, ymax)
        hmin: Minimum edgelength populating domain, (meters)
        segy: Segy file containing velocity model, (assumes velocity is in KM PER SECOND)
                                **kwargs 
        wl:   number of nodes per wavelength for given max. freq, (num. of nodes per wl.)
        freq: maximum source frequency for which to estimate wl (hertz)
        hmax: maximum edgelength in the domain (meters) 


    Returns
    -------
        fh: scipy.inerpolate.RegularGridInterpolater representing isotropic mesh sizes in domain


    Example
    ------

    """
    # call the desired mesh size functions
    for key, value in kwargs.items():
        # read in velocity model as a segy file 
        width = max(bbox) 
        depth = min(bbox) 
        vp,nz,nx = utils.ReadVelocityModel(segy, width, depth)
        # build the wavelength mesh size function
        if(key == "wl"):
            alpha_wl = value
            print('INFO: wavelength sizing function is activated...')
            print(' Mesh sizes with be built to resolve the wavelength with ' +
                  str(alpha_wl)+' vertices...')
            maxFreq = 2  # 2 hz by default
            for key2, value2 in kwargs.items():
                if(key2 == "freq"): 
                    maxFreq = value2
                    print(' Mesh sizes will be built for a max. source freq. of ' + str(maxFreq)+" hz...")
            hh_m=np.zeros(shape=(nz, nx))
            hh_m=(vp*maxFreq)/alpha_wl

        # Grade the mesh sizes 

        # Adjust based on the CFL limit 

        # Construct a interpolator object to be queried during mesh generation
        z_vec,x_vec = utils.CreateDomainVectors(nz,nx,depth,width)
        interpolant = RegularGridInterpolator((z_vec,x_vec),hh_m)
    return interpolant,nz,nx


def PlotMeshSizes(nz,nx,depth,width,size_fx,stride=5):
    ''' Plot the isotropic mesh size function '''
    zg,xg = utils.CreateDomainMatrices(nz,nx,depth,width)
    sz1z,sz1x = zg.shape
    sz2 = sz1z*sz1x
    _zg = np.reshape(zg,(sz2,1))
    _xg = np.reshape(xg,(sz2,1))
    hh_m  = size_fx((_zg,_xg))
    hh_m  = np.reshape(hh_m,(sz1z,sz1x))
    plt.pcolormesh(xg[0::stride],zg[0::stride],hh_m[0::stride],edgecolors='none')
    plt.title('Isotropic mesh sizes')
    plt.colorbar(label='mesh size (m)')
    plt.xlabel('x-direction (km)')
    plt.ylabel('z-direction (km)')
    plt.show()
    return plt
