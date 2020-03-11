import sys
import io

import numpy as np
import matplotlib.pyplot as plt
import segyio
import distmesh as dm # for signed distnace function
from scipy.interpolate import RegularGridInterpolator

import FastHJ 


def PlotMeshSizes(nz,nx,depth,width,fh,stride=5):
    ''' Plot the isotropic mesh size function

    Usage
    -------
    >>>> PlotMeshSizes(nz,nx,depth,width,fh)


    Parameters
    -------
        nz:  number of points in the depth dimension
        nx:  number of points in the width dimension
        depth: min. z (- in km)
        width: max. x (+ in km)
        fh: mesh size function (SciPy.RegularGriddedInterpolant)
                                **kwargs
        stride: downsample the image by n (n=5 by default)

     Returns
    -------
        none
    '''
    zg,xg = utils.CreateDomainMatrices(nz,nx,depth,width)
    sz1z,sz1x = zg.shape
    sz2 = sz1z*sz1x
    _zg = np.reshape(zg,(sz2,1))
    _xg = np.reshape(xg,(sz2,1))
    hh_m = fh((_zg,_xg))
    hh_m = np.reshape(hh_m,(sz1z,sz1x))
    plt.pcolormesh(xg[0::stride],zg[0::stride],hh_m[0::stride],edgecolors='none')
    plt.title('Isotropic mesh sizes')
    plt.colorbar(label='mesh size (m)')
    plt.xlabel('x-direction (km)')
    plt.ylabel('z-direction (km)')
    plt.show()
    return

class EdgeFunction:
    """
    EdgeFunction: build an isotropic mesh size function for seismic domains.
    
    Usage
    -------
    >>>> fh,fd,nz,nx = EdgeFunction(bbox,hmin,segy,**kwargs)
    
    
    Parameters
    -------
        bbox: Bounding box, (xmin, xmax, ymin, ymax)
        hmin: Minimum edgelength populating domain, (meters)
        segy: Segy file containing velocity model, (assumes velocity is in METERS PER SECOND)
    
                                **kwargs
        wl:   number of nodes per wavelength for given max. freq, (num. of nodes per wl.)
        freq: maximum source frequency for which to estimate wl (hertz)
        hmax: maximum edgelength in the domain (meters)
        dt: maximum stable timestep (in seconds given Courant number cr)
        cr_max: dt is theoretically stable with this Courant number (default 0.2)
        grade: maximum allowable variation in mesh size (default 0.15)
    
    
    Returns
    -------
        fh: lambda function w/ scipy.inerpolate.RegularGridInterpolater representing isotropic mesh sizes in domain
        fd: lambda function representing the signed distance function of domain
    
    
    Example
    ------
    
    fh,nz,nx = ef.edgefx(
            bbox=(-12e3,0,0,67e3),
            ,segy=fname,
            hmin=2e3
            #,wl=5,freq=5 
            #,hmax=4e3
            #,grade=10.0
            #,dt=0.001
            #,cr_max=0.1
            )
    
    """

    def __init__(
            self, bbox, hmin, segy, 
            wl=0.0, freq=5.0, hmax=np.inf, 
            dt=0.0, cr_max=0.2, grade=0.0
            ):
        self.__bbox=bbox 
        self.__hmin=hmin 
        self.__segy=segy 
        self.__wl=wl 
        self.__freq=freq 
        self.__hmax=hmax 
        self.__dt=dt 
        self.__cr_max=cr_max 
        self.__grade=grade 

    # TODO add sanity checks for getters/setters 

    @property 
    def bbox(self):
        return self.__bbox

    @bbox.setter 
    def bbox(self, value): 
        self.__bbox = value

    @property
    def hmin(self): 
        return self.__hmin

    @hmin.setter 
    def hmin(self, value):
        self.hmin = value

    @property 
    def segy(self): 
        return self.__segy 

    @segy.setter 
    def segy(self, value): 
        self.segy = value 

    @property 
    def wl(self): 
        return self.__wl 

    @wl.setter 
    def wl(self, value): 
        self.wl = value 

    @property 
    def freq(self): 
        return self.__freq 

    @freq.setter 
    def freq(self,value): 
        self.__freq = value 

    @property 
    def hmax(self): 
        return self.__hmax 
    
    @hmax.setter 
    def hmax(self,value): 
        self.__hmax = value 

    @property 
    def dt(self): 
        return self.__dt 

    @dt.setter 
    def dt(self,value): 
        self.__dt = value 

    @property 
    def cr_max(self): 
        return self.__cr_max 

    @cr_max.setter 
    def cr_max(self,value): 
        self.__cr_max = value 

    @property 
    def grade(self): 
        return self.__grade 

    @grade.setter 
    def grade(self,value): 
        self.__grade = value 

    ### PUBLIC METHODS ###

    ### PRIVATE METHODS ###

    def __read_velocity_model(self):
        """ Read a velocity model stored as a segy file. 
            Uses the python package segyio. 
        """
        _fname = self.segy 
        _bbox = self.bbox 

        width = max(_bbox)
        depth = min(_bbox)

        found = False
        with segyio.open(_fname, ignore_geometry=True) as f:
            found = True
            # determine length of velocity model from file
            nz = len(f.samples)
            nx = len(f.trace)
            vp = np.zeros(shape=(nz, nx))
            index = 0
            for trace in f.trace:
                vp[:,index]=trace # convert to m-s?
                index += 1
            vp = np.flipud(vp)
            if not FOUND:
                print('Exiting...segy file called '+_fname+'not found...')
                sys.exit(1)
            return vp, nz, nx 

    def _CreateDomainVectors(nz,nx,depth,width):
        xvec = np.linspace(0, width, nx)
        zvec = np.linspace(depth, 0, nz)
        return zvec,xvec
    
    def _CreateDomainMatrices(nz,nx,depth,width):
        zvec,xvec = _CreateDomainVectors(nz,nx,depth,width)
        zg, xg = np.meshgrid(zvec, xvec, indexing='ij')
        return zg,xg


    #    # call the desired mesh size functions
    #    hh_m=np.zeros(shape=(nz, nx)) + hmin
    #    if(alpha_wl > 0):
    #        print('Mesh sizes with be built to resolve an estimate of wavelength with ' +
    #          str(alpha_wl)+' vertices...')
    #        hh_m=vp/(max_freq*alpha_wl)
    #    # enforce min (and optionally max) sizes
    #    hh_m=np.where(hh_m<hmin, hmin, hh_m)
    #    if(hmax < np.inf):
    #        print('Enforcing maximum mesh resolution...')
    #        hh_m=np.where(hh_m>hmax, hmax, hh_m)
    #    # grade the mesh sizes
    #    if(grade > 0):
    #        print('Enforcing mesh gradation...')
    #        elen=width/nx
    #        imax=10000
    #        ffun=hh_m.flatten('F')
    #        ffun_list = ffun.tolist()
    #        tmp = FastHJ.limgrad([nz,nx,1],elen,grade,imax,ffun_list)
    #        tmp = np.asarray(tmp)
    #        hh_m = np.reshape(tmp,(nz,nx),'F')
    #    # adjust based on the CFL limit so cr < cr_max
    #    if( dt > 0 ):
    #        print('Enforcing timestep of '+str(dt)+' seconds...')
    #        cr_old = (vp*dt)/hh_m
    #        dxn = (vp*dt)/cr_max
    #        hh_m = np.where( cr_old > cr_max, dxn, hh_m)
    #    # construct a interpolator object to be queried during mesh generation
    #    z_vec,x_vec = utils.CreateDomainVectors(nz,nx,depth,width)
    #    assert np.all(hh_m > 0.0),"edge_size_function must be strictly positive."
    #    interpolant = RegularGridInterpolator((z_vec,x_vec),hh_m,bounds_error=False)
    #    fh = lambda p: fh(p)
    #
    #
    #    return fh, fd nz, nx 
    #
    

