import sys
import io

from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
import segyio
import distmesh as dm  # for signed distnace function

import FastHJ


class MeshSizeFunction:
    """
    MeshSizeFunction: build an isotropic mesh size function for seismic domains.

    Usage
    -------
    >>>> obj = MeshSizeFunction(bbox,hmin,segy,**kwargs)


    Parameters
    -------
        bbox: Bounding box, (xmin, xmax, ymin, ymax)
        hmin: Minimum triangular edgelength populating domain, (meters)
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
        MeshSizeFunction object


    Example
    ------
    from MeshSizeFunction import MeshSizeFunction
    ef = MeshSizeFunction(
            bbox=(-12e3,0,0,67e3),
            ,segy=fname,
            hmin=2e3
            # ,wl=5,freq=5
            # ,hmax=4e3
            # ,grade=10.0
            # ,dt=0.001
            # ,cr_max=0.1
            )

    """

    def __init__(
            self, bbox, hmin, segy,
            wl=0.0, freq=5.0, hmax=np.inf,
            dt=0.0, cr_max=0.2, grade=0.0
    ):
        self.bbox = bbox
        self.hmin = hmin
        self.segy = segy
        self.wl = wl
        self.freq = freq
        self.hmax = hmax
        self.dt = dt
        self.cr_max = cr_max
        self.grade = grade
        self.vp = None
        self.nz = None
        self.nx = None
        self.fh = None
        self.fd = None

    @property
    def bbox(self):
        return self.__bbox

    @bbox.setter
    def bbox(self, value):
        assert len(value) >= 4 and len(
            value) <= 6, "bbox has wrong number of args"
        self.__bbox = value

    @property
    def hmin(self):
        return self.__hmin

    @hmin.setter
    def hmin(self, value):
        assert value > 0.0, "hmin must be non-zero"
        self.__hmin = value

    @property
    def segy(self):
        return self.__segy

    @segy.setter
    def segy(self, value):
        assert isinstance(value, str) is True, "segy must be a filename"
        self.__segy = value
        _vp, _nz, _nx = self.__ReadVelocityModel()
        self.__vp = _vp
        self.__nz = _nz
        self.__nx = _nx

    @property 
    def vp(self): 
        return self.__vp 

    @vp.setter 
    def vp(self,value): 
        # PROVIDE CHECKS ON VP 
        self.__vp = value

    @property 
    def nz(self): 
        return self.__nz 

    @nz.setter
    def nz(self,value):
        assert value is None or value > 0, "nz must be non-zero"
        self.__nz = value

    @property 
    def nx(self): 
        return self.__nx

    @nx.setter
    def nx(self,value):
        assert value is None or value > 0, "nx must be non-zero"
        self.__nx = value 

    @property
    def wl(self):
        return self.__wl

    @wl.setter
    def wl(self, value):
        self.__wl = value

    @property
    def freq(self):
        return self.__freq

    @freq.setter
    def freq(self, value):
        self.__freq = value

    @property
    def hmax(self):
        return self.__hmax

    @hmax.setter
    def hmax(self, value):
        self.__hmax = value

    @property
    def dt(self):
        return self.__dt

    @dt.setter
    def dt(self, value):
        self.__dt = value

    @property
    def cr_max(self):
        return self.__cr_max

    @cr_max.setter
    def cr_max(self, value):
        self.__cr_max = value

    @property
    def grade(self):
        return self.__grade

    @grade.setter
    def grade(self, value):
        self.__grade = value

    ### PUBLIC METHODS ###

    def build(self):
        '''Builds the isotropic mesh size function according
            to the user arguments that were pased.

        Usage
        -------
        >>>> obj = build(self)


        Parameters
        -------
            MeshSizeFunction object

         Returns
        -------
            MeshSizeFunction object with specific fields populated:
                fh: lambda function w/ scipy.inerpolate.RegularGridInterpolater representing isotropic mesh sizes in domain
                fd: lambda function representing the signed distance function of domain

        '''
        _bbox = self.bbox
        _nz = self.nz
        _nx = self.nx
        _vp = self.vp

        _hmax = self.hmax
        _hmin = self.hmin
        _grade = self.grade

        _wl = self.wl
        _freq = self.freq

        _dt = self.dt
        _cr_max = self.cr_max

        width = max(_bbox)
        depth = min(_bbox)

        print(_hmin)
        print(_nz,_nx)
        hh_m = np.zeros(shape=(_nz, _nx)) + _hmin
        if(_wl > 0):
            print('Mesh sizes with be built to resolve an estimate of wavelength with ' +
                  str(_wl)+' vertices...')
            hh_m = _vp/(_freq*_wl)
        # enforce min (and optionally max) sizes
        hh_m = np.where(hh_m < _hmin, _hmin, hh_m)
        if(_hmax < np.inf):
            print('Enforcing maximum mesh resolution...')
            hh_m = np.where(hh_m > _hmax, _hmax, hh_m)
        # grade the mesh sizes
        if(_grade > 0):
            print('Enforcing mesh gradation...')
            elen = width/_nx
            imax = 10000
            ffun = hh_m.flatten('F')
            ffun_list = ffun.tolist()
            tmp = FastHJ.limgrad([_nz, _nx, 1], elen, _grade, imax, ffun_list)
            tmp = np.asarray(tmp)
            hh_m = np.reshape(tmp, _nz, _nx, 'F')
        # adjust based on the CFL limit so cr < cr_max
        if(_dt > 0):
            print('Enforcing timestep of '+str(dt)+' seconds...')
            cr_old = (_vp*_dt)/hh_m
            dxn = (_vp*_dt)/_cr_max
            hh_m = np.where(cr_old > _cr_max, dxn, hh_m)
        # construct a interpolator object to be queried during mesh generation
        z_vec, x_vec = __CreateDomainVectors(_nz, _nx, depth, width)
        assert np.all(hh_m > 0.0), "edge_size_function must be strictly positive."
        interpolant = RegularGridInterpolator(
            (z_vec, x_vec), hh_m, bounds_error=False)

        # create a mesh size function interpolant
        self.fh = lambda p: interpolant(p)

        # create a signed distance function

    def plot(self, stride=5):
        ''' Plot the isotropic mesh size function

        Usage
        -------
        >>>> plot(self)


        Parameters
        -------
            self: MeshSizeFunction object
                        **kwargs
            stride: downsample the image by n (n=5 by default)

         Returns
        -------
            none
            '''
        _nz = self.nz
        _nx = self.nx
        _bbox = self.bbox
        width = max(_bbox)
        depth = min(_bbox)
        fh = self.fh 

        zg, xg=__CreateDomainMatrices(_nz, _nx, depth, width)
        sz1z, sz1x=zg.shape
        sz2=sz1z*sz1x
        _zg=np.reshape(zg, (sz2, 1))
        _xg=np.reshape(xg, (sz2, 1))
        hh_m=fh((_zg, _xg))
        hh_m=np.reshape(hh_m, (sz1z, sz1x))
        plt.pcolormesh(xg[0::stride], zg[0::stride],
                       hh_m[0::stride], edgecolors = 'none')
        plt.title('Isotropic mesh sizes')
        plt.colorbar(label = 'mesh size (m)')
        plt.xlabel('x-direction (km)')
        plt.ylabel('z-direction (km)')
        plt.show()

    ### PRIVATE METHODS ###
    def __ReadVelocityModel(self):
        """ uses the python package segyio."""
        _fname=self.segy
        _bbox=self.bbox

        width=max(_bbox)
        depth=min(_bbox)

        found=False
        with segyio.open(_fname, ignore_geometry = True) as f:
            found=True
            # determine length of velocity model from file
            nz=len(f.samples)
            nx=len(f.trace)
            vp=np.zeros(shape = (nz, nx))
            index=0
            for trace in f.trace:
                vp[:, index]=trace  # convert to m-s?
                index += 1
            vp=np.flipud(vp)
        if not found:
            print('Exiting...segy file called '+_fname+'not found...')
            sys.exit(1)
        return vp, nz, nx

    def __CreateDomainVectors(nz, nx, depth, width):
        xvec=np.linspace(0, width, nx)
        zvec=np.linspace(depth, 0, nz)
        return zvec, xvec

    def __CreateDomainMatrices(nz, nx, depth, width):
        zvec, xvec=_CreateDomainVectors(nz, nx, depth, width)
        zg, xg=np.meshgrid(zvec, xvec, indexing = 'ij')
        return zg, xg
