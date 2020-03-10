import numpy as np
import utils


def edgefx(bbox, hmin, segy,  **kwargs):
    from scipy.interpolate import RegularGridInterpolator
    import FastHJ 
    """
    edgefx: build a mesh size function for seismic problems

    Usage
    -------
    >>>> fh = edgefx(bbox,hmin,segy,**kwargs)


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
        fh: scipy.inerpolate.RegularGridInterpolater representing isotropic mesh sizes in domain

    Example
    ------

    fh,nz,nx = ef.edgefx(bbox=(-12.0,0,0,67),segy=fname,
            wl=5,freq=5,
            hmin=2e3,hmax=4e3,
            grade=10.0)

    """
    # set reasonable default values
    hmax = np.inf # meters
    maxFreq = 5  # hz
    alpha_wl = 0.0 # no. of nodes per wavelength
    dt = 0.0 # maximum stable timestep
    cr_max = 0.2 # dt is enforced for this Courant number
    grade = 0.0 # maximum variation in mesh size
    # set the values defined by the user
    for key, value in kwargs.items():
        if(key == "wl"):
            alpha_wl = value
        elif(key == "freq"):
            maxFreq = value
        elif(key == "hmax"):
            hmax = value
        elif(key == "dt"):
            dt = value
        elif(key == "cr_max"):
            cr_max = value
        elif(key == "grade"):
            grade = value
    # read in velocity model as a segy file
    width = max(bbox)
    depth = min(bbox)
    vp,nz,nx = utils.ReadVelocityModel(segy, depth, width)
    # call the desired mesh size functions
    hh_m=np.zeros(shape=(nz, nx)) + hmin
    if(alpha_wl > 0):
        print('Mesh sizes with be built to resolve an estimate of wavelength with ' +
          str(alpha_wl)+' vertices...')
        hh_m=(vp*maxFreq)/alpha_wl
    # enforce min (and optionally max) sizes
    hh_m=np.where(hh_m<hmin, hmin, hh_m)
    if(hmax < np.inf):
        print('Enforcing maximum mesh resolution...')
        hh_m=np.where(hh_m>hmax, hmax, hh_m)
    # grade the mesh sizes
    if(grade > 0):
        print('Enforcing mesh gradation...')
        elen=width/nx
        imax=10000
        ffun=hh_m.flatten('F')
        ffun_list = ffun.tolist()
        tmp = FastHJ.limgrad([nz,nx,1],elen,grade,imax,ffun_list)
        tmp = np.asarray(tmp)
        hh_m = np.reshape(tmp,(nz,nx),'F')
    # adjust based on the CFL limit so cr < cr_max
    if( dt > 0 ):
        print('Enforcing timestep of '+str(dt)+' seconds...')
        cr_old = (vp*dt)/hh_m
        dxn = (vp*dt)/cr_max
        hh_m = np.where( cr_old > cr_max, dxn, hh_m)
    # construct a interpolator object to be queried during mesh generation
    z_vec,x_vec = utils.CreateDomainVectors(nz,nx,depth,width)
    assert np.all(hh_m > 0.0),"edge_size_function must be strictly positive."
    interpolant = RegularGridInterpolator((z_vec,x_vec),hh_m,bounds_error=False)
    return interpolant,nz,nx


def PlotMeshSizes(nz,nx,depth,width,fh,stride=5):
    import matplotlib.pyplot as plt
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
        stride: Downsample the image by n (n=5 by default)

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
