import segyio
import numpy as np
import sys,io

def ReadVelocityModel(segy,depth,width):
    ''' Read a velocity model from disk in a segy format '''
    FOUND = False
    with segyio.open(segy, ignore_geometry=True) as f:
        FOUND = True
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
            print('Exiting...segy file not found...')
            sys.exit(1)
    return vp, nz, nx

def CreateDomainVectors(nz,nx,depth,width):
    xvec = np.linspace(0, width, nx)
    zvec = np.linspace(depth, 0, nz)
    return zvec,xvec

def CreateDomainMatrices(nz,nx,depth,width):
    zvec,xvec = CreateDomainVectors(nz,nx,depth,width)
    zg, xg = np.meshgrid(zvec, xvec, indexing='ij')
    return zg,xg

