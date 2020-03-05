import numpy as np 

def hj(dims,elen,dfdx,imax,ffun):
    ''' PDE-gradient limit a function defined on a 
        structured grid
    '''
    eps = np.finfo(float).eps
    ftol = np.min(ffun)*np.sqrt(eps) 
    aset = np.zeros(shape=(np.prod(dims)))-1
    stencil = np.zeros(4,dtype=int)
    for iter in range(imax):
        active = np.where(aset == iter-1)
        for inod in np.nditer(active):
            (i,j) = np.unravel_index(inod,dims)
            stencil[0] = np.ravel_multi_index((i+1,j),dims,mode='wrap')
            stencil[1] = np.ravel_multi_index((i,j+1),dims,mode='wrap')
            stencil[2] = np.ravel_multi_index((i-1,j),dims,mode='wrap')
            stencil[3] = np.ravel_multi_index((i,j-1),dims,mode='wrap')
            for jnod in stencil: 
                if(ffun[inod] > ffun[jnod]):
                    fun1 = ffun[jnod] + elen * dfdx
                    if(ffun[inod] > fun1 + ftol): 
                        ffun[inod] = fun1
                        aset[inod] = iter 
                else:
                    fun2 = ffun[inod] + elen*dfdx
                    if(ffun[jnod] > fun2 + ftol):
                        ffun[jnod] = fun2
                        aset[jnod] = iter 

        print("ITERATION : "+str(iter))
    return ffun
