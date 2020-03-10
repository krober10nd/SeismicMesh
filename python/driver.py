import edgefx as ef
import dmsh 

fname = 'vel_z6.25m_x12.5m_exact.segy'

fh,nz,nx = ef.edgefx(bbox=(-12.0,0,0,67),segy=fname,
            wl=5,freq=5,
            hmin=2e3,hmax=4e3)
ef.PlotMeshSizes(nz,nx,-12.0,67.0,fh )

geo = dmsh.Rectangle(0.0,67.0,-12.0,0.0)

def edge_size(x):
    import matplotlib.pyplot as plt
    import numpy as np
    sizes = fh((x[1,:],x[0,:]))
    sizes = np.nan_to_num(sizes,nan=1e3)
    print(np.count_nonzero(np.isnan(sizes)))
    plt.scatter(x[0,:],x[1,:],c = sizes)
    plt.show()
    return sizes

X, cells = dmsh.generate(geo, edge_size, show=True, verbose=True)
    

