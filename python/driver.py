import edgefx as ef
import dmsh 

fname = 'vel_z6.25m_x12.5m_exact.segy'

fh,nz,nx = ef.edgefx(bbox=(-12.0,0,0,67),segy=fname,
            hmin=2e3,hmax=4e3)
#ef.PlotMeshSizes(nz,nx,-12.0,67.0,fh )

geo = dmsh.Rectangle(-12.0,0.0,0,67.0)

def edge_size(x):
    import matplotlib.pyplot as plt
    print(x[0,:].min())
    print(x[0,:].max())
    print(x[1,:].min())
    print(x[1,:].max())
    plt.plot(x[0,:],x[1,:])
    plt.show()
    sizes = fh((x[0,:],x[1,:]))
    print(sizes)
    return sizes

X, cells = dmsh.generate(geo, edge_size, show=True, verbose=True)
    

