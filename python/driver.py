import edgefx as ef

size_fx,nz,nx = ef.edgefx(bbox=(-3.5,0,0,17),segy='vp_marmousi-ii.segy',
            wl=5,freq=5,
            hmin=1e3,hmax=10e3)

ef.PlotMeshSizes(nz,nx,-3.5,17,size_fx)
