import edgefx as ef

fname = 'vel_z6.25m_x12.5m_exact.segy'
fh,nz,nx = ef.edgefx(bbox=(-12.0,0,0,67),segy=fname,
            wl=5,freq=5,
            cr_max = 0.1, dt=0.10,
            hmin=2e3,hmax=4e3)
ef.PlotMeshSizes(nz,nx,-12.0,67.0,fh)
