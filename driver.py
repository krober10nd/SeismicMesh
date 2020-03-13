import SeismicMesh as sm

fname = 'velocity_models/vel_z6.25m_x12.5m_exact.segy'

# Construct mesh sizing object
ef = sm.MeshSizeFunction(
    bbox=(-12e3, 0, 0, 67e3),
    wl=5,
    segy=fname,
    hmin=10)

# Build mesh size function
ef = ef.build()

# Visualize mesh size function
ef.plot()

# Construct mesh generator
mshgen = sm.MeshGenerator(SizingFunction=ef)

mshgen = mshgen.build()


