from MeshSizeFunction import MeshSizeFunction
from MeshGenerator import MeshGenerator 

fname = 'velocity_models/vel_z6.25m_x12.5m_exact.segy'

# Construct mesh sizing object 
ef = MeshSizeFunction(
    bbox=(-12e3, 0, 0, 67e3),
    wl=5,
    segy=fname,
    hmin=10)

# Build mesh size function
ef = ef.build()

# Visualize mesh size function
#ef.plot()

# Consturct mesh generator 
mshgen = MeshGenerator(SizingFunction=ef)

mshgen = mshgen.build() 


