import numpy as np
import numpy.matlib as npm

# Dimensions
nx, ny, nz = 10, 20, 10
vp_dummy = np.zeros((nz, nx, ny))

# cross sectional profile
ca = npm.repmat(np.linspace(1500, 4500, 20), 10, 1)
vp_dummy[:, :, 0::] = ca
print(vp_dummy.shape)
vp_dummy.astype(dtype="float32", order="F").tofile("test/test3D.bin")

_vp = np.fromfile("test/test3D.bin", dtype=np.dtype("float32")).newbyteorder(">")
_vp = _vp.byteswap()
print(_vp.shape)
print(_vp)
