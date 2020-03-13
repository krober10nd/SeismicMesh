import os

import numpy as np

import SeismicMesh


def test_hmin():
    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    hmin = 100
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-1e3, 0, 0, 1e3), wl=5, segy=fname, hmin=hmin
    )
    ef = ef.build()

    fh = ef.fh
    zg, xg = ef.GetDomainMatrices()
    sz1z, sz1x = zg.shape
    sz2 = sz1z * sz1x
    _zg = np.reshape(zg, (sz2, 1))
    _xg = np.reshape(xg, (sz2, 1))
    hh = fh((_zg, _xg))
    hh = np.reshape(hh, (sz1z, sz1x))
    assert np.amin(hh) == hmin


if __name__ == "__main__":
    test_hmin()
