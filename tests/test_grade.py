import os

import numpy as np
import pytest

import SeismicMesh


def CheckGrade(hh, nz, nx, elen):
    max_grade = 0.0
    for i in range(1, nz - 1):
        for j in range(1, nx - 1):
            val = hh[i, j]
            nval = hh[i, j + 1]
            max_grade = np.max([np.abs(nval - val) / elen, max_grade])
            nval = hh[i + 1, j]
            max_grade = np.max([np.abs(nval - val) / elen, max_grade])
            nval = hh[i, j - 1]
            max_grade = np.max([np.abs(nval - val) / elen, max_grade])
            nval = hh[i - 1, j]
            max_grade = np.max([np.abs(nval - val) / elen, max_grade])
    return max_grade


@pytest.mark.serial
def test_grade():
    fname = os.path.join(os.path.dirname(__file__), "testing.segy")
    vp = SeismicMesh.ReadSegy(fname)
    wl = 5
    hmin = 10
    grade = 0.005
    elen = 1e3
    ef = SeismicMesh.MeshSizeFunction(
        bbox=(-10e3, 0, 0, 10e3), grade=grade, wl=wl, velocity_grid=vp, hmin=hmin
    )
    ef = ef.build()

    fh = ef.fh
    zg, xg = ef.GetDomainMatrices()
    nz, nx = zg.shape
    _zg = np.reshape(zg, (nz * nx, 1))
    _xg = np.reshape(xg, (nz * nx, 1))
    hh = fh((_zg, _xg))
    hh = np.reshape(hh, (nz, nx))
    max_grade = CheckGrade(hh, nz, nx, elen=elen)
    assert max_grade <= grade


if __name__ == "__main__":
    test_grade()
