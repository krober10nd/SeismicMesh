import pytest
import os
import numpy as np

import SeismicMesh
from SeismicMesh.geometry import (
    SignedDistanceFunctionGenerator as SdfGen,
)  # import the tool used to generate the SDF


@pytest.mark.serial
def test_3dmesher():

    fname = os.path.join(os.path.dirname(__file__), "test3D.bin")
    nz, nx, ny = 20, 10, 10

    # Load data
    with open(fname, "r") as file:
        vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder("<"))
        vp = vp.reshape(nx, ny, nz, order="F")
        vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y

    wl = 10
    hmin = 50
    freq = 2
    grade = 0.005
    bbox = (-2000.0, 0.0, 0.0, 1000.0, 0.0, 1000.0)
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox, grade=grade, freq=freq, wl=wl, velocity_grid=vp, hmin=hmin,
    )
    ef = ef.build()

    SDF = SdfGen(
        bbox=bbox, field=vp, min_threshold=2000.0, gridspacing=(100.0, 100.0, 100.0),
    ).SDF

    mshgen = SeismicMesh.MeshGenerator(fh=ef.fh, fd=SDF, hmin=hmin, bbox=bbox)

    points, cells = mshgen.build(nscreen=1, max_iter=50, seed=1)
    # 16479 89862

    points, cells = mshgen.build(points=points, mesh_improvement=True)
    print(len(points), len(cells))

    assert np.abs(len(points) - 16479) < 200
    assert np.abs(len(cells) - 89862) < 200


if __name__ == "__main__":
    test_3dmesher()
