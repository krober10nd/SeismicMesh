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
    wl = 10
    hmin = 50
    freq = 2
    grade = 0.005
    nz, nx, ny = 20, 10, 10
    bbox = (-2000.0, 0.0, 0.0, 1000.0, 0.0, 1000.0)
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        nx=nx,
        ny=ny,
        nz=nz,
        endianness="big",
        grade=grade,
        freq=freq,
        wl=wl,
        model=fname,
        hmin=hmin,
    )
    ef = ef.build()

    SDF = SdfGen(
        bbox=bbox, field=ef.vp, min_threshold=2000.0, gridspacing=(100.0, 100.0, 100.0),
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
