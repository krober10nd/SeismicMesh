import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_2d_SDF():
    """Unit circle"""
    hmin = 0.01

    def SDF(p):
        return np.sqrt((p ** 2).sum(1)) - 1.0

    def EF(p):
        d = SDF(p)
        h = hmin - d * 0.15
        return h

    mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=(-1, 1, -1, 1), fd=SDF, fh=EF)

    points, cells = mshgen.build(max_iter=60, perform_checks=True)
    print(len(points), len(cells))
    # 3478 6406
    assert np.abs(len(points) - 3478) < 20
    assert np.abs(len(cells) - 6406) < 20


if __name__ == "__main__":
    test_2d_SDF()
