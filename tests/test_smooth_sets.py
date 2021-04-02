import numpy as np
import pytest

import SeismicMesh as sm


@pytest.mark.serial
def test_smooth_diff():
    cube1 = sm.Cube((-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
    ball1 = sm.Ball((0.0, 0.0, 0.5), 0.85)

    domain = sm.Difference([ball1, cube1], smoothness=0.20)
    points, cells = sm.generate_mesh(
        domain=domain,
        edge_length=0.10,
    )
    points, cells = sm.sliver_removal(
        points=points,
        domain=domain,
        edge_length=0.10,
    )
    assert np.abs(cells.shape[0] - 9004) < 100


if __name__ == "__main__":
    test_smooth_diff()
