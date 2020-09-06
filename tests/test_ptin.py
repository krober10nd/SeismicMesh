import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_ptin():
    # a perfect tetra
    points = np.array(
        [[0, 0, 0], [0, 1, 0], [np.sqrt(2), 0.5, 0], [0.5, 0.5, np.sqrt(2)]],
        dtype=np.float,
    )

    cells = np.array(
        [[0, 1, 2, 3]],
        dtype=np.int,
    )

    # test point if in tetraheadral
    queryIn = (0.5, 0.5, 1.0)  # should be in
    queryOut = (0.0, 0.0, 1.0)  # should be out

    isIn = SeismicMesh.geometry.utils.vertex_in_entity3(
        queryIn,
        (
            *points[cells[0][0], 0:3],
            *points[cells[0][1], 0:3],
            *points[cells[0][2], 0:3],
            *points[cells[0][3], 0:3],
        ),
    )
    assert isIn != 0

    isIn = SeismicMesh.geometry.utils.vertex_in_entity3(
        queryOut,
        (
            *points[cells[0][0], 0:3],
            *points[cells[0][1], 0:3],
            *points[cells[0][2], 0:3],
            *points[cells[0][3], 0:3],
        ),
    )
    assert isIn != 1


if __name__ == "__main__":
    test_ptin()
