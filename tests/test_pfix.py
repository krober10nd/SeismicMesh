import numpy as np
import pytest

import SeismicMesh


@pytest.mark.serial
def test_pfix():
    hmin = 0.05
    bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    pfix = np.linspace((0.0, 0.0, 0.0), (1.0, 0.0, 1.0), int(np.sqrt(2) / hmin))

    pfix = np.vstack((pfix, SeismicMesh.geometry.corners(bbox)))

    cube = SeismicMesh.Cube(bbox)

    points, cells = SeismicMesh.generate_mesh(domain=cube, edge_length=hmin, pfix=pfix)

    def _closest_node(node, nodes):
        nodes = np.asarray(nodes)
        deltas = nodes - node
        dist_2 = np.einsum("ij,ij->i", deltas, deltas)
        return dist_2

    for p in pfix:
        d = _closest_node(p, points)
        assert np.isclose(np.min(d), 0.0)


if __name__ == "__main__":
    test_pfix()
