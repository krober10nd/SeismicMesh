import numpy as np

import pytest

import SeismicMesh

radii = [0.25, 0.30, 0.35]
hmin = 0.05


@pytest.mark.serial
def test_immersion():
    for radius in radii:
        box0 = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))
        disk0 = SeismicMesh.Disk([0.5, 0.5], radius)

        def fh(p):
            return 0.05 * np.abs(disk0.eval(p)) + hmin

        points, cells = SeismicMesh.generate_mesh(
            domain=box0,
            edge_length=fh,
            h0=hmin,
            subdomains=[disk0],
        )

        pmid = points[cells].sum(1) / 3  # Compute centroids
        sd = disk0.eval(pmid)

        # measure the subdomain volume
        subdomain_vol = np.sum(SeismicMesh.geometry.simp_vol(points, cells[sd < 0]))

        assert np.isclose(subdomain_vol, np.pi * radius ** 2, rtol=1e-2)


if __name__ == "__main__":
    test_immersion()
