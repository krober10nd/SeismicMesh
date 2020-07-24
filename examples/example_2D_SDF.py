import meshio
from scipy.interpolate import RegularGridInterpolator
import numpy as np

import SeismicMesh


def example_2D_SDF():
    bbox = (0, 0.6, 0, 1)
    hmin = 0.005

    # build a simple signed distance function
    xvec = np.linspace(0, 0.6, 110)
    yvec = np.linspace(0, 1, 180)
    # load in the SDF data from a file.
    # vals = np.loadtxt("simple_sdf.txt", delimiter=",")
    vals = np.load("last_phi_mat.npy")
    interpolant = RegularGridInterpolator(
        (xvec, yvec), vals, bounds_error=False, fill_value=None
    )

    def SDF(p):
        return interpolant(p)

    # build a simple sizing function
    def EF(p):
        d = SDF(p)
        h = hmin - d * 0.15
        return h

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(hmin=hmin, bbox=bbox, fd=SDF, fh=EF)

    # Build the mesh (note the seed makes the result deterministic)
    points, facets = mshgen.build(max_iter=50)

    meshio.write_points_cells(
        "simple_irregular.vtk",
        points / 1000,
        [("triangle", facets)],
        file_format="vtk",
    )


if __name__ == "__main__":

    example_2D_SDF()
