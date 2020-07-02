import meshio
from scipy.interpolate import RegularGridInterpolator
import numpy as np

import SeismicMesh


def example_2D_SDF():
    bbox = (0, 1, 0, 1)
    hmin = 0.01

    # build a simple signed distance function
    xvec = np.linspace(0, 1, 55)
    yvec = np.linspace(0, 1, 90)
    # load in the SDF data from a file.
    vals = np.loadtxt("simple_sdf.txt", delimiter=",")
    interpolant = RegularGridInterpolator(
        (xvec, yvec), vals, bounds_error=False, fill_value=None
    )

    def SDF(p):
        return interpolant(p)

    # build a simple sizing function
    def EF(p):
        print(len(p))
        return np.repeat(hmin, len(p))

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(method="cgal")

    # Build the mesh (note the seed makes the result deterministic)
    points, facets = mshgen.build(h0=hmin, bbox=bbox, fd=SDF, fh=EF, max_iter=20)

    meshio.write_points_cells(
        "simple_irregular.vtk",
        points / 1000,
        [("triangle", facets)],
        file_format="vtk",
    )


if __name__ == "__main__":

    example_2D_SDF()
