import meshio

import SeismicMesh


def example_3D():
    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/EGAGE_Salt.bin"
    bbox = (-4200, 0, 0, 13520, 0, 13520)
    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        nx=676,
        ny=676,
        nz=210,
        dt=0.001,
        freq=2,
        wl=3,
        hmin=250,
        grade=0.35,
        grad=250,
        domain_ext=1000,
    )

    # Build mesh size function
    ef = ef.build()

    # Save your options so you have a record
    ef.SaveMeshSizeFunctionOptions("EGAGE_Salt")

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(ef, method="qhull")  # or cgal if installed

    # Build the mesh
    points, cells = mshgen.build(nscreen=1, max_iter=30, seed=0)

    # Write to disk (see meshio for more details)
    meshio.write_points_cells(
        "foo3D_V3.vtk", points, [("tetra", cells)],
    )


if __name__ == "__main__":

    example_3D()
