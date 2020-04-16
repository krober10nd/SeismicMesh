import meshio
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def example_2D_parallel():
    """
    Build the BP-2004 model in parallel
    """

    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        domain_ext=500,
        dt=0.001,
        freq=5,
        wl=5,
        hmax=1e3,
        hmin=50.0,
        grade=0.05,
    )

    # Build mesh size function
    ef = ef.build(rank=rank)

    if rank == 0:
        ef.WriteVelocityModel("BP2004")

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(
        ef, method="qhull"
    )  # parallel currently only works in qhull

    # Build the mesh (note the seed makes the result deterministic)
    points, facets = mshgen.build(max_iter=50, nscreen=1, seed=0, COMM=comm)

    # Write to disk (see meshio for more details)
    # Write as a vtk format for visualization in Paraview
    if rank == 0:
        meshio.write_points_cells(
            "BP2004.vtk", points / 1000, [("triangle", facets)], file_format="vtk",
        )
        # Write to gmsh22 format
        meshio.write_points_cells(
            "BP2004.msh",
            points / 1000,
            [("triangle", facets)],
            file_format="gmsh22",
            binary=False,
        )


if __name__ == "__main__":

    example_2D_parallel()
