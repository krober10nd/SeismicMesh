import meshio
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def example_3D_parallel():
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
        wl=10,
        hmin=100,
    )

    # Build mesh size function (in parallel)
    ef = ef.build(rank=rank)

    # if rank == 0:
    #    # Save your options so you have a record
    #    ef.SaveMeshSizeFunctionOptions("EGAGE_Salt")

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(
        ef, method="qhull"
    )  # parallel currently only works in qhull

    # Build the mesh (note the seed makes the result deterministic)
    points, cells = mshgen.build(max_iter=50, nscreen=1, seed=10, COMM=comm, axis=2)

    if rank == 0:
        # Write to disk (see meshio for more details)
        meshio.write_points_cells(
            "foo3D_V3.vtk", points, [("tetra", cells)],
        )


if __name__ == "__main__":

    example_3D_parallel()
