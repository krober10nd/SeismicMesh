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
        wl=5,
        grade=0.25,
        hmin=100,
        hmax=5e3,
        domain_ext=250,
        padstyle="linear_ramp",
    )

    # Build mesh size function (in parallel)
    ef = ef.build(comm=comm)
    # Build lambda functions
    ef = ef.construct_lambdas(comm)

    if rank == 0:
        ef.WriteVelocityModel("EGAGE_Salt")

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(ef, method="cgal")

    # Build the mesh (note the seed makes the result deterministic)
    points, cells = mshgen.build(max_iter=75, seed=0, COMM=comm)

    if rank == 0:

        # Do mesh improvement in serial to bound lower dihedral angle
        # mshgen.method = "cgal"
        #points, cells = mshgen.build(
        #    mesh_improvement=True,
        #    points=points,
        #    max_iter=50,
        #    min_dh_bound=5,
        #)
        # Write to disk (see meshio for more details)
        meshio.write_points_cells(
            "EGAGE_Salt.vtk", points / 1000.0, [("tetra", cells)],
        )
        ## Write to gmsh22 format (quite slow)
        meshio.write_points_cells(
            "EGAGE_Salt.msh",
            points / 1000,
            [("tetra", cells)],
            file_format="gmsh22",
            binary=False,
        )


if __name__ == "__main__":

    example_3D_parallel()
