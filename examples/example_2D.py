import meshio
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def example_2D():
    """
    Build a mesh of the BP2004 benchmark velocity model in serial or parallel
    Takes roughly 1 minute with 2 processors.
    The velocity model can be downloaded from here: https://wiki.seg.org/wiki/2004_BP_velocity_estimation_benchmark_model
    """

    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    # Read in it
    vp = SeismicMesh.ReadSegy(fname)

    # Bounding box describing domain extents (corner coordinates)
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        velocity_grid=vp,
        freq=2,
        wl=10,
        dt=0.001,
        hmin=75.0,
        grade=0.15,
        domain_ext=1e3,
        padstyle="linear_ramp",
    )

    # Build mesh size function
    ef = ef.build()

    # Write to disk for later use
    ef.WriteVelocityModel("BP2004_w1KM_EXT")

    # Visualize the mesh size function
    ef.plot()

    # Construct a mesh generator object
    mshgen = SeismicMesh.MeshGenerator(ef)

    # Build the mesh
    points, facets = mshgen.build(axis=1)

    if rank == 0:

        # Write the mesh as a vtk format for visualization in Paraview

        # NOTE: SeismicMesh outputs assumes (z,x) so for visualization
        # in ParaView we swap the axes so it appears as in the (x,z) plane

        meshio.write_points_cells(
            "BP2004.vtk",
            points[:, [1, 0]] / 1000,
            [("triangle", facets)],
            file_format="vtk",
        )

        # Write to gmsh22 format (quite slow) used by many numerical solvers
        meshio.write_points_cells(
            "BP2004.msh",
            points[:, [1, 0]] / 1000,
            [("triangle", facets)],
            file_format="gmsh22",
            binary=False,
        )


if __name__ == "__main__":

    example_2D()
