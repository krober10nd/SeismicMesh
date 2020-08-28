import meshio
import numpy as np
from mpi4py import MPI

import SeismicMesh
from SeismicMesh.geometry import (
    SignedDistanceFunctionGenerator as SdfGen,
)  # import the tool used to generate the SDF

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if size > 1:
    print("This is a serial example!")


def example_2D_Irregular():
    """
    Build a 2D mesh of a benchmark velocity  model of a Monuntainous thrust region.
    Builds in serial or parallel. Demonstrates how to ensure the topography is respected
    in the mesh.

    Velocity model was downloaded from here: https://wiki.seg.org/wiki/1994_BP_migration_from_topography
    """

    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/velocity.segy"
    vp = SeismicMesh.ReadSegy(fname)

    # Bounding box describing domain extents (corner coordinates)
    bbox = (-8e3, 2e3, 0, 25e3)
    hmin = 37.5

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        velocity_grid=vp,
        freq=4,
        wl=10,
        dt=0.001,
        hmin=hmin,
        grade=0.15,
    )

    # Build mesh size function
    ef = ef.build()

    # Write to disk for later use
    ef.WriteVelocityModel("Irregular2D")

    # Visualize the mesh size function
    ef.plot()

    # Bulid a signed distance function from the seismic velocity model
    # Some pockets of velocity < 4000 exist, fill those in.
    vp2 = vp.copy()
    vp2 = np.where(vp2 < 4000, 4001, vp2)
    SDF = SdfGen(
        bbox=bbox,
        field=vp2,
        min_threshold=4000.0,
        gridspacing=(10.0, 15.0),
    ).SDF

    # Construct a mesh generator object
    mshgen = SeismicMesh.MeshGenerator(bbox=bbox, hmin=hmin, fh=ef.fh, fd=SDF)

    # Build the mesh
    points, facets = mshgen.build(perform_checks=True)

    if rank == 0:

        # Write the mesh as a vtk format for visualization in ParaView

        # NOTE: SeismicMesh outputs assumes (z,x) so for visualization
        # in ParaView we swap the axes so it appears as in the (x,z) plane

        meshio.write_points_cells(
            "Irregular2D.vtk",
            points[:, [1, 0]] / 1000,
            [("triangle", facets)],
            file_format="vtk",
        )

        # Write to gmsh22 format (quite slow) used by many numerical solvers
        meshio.write_points_cells(
            "Irregular2D.msh",
            points[:, [1, 0]] / 1000,
            [("triangle", facets)],
            file_format="gmsh22",
            binary=False,
        )


if __name__ == "__main__":

    example_2D_Irregular()
