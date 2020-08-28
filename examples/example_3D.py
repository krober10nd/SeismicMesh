import zipfile

import meshio
import numpy as np
from mpi4py import MPI

import SeismicMesh

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# Serial or parallel 3d mesh generation building a mesh roughly 1.6 million cells.
# Warning: In serial, this example takes roughly 20 minutes...


def example_3D():

    if rank == 0:
        # The velocity model was downloaded from here:
        # https://s3.amazonaws.com/open.source.geoscience/open_data/seg_eage_models_cd/Salt_Model_3D.tar.gz

        # Dimensions of model
        nx, ny, nz = 676, 676, 210

        path = "velocity_models/Salt_Model_3D/3-D_Salt_Model/VEL_GRIDS/"
        # Extract Saltf@@ from SALTF.ZIP
        zipfile.ZipFile(path + "SALTF.ZIP", "r").extract("Saltf@@", path=path)

        # Load data
        with open(path + "Saltf@@", "r") as file:
            vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder(">"))
            vp = vp.reshape(nx, ny, nz, order="F")
            vp = np.flipud(vp.transpose((2, 0, 1)))  # z, x and then y
    else:
        vp = np.zeros(shape=(1, 1, 1))
        vp[:] = 1500.0

    # Bounding box describing domain extents (corner coordinates)
    bbox = (-4200, 0, 0, 13520, 0, 13520)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        velocity_grid=vp,
        dt=0.001,
        freq=2,
        wl=5,
        grade=0.25,
        hmin=150,
        hmax=5e3,
        domain_ext=250,
        padstyle="linear_ramp",
    )

    # Build mesh size function (in parallel)
    ef = ef.build()

    # Write to disk for later use
    ef.WriteVelocityModel("EAGE_Salt")

    # Construct a mesh generator object
    mshgen = SeismicMesh.MeshGenerator(ef)

    # Build the mesh
    points, cells = mshgen.build(max_iter=75, axis=1)

    # Do mesh improvement in serial to bound lower dihedral angle to >= 5 degrees
    points, cells = mshgen.build(
        points=points,
        mesh_improvement=True,
        max_iter=50,
        min_dh_bound=5,
    )

    if rank == 0:
        # Write to disk (see meshio package for more details)
        meshio.write_points_cells(
            "EAGE_Salt.vtk",
            points[:, [1, 2, 0]] / 1000.0,
            [("tetra", cells)],
        )

        # Write to gmsh22 format (quite slow)
        # NOTE: SeismicMesh outputs assumes the domain is (z,x,y) so for visualization
        # in ParaView, we swap the axes so it appears as in the (x,y,z) plane.
        meshio.write_points_cells(
            "EAGE_Salt.msh",
            points[:, [1, 2, 0]] / 1000,
            [("tetra", cells)],
            file_format="gmsh22",
            binary=False,
        )


if __name__ == "__main__":

    example_3D()
