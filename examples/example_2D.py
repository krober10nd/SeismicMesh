import meshio
import numpy as np

import SeismicMesh


def example_2D():
    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        model=fname,
        domain_ext=1e3,
        dt=0.001,
        grade=0.15,
        grad=50,
        freq=2,
        wl=3,
        hmax=1e3,
        hmin=50.0,
        padstyle="linear_ramp",
    )

    # Build mesh size function
    ef = ef.build()

    ef.WriteVelocityModel("BP2004")

    # Visualize mesh size function
    ef.plot()

    # Construct mesh generator
    mshgen = SeismicMesh.MeshGenerator(
        ef, method="cgal"
    )  # if you have cgal installed, you can use method="cgal"

    # Build the mesh (note the seed makes the result deterministic)
    points, facets = mshgen.build(max_iter=50, nscreen=1, seed=0)

    # Apply boundary conditions
    # 1) top of domain is labeled 11
    # 2) the three sides are labeled 10
    ordered_bnde = SeismicMesh.geometry.get_winded_boundary_edges(facets)
    ordered_bpts = points[ordered_bnde.flatten(), :]
    tmp = np.argwhere(ordered_bpts[:, 0] > -50)
    ix1, ix2 = int(tmp[0] / 2), int(tmp[-1] / 2)
    top = ordered_bnde[ix1:ix2, :]
    sides_bot = np.concatenate(
        (ordered_bnde[ix2 + 1 : -1, :], ordered_bnde[0 : ix1 - 1, :])
    )

    meshio.write_points_cells(
        "BP2004_w1KM_EXT.vtk", points / 1000, [("triangle", facets)], file_format="vtk",
    )

    # Write to gmsh22 format with boundary conditions
    meshio.write_points_cells(
        "BP2004_w1KM_EXT.msh",
        points / 1000,
        cells=[("triangle", facets)],
        file_format="gmsh22",
        binary=False,
    )
    quit()
    meshio.write_points_cells(
        "BP2004_w1KM_EXT.msh",
        points / 1000,
        cells=[
            ("triangle", facets),
            ("line", np.array(top)),
            ("line", np.array(sides_bot)),
        ],
        field_data={
            "TopEdges": np.array([11, 1]),
            "SidesBottomEdges": np.array([10, 1]),
        },
        cell_data={
            "gmsh:physical": np.array(
                [
                    np.repeat(3, len(facets)),
                    np.repeat(11, len(top)),
                    np.repeat(10, len(sides_bot)),
                ]
            ),
            "gmsh:geometrical": np.array(
                [
                    np.repeat(1, len(facets)),
                    np.repeat(1, len(top)),
                    np.repeat(1, len(sides_bot)),
                ]
            ),
        },
        file_format="gmsh22",
        binary=False,
    )


if __name__ == "__main__":

    example_2D()
