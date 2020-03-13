import meshio

import SeismicMesh

fname = "vel_z6.25m_x12.5m_exact.segy"

# Construct mesh sizing object
ef = SeismicMesh.MeshSizeFunction(
    bbox=(-12e3, 0, 0, 67e3), segy=fname, wl=5, grade=5.0, hmin=100
)

# Build mesh size function
ef = ef.build()

# Visualize mesh size function
# ef.plot()

# Construct mesh generator
mshgen = SeismicMesh.MeshGenerator(ef)

points, facets = mshgen.build(max_iter=100)

meshio.write_points_cells(
    "foo.vtk",
    points,
    [("triangle", facets)],
    # Optionally provide extra data on points, cells, etc.
    # point_data=point_data,
    # cell_data=cell_data,
    # field_data=field_data
)
