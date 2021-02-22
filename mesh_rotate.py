import meshio
import numpy as np
import SeismicMesh

# bbox = (0.0, 1.0, 0.0, 1.0)
# rotations = np.linspace(-3.14, 3.14, 40)
# squares = []
# for _, rotate in enumerate(rotations):
#    squares.append(SeismicMesh.Rectangle(bbox, rotate=rotate))
#
# rotated_squares = SeismicMesh.Union(squares)
#
# points, cells = SeismicMesh.generate_mesh(domain=rotated_squares, edge_length=0.05)
# meshio.write_points_cells(
#    "square" + str(rotate) + ".vtk",
#    points,
#    [("triangle", cells)],
#    file_format="vtk",
# )

bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
rotations = np.linspace(-3.14, 3.14, 40)
cubes = []
for _, rotate in enumerate(rotations):
    cubes.append(SeismicMesh.Cube(bbox, rotate=rotate))

rotated_cubes = SeismicMesh.Union(cubes)

points, cells = SeismicMesh.generate_mesh(domain=rotated_cubes, edge_length=0.10)
meshio.write_points_cells(
    "rotated_cubes.vtk",
    points,
    [("tetra", cells)],
    file_format="vtk",
)
