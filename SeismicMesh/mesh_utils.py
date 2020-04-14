import numpy as np

""" Utilities for meshes """


def get_boundary_edges_of_mesh2(faces):
    """
    Get the boundary edges of the mesh
    """
    edges = get_edges_of_mesh2(faces)
    unq, cnt = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_edges


def get_edges_of_mesh2(faces):
    """
    Get the edges of 2D triangular mesh in no order.
    """
    faces = np.array(faces)
    edges = faces[:, [0, 1]]
    edges = np.append(edges, faces[:, [0, 2]], axis=0)
    edges = np.append(edges, faces[:, [1, 2]], axis=0)
    return edges
