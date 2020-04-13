import numpy as np

""" Utilities for meshes """


def get_boundary_edges_of_mesh2(faces):
    """
    Get the boundary edges of the mesh
    """
    edges = get_edges_of_mesh2(faces)
    edges = np.sort(edges, axis=1)
    edges = edges[edges[:, 0].argsort(), :]
    diff = np.diff(edges, axis=0)
    idx = [True if cond1 and cond2 else False for cond1, cond2 in diff]
    tmp1 = np.append(idx, False)
    tmp2 = np.insert(idx, 0, False)
    idx = [True if cond1 or cond2 else False for cond1, cond2 in zip(tmp1, tmp2)]
    boundary_edges = np.array([edge for edge, res in zip(edges, idx) if not res])
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
