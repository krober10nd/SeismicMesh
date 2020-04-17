import numpy as np


def get_edges_of_mesh2(faces):
    """
    Get the edges of 2D triangular mesh in no order.
    """
    faces = np.array(faces)
    edges = faces[:, [0, 1]]
    edges = np.append(edges, faces[:, [0, 2]], axis=0)
    edges = np.append(edges, faces[:, [1, 2]], axis=0)
    return edges


def get_boundary_edges_of_mesh2(faces):
    """
    Get the boundary edges of the mesh
    """
    edges = get_edges_of_mesh2(faces)
    edges = np.sort(edges, axis=1)
    unq, cnt = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_edges


def get_winded_boundary_edges_of_mesh2(faces):
    """
    Order the boundary edges of the mesh in a winding order.
    """
    boundary_edges = get_boundary_edges_of_mesh2(faces)
    _bedges = boundary_edges.copy()

    choice = 0
    isVisited = np.zeros((len(_bedges)))
    ordering = np.array([choice])
    isVisited[choice] = 1

    vStart, vNext = _bedges[choice, :]
    while True:
        locs = np.column_stack(np.where(_bedges == vNext))
        rows = locs[:, 0]
        choices = [row for row in rows if isVisited[row] == 0]
        if len(choices) == 0:
            break
        choice = choices[0]
        ordering = np.append(ordering, [choice])
        isVisited[choice] = 1
        nextEdge = _bedges[choice, :]
        tmp = [v for v in nextEdge if v != vNext]
        vNext = tmp[0]
    boundary_edges = boundary_edges[ordering, :]
    return boundary_edges
