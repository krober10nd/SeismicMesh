import numpy as np


def simpvol(p, t):
    """Signed volumes of the simplex elements in the mesh."""
    dim = p.shape[1]
    if dim == 1:
        d01 = p[t[:, 1]] - p[t[:, 0]]
        return d01
    elif dim == 2:
        d01 = p[t[:, 1]] - p[t[:, 0]]
        d02 = p[t[:, 2]] - p[t[:, 0]]
        return (d01[:, 0] * d02[:, 1] - d01[:, 1] * d02[:, 0]) / 2
    else:
        raise NotImplementedError


def simpqual(p, t):
    """Simplex quality.
    Usage
    -----
    q = simpqual(p, t)
    Parameters
    ----------
    p : array, shape (np, dim)
        nodes
    t : array, shape (nt, dim+1)
        triangulation
    Returns
    -------
    q : array, shape (nt, )
        qualities
    """
    assert p.ndim == 2 and t.ndim == 2 and p.shape[1] + 1 == t.shape[1]

    def length(p1):
        return np.sqrt((p1 ** 2).sum(1))

    a = length(p[t[:, 1]] - p[t[:, 0]])
    b = length(p[t[:, 2]] - p[t[:, 0]])
    c = length(p[t[:, 2]] - p[t[:, 1]])
    r = 0.5 * np.sqrt((b + c - a) * (c + a - b) * (a + b - c) / (a + b + c))
    R = a * b * c / np.sqrt((a + b + c) * (b + c - a) * (c + a - b) * (a + b - c))
    return 2 * r / R


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


def get_boundary_vertices2(faces):
    """
    Get the unique boundary indices of the mesh
    """
    bedges = get_boundary_edges_of_mesh2(faces)
    indices = np.unique(bedges.reshape(-1))
    return indices


def are_boundary_vertices2(points, faces):
    """
    Return array of 1 or 0 if vertex is boundary vertex or not
    """
    ix = get_boundary_vertices2(faces)
    result = np.zeros((len(points), 1), dtype=int)
    result[ix] = 1
    return result
