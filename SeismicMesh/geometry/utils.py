import numpy as np

from . import signed_distance_functions as sdf

"""
Routines to perform geometrical operations on mesh

"""

def remove_external_faces(points, faces, extents):
    """
    Remove faces with all three vertices outside block (external)
    and points that are very far from local domain
    """
    signed_distance = sdf.drectangle(
        points[faces.flatten(), :],
        x1=extents[0],
        x2=extents[2],
        y1=extents[1],
        y2=extents[3],
    )
    isOut = np.reshape(signed_distance > 0, (-1, 3))
    isFar = np.reshape(signed_distance > 1000, (-1, 3))
    faces_new = faces[(np.sum(isOut, axis=1) != 3) & (np.any(isFar, axis=1) != 1), :]

    points_new, faces_new, jx = fixmesh(points, faces_new)

    return points_new, faces_new, jx


def vertex_to_elements(points, faces):
    """
    determine which elements connected to vertices
    """
    num_faces = len(faces)

    ext = np.tile(np.arange(0, num_faces), (3, 1)).reshape(-1, order="F")
    ve = np.reshape(faces, (-1,))
    ve = np.vstack((ve, ext)).T
    ve = ve[ve[:, 0].argsort(), :]

    idx = np.insert(np.diff(ve[:, 0]), 0, 0)
    vtoe_pointer = np.argwhere(idx)
    vtoe_pointer = np.insert(vtoe_pointer, 0, 0)
    vtoe_pointer = np.append(vtoe_pointer, num_faces * 3)

    vtoe = ve[:, 1]

    return vtoe, vtoe_pointer


def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]
    Returns I if return_index is True
    Returns J if return_inverse is True
    """
    A = np.require(A, requirements="C")
    assert A.ndim == 2, "array must be 2-dim'l"

    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize * ncolumns))
    B, I, J = np.unique(A.view(dtype), return_index=True, return_inverse=True)

    B = B.view(orig_dtype).reshape((-1, ncolumns), order="C")

    # There must be a better way to do this:
    if return_index:
        if return_inverse:
            return B, I, J
        else:
            return B, I
    else:
        if return_inverse:
            return B, J
        else:
            return B


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

def fixmesh(p, t, ptol=2e-13):
    """Remove duplicated/unused nodes and
       ensure orientation of elements is CCW
    Parameters
    ----------
    p : array, shape (np, dim)
    t : array, shape (nt, nf)
    Usage
    -----
    p, t = fixmesh(p, t, ptol)
    """
    snap = (p.max(0) - p.min(0)).max() * ptol
    _, ix, jx = unique_rows(np.round(p / snap) * snap, True, True)

    p = p[ix]
    t = jx[t]

    t = np.sort(t, axis=1)
    t = unique_rows(t)

    flip = simpvol(p,t)<0
    t[flip, :2] = t[flip, 1::-1]

    return p, t, jx

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
