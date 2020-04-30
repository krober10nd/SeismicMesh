import scipy.sparse as spsparse
import numpy as np

from . import signed_distance_functions as sdf

"""
Routines to perform geometrical operations on meshes
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
    Determine which elements connected to vertices
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


def fixmesh(p, t, ptol=2e-13, deldup=False):
    """
    Remove duplicated/unused nodes and
    ensure orientation of elements is CCW
    Parameters
    ----------
    p : array, shape (np, dim)
    t : array, shape (nt, nf)
    Usage
    -----
    p, t = fixmesh(p, t, ptol)
    """
    # duplicate vertices
    snap = (p.max(0) - p.min(0)).max() * ptol
    _, ix, jx = unique_rows(np.round(p / snap) * snap, True, True)

    p = p[ix]
    t = jx[t]

    # duplicate elements
    t = np.sort(t, axis=1)
    t = unique_rows(t)

    # delete unused vertices
    if deldup:
        pix,_,jx = np.unique(t,return_index=True, return_inverse=True)
        t=np.reshape(jx,(t.shape))
        p=p[pix,:]
        pix=ix[pix]

    # element orientation
    flip = simpvol(p, t) < 0
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
    Get the boundary edges of the mesh.
    """
    edges = get_edges_of_mesh2(faces)
    edges = np.sort(edges, axis=1)
    unq, cnt = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_edges


def get_winded_boundary_edges_of_mesh2(faces):
    """
    Order the boundary edges of the mesh in a winding fashion.
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
    Get the indices of the mesh representing boundary vertices.
    """
    bedges = get_boundary_edges_of_mesh2(faces)
    indices = np.unique(bedges.reshape(-1))
    return indices


def are_boundary_vertices2(points, faces):
    """
    Return array of 1 or 0 if vertex is boundary vertex or not
    """
    ix = get_boundary_vertices2(faces)
    areBoundaryVertices = np.zeros((len(points), 1), dtype=int)
    areBoundaryVertices[ix] = 1
    return areBoundaryVertices


def get_boundary_elements2(points, faces):
    """
    Determine the boundary elements of the mesh.
    """
    boundary_vertices = get_boundary_vertices2(faces)
    vtoe, ptr = vertex_to_elements(points, faces)
    bele = np.array([], dtype=int)
    for vertex in boundary_vertices:
        for ele in zip(vtoe[ptr[vertex] : ptr[vertex + 1]]):
            bele = np.append(bele, ele)
    bele = np.unique(bele)
    return bele


def delete_boundary_elements(points, faces, minqual=0.10):
    """
    Delete boundary elements with poor quality (i.e., < minqual)
    """
    qual = simpqual(points, faces)
    bele = get_boundary_elements2(points, faces)
    qualBou = qual[bele]
    delete = qualBou < minqual
    print(
        "Deleting " + str(np.sum(delete)) + " poor quality boundary elements...",
        flush=True,
    )
    print(faces.shape, flush=True)
    faces = np.delete(faces, delete == 1, axis=0)
    print(faces.shape, flush=True)
    points, faces, _ = fixmesh(points, faces)
    return points, faces


def collapse_edges(points, faces, minqual=0.10):
    """
    Collapse triangles that are exceedingly thin incrementally modifying the
    connectivity.
    Thin triangles are identified by containing highly actue angles.
    The shortest edge of the identified triangle is collapsed to a point.
    and the triangle table and point matrix are updated accordingly.
    """
    qual = simpqual(points, faces)
    kount = 0
    while np.any(qual < minqual):
        ixx = np.argwhere(qual < minqual)
        ix = ixx[0]
        tmpf = faces[ix, :]
        ee = np.array([tmpf[0, [0, 1]], tmpf[0, [0, 2]], tmpf[0, [1, 2]]])
        evec = points[ee[:, 0], :] - points[ee[:, 1], :]
        meid = np.argmin(np.sum(evec ** 2, axis=1))
        PIDToRemove = ee[meid, 0]  # remove this id for ReplaceID
        PIDToReplace = ee[meid, 1]  # replace remove ID with this id
        #  Delete triangle that has the thin edge
        ofaces = faces
        faces = np.delete(faces, ix, axis=0)
        faces = np.where(faces == PIDToRemove, PIDToReplace, faces)
        isOk = (
            ((faces[:, 0] - faces[:, 2]) == 0)
            + ((faces[:, 1] - faces[:, 2]) == 0)
            + ((faces[:, 0] - faces[:, 1]) == 0)
        )
        if np.any(isOk) > 0:
            # cannot delete this element, make the quality "good".
            faces = ofaces
            qual[ix] = 1
            continue
        # recompute qualities for elements with collapsed edges.
        qual = np.delete(qual, ix)
        # recompute qualities that were affected by collapse
        sel = np.argwhere(
            (faces[:, 0] == PIDToReplace)
            | (faces[:, 1] == PIDToReplace)
            | (faces[:, 2] == PIDToReplace)
        )
        qual = simpqual(points, faces)
        kount += 1
    points, faces, _ = fixmesh(points, faces)
    print("There were " + str(kount) + " thin triangles collapsed...", flush=True)
    return points, faces


def sparse(I, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...)

    Usage
    -----
    >>> shape = (m, n)
    >>> A = sparse(I, J, S, shape, dtype)
    """

    # Advanced usage: allow J and S to be scalars.
    if np.isscalar(J):
        x = J
        J = np.empty(I.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(I.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat
    II = I.flat
    J = J.flat
    return spsparse.coo_matrix((S, (II, J)), shape, dtype)


def laplacian2(points, faces, max_iter=20, tol=0.01):
    """
    Move points to the average position of their connected neighbors
    to hopefully improve triangular element quality
    """
    eps = np.finfo(float).eps

    # unused vertices
    #points = points[np.unique(faces.reshape(-1)),:]

    n = len(points)

    S = sparse(
        faces[:, [0, 0, 1, 1, 2, 2]], faces[:, [1, 2, 0, 2, 0, 1]], 1, shape=(n, n)
    )
    bnd = get_boundary_vertices2(faces)
    edge = get_edges_of_mesh2(faces)
    W = np.sum(S, 1)
    if np.any(W == 0):
        print("Invalid mesh. Disjoint vertices found. Returning", flush=True)
        print(np.argwhere(W==0),flush=True)
        return points, faces

    L = np.sqrt(
        np.sum(np.square(points[edge[:, 0], :] - points[edge[:, 1], :]), axis=1)
    )
    L[L < eps] = eps
    L = L[:, None]
    for it in range(max_iter):
        pnew = np.divide(S * np.matrix(points), np.hstack((W, W)))
        pnew[bnd, :] = points[bnd, :]
        points = pnew
        Lnew = np.sqrt(
            np.sum(np.square(points[edge[:, 0], :] - points[edge[:, 1], :]), axis=1)
        )
        Lnew[Lnew < eps] = eps
        move = np.amax(np.divide((Lnew - L), Lnew))
        if move < tol:
            print(
                "Movement tolerance reached after " + str(it) + " iterations..exiting",
                flush=True,
            )
            break
        L = Lnew
    points = np.array(points)
    return points, faces


def linter(points, faces, minqual=0.10):
    """
    Remove and check mesh for defects running a sequence of mesh improvement
    strategies.
    """
    # smooth mesh
    points, faces = laplacian2(points, faces, max_iter=30, tol=0.01)
    # delete low quality boundary elements
    points, faces = delete_boundary_elements(points, faces, minqual=minqual)
    # check for non-manifold boundaries
    bedges = get_boundary_edges_of_mesh2(faces)
    if bedges.size != points[np.unique(bedges), :].size:
        print("mesh has a non-manifold boundary...")
    # collapse thin interior triangles
    points, faces = collapse_edges(points, faces, minqual=minqual)
    # calculate final minimum simplex quality
    qual = simpqual(points, faces)
    minimum_quality = np.amin(qual)
    return points, faces, minimum_quality
