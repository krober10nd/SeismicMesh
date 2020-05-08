import scipy.sparse as spsparse
import numpy as np

from . import signed_distance_functions as sdf

"""
Routines to perform geometrical/topological operations on meshes
"""


def remove_duplicates(data):
    """
    removes duplicate rows for int data
    such as the face table
    """
    data = np.sort(data, axis=1)
    return np.unique(data, axis=1)


def remove_external_faces(points, faces, extents, dim=2):
    """
    Remove faces with all three vertices outside block (external)
    and points that are very far from local domain
    """
    if dim == 2:
        signed_distance = sdf.drectangle(
            points[faces.flatten(), :],
            x1=extents[0],
            x2=extents[2],
            y1=extents[1],
            y2=extents[3],
        )
    elif dim == 3:
        signed_distance = sdf.dblock(
            points[faces.flatten(), :],
            x1=extents[0],
            x2=extents[3],
            y1=extents[1],
            y2=extents[4],
            z1=extents[2],
            z2=extents[5],
        )
    # keep faces that don't have all their nodes "out" of the local domain
    # and
    # faces that have all their nodes in the "close" to interior of the local block
    isOut = np.reshape(signed_distance > 0, (-1, (dim + 1)))
    # todo: this needs to be more objective
    isFar = np.reshape(signed_distance > 1000, (-1, (dim + 1)))
    faces_new = faces[
        (np.sum(isOut, axis=1) != (dim + 1)) & (np.any(isFar, axis=1) != 1), :
    ]
    points_new, faces_new, jx = fixmesh(points, faces_new, dim=dim)
    return points_new, faces_new, jx


def vertex_to_elements(points, faces, dim=2):
    """
    Determine which elements connected to vertices
    """
    num_faces = len(faces)

    ext = np.tile(np.arange(0, num_faces), (dim + 1, 1)).reshape(-1, order="F")
    ve = np.reshape(faces, (-1,))
    ve = np.vstack((ve, ext)).T
    ve = ve[ve[:, 0].argsort(), :]

    idx = np.insert(np.diff(ve[:, 0]), 0, 0)
    vtoe_pointer = np.argwhere(idx)
    vtoe_pointer = np.insert(vtoe_pointer, 0, 0)
    vtoe_pointer = np.append(vtoe_pointer, num_faces * (dim + 1))

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
    elif dim == 3:
        d01 = p[t[:, 1], :] - p[t[:, 0], :]
        d02 = p[t[:, 2], :] - p[t[:, 0], :]
        d03 = p[t[:, 3], :] - p[t[:, 0], :]
        return np.einsum("ij,ij->i", np.cross(d01, d02), d03) / 6
    else:
        raise NotImplementedError


def fixmesh(p, t, ptol=2e-13, delunused=False, dim=2):
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
    if delunused:
        pix, _, jx = np.unique(t, return_index=True, return_inverse=True)
        t = np.reshape(jx, (t.shape))
        p = p[pix, :]
        pix = ix[pix]

    # element orientation is CCW
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


def get_edges(faces, dim=2):
    """
    Get the edges of mesh in no order
    and are repeated.
    """
    num_faces = len(faces)
    faces = np.array(faces)
    if dim == 2:
        edges = faces[:, [[0, 1], [0, 2], [1, 2]]]
        edges = edges.reshape((num_faces * 3, 2))
    elif dim == 3:
        edges = faces[:, [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]]]
        edges = edges.reshape((num_faces * 6, 2))
    return edges


def get_boundary_edges(faces):
    """
    Get the boundary edges of the mesh.
    """
    edges = get_edges(faces)
    edges = np.sort(edges, axis=1)
    unq, cnt = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_edges


def get_winded_boundary_edges(faces):
    """
    Order the boundary edges of the mesh in a winding fashion.
    only works in 2D
    """
    boundary_edges = get_boundary_edges(faces)
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


def get_boundary_vertices(faces, dim=2):
    """
    Get the indices of the mesh representing boundary vertices.
    works in 2d and 3d
    """
    if dim == 2:
        b = get_boundary_edges(faces)
    elif dim == 3:
        b = get_boundary_facets(faces)
    indices = np.unique(b.reshape(-1))
    return indices


def get_boundary_elements(points, faces):
    """
    Determine the boundary elements of the mesh.
    """
    boundary_vertices = get_boundary_vertices(faces)
    vtoe, ptr = vertex_to_elements(points, faces)
    bele = np.array([], dtype=int)
    for vertex in boundary_vertices:
        for ele in zip(vtoe[ptr[vertex] : ptr[vertex + 1]]):
            bele = np.append(bele, ele)
    bele = np.unique(bele)
    return bele


def get_facets(cells):
    """
    Gets the 4 facets of each cell in no order
    """
    ix = [[0, 1, 3], [1, 2, 3], [2, 0, 3], [1, 2, 0]]
    return np.array(cells[:, ix]).reshape((len(cells) * 4, 3))


def get_boundary_facets(cells):
    """
    Get the facets shared by owned 1 cell
    """
    facets = get_facets(cells)
    facets = np.sort(facets, axis=1)
    unq, cnt = np.unique(facets, axis=0, return_counts=True)
    boundary_facets = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_facets


def delete_boundary_elements(points, faces, minqual=0.10, dim=2):
    """
    Delete boundary elements with poor quality (i.e., < minqual)
    """
    qual = simpqual(points, faces)
    bele = get_boundary_elements(points, faces)
    qualBou = qual[bele]
    delete = qualBou < minqual
    print(
        "Deleting " + str(np.sum(delete)) + " poor quality boundary elements...",
        flush=True,
    )
    delete = np.argwhere(delete == 1)
    faces = np.delete(faces, bele[delete], axis=0)
    points, faces, _ = fixmesh(points, faces, delunused=True, dim=dim)
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

    n = len(points)

    S = sparse(
        faces[:, [0, 0, 1, 1, 2, 2]], faces[:, [1, 2, 0, 2, 0, 1]], 1, shape=(n, n)
    )
    bnd = get_boundary_vertices(faces)
    edge = get_edges(faces)
    W = np.sum(S, 1)
    if np.any(W == 0):
        print("Invalid mesh. Disjoint vertices found. Returning", flush=True)
        print(np.argwhere(W == 0), flush=True)
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


def isManifold(points, faces, dim=2):
    """
    Determine if mesh is manifold.
    1. A boundary edge should have one element
    2. A non-boundary edge should have two elements
    3. The number of boundary vertices == number of boundary edges
    """
    bedges = get_boundary_edges(faces, dim=dim)
    if bedges.size != points[np.unique(bedges), :].size:
        print("Mesh has a non-manifold boundary...", flush=True)
        return False
    return True


def ptInFace2(point, face):
    """
    Does the 2D point lie in the face with vertices (x1,y1,x2,y2,x3,y3)
    2D face?
    """
    (x, y) = point
    (x1, y1, x2, y2, x3, y3) = face
    a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / (
        (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    )
    b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / (
        (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    )
    c = 1 - a - b
    # pt lies in T if and only if 0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1
    return 0 <= a and a <= 1 and 0 <= b and b <= 1 and 0 <= c and c <= 1


def getCentroids(points, faces, dim=2):
    """
    Calculate the centroids of all the faces
    """
    return points[faces].sum(1) / (dim + 1)


def doAnyFacesOverlap(points, faces):
    """
    Check if any faces connected to boundary of the mesh overlap
    with another face connected to the boundary ignoring self-intersections.
    Checks only the 1-ring around each boundary element for potential intersections
    """
    vtoe, ptr = vertex_to_elements(points, faces)
    # all elements that have a boundary vertex
    beles = get_boundary_elements(points, faces)
    # centroids of these elements beles
    bcents = getCentroids(points, faces[beles, :])
    # store intersection pairs
    intersections = []
    # for all boundary triangles
    for ie, cent in enumerate(bcents):
        # collect all elements neis around boundary element ie
        neis = np.array([], dtype=int)
        for vertex in faces[beles[ie], :]:
            for ele in zip(vtoe[ptr[vertex] : ptr[vertex + 1]]):
                neis = np.append(neis, ele)
        # for all neighboring elements  to boundary element ie
        for iee, bele in enumerate(neis):
            # centroid ie should be in face ie by definition!
            if beles[ie] == bele:
                continue
            # does a centroid live inside another face?
            x1, x2, x3 = points[faces[bele, :], 0]
            y1, y2, y3 = points[faces[bele, :], 1]
            if ptInFace2((cent[0], cent[1]), (x1, y1, x2, y2, x3, y3)):
                # centroid ie is inside face iee
                print(
                    "Alert: face "
                    + str(beles[ie])
                    + " intersects with face "
                    + str(bele)
                    + ". These will be adjusted.",
                    flush=True,
                )
                # record all intersection pairs
                intersections.append((beles[ie], bele))
    return intersections


def linter(points, faces, minqual=0.10):
    """
    Remove and check mesh for defects
    """
    qual = simpqual(points, faces)
    # determine if there's degenerate overlapping elements
    intersections = doAnyFacesOverlap(points, faces)
    # delete the lower quality in the pair
    delete = []
    for intersect in intersections:
        ix = [i for i in intersect]
        sel = np.argmin(qual[ix])
        delete = np.append(delete, intersect[sel])
    delete = np.unique(delete)
    print("Deleting " + str(len(delete)) + " overlapped faces", flush=True)
    faces = np.delete(faces, delete, axis=0)
    # clean up
    points, faces, _ = fixmesh(points, faces, delunused=True)
    # delete remaining low quality boundary elements
    points, faces = delete_boundary_elements(points, faces, minqual=minqual)
    # check for non-manifold boundaries and alert
    _ = isManifold(points, faces)
    # calculate final minimum simplex quality
    qual = simpqual(points, faces)
    minimum_quality = np.amin(qual)
    print(
        "There are "
        + str(len(points))
        + " vertices and "
        + str(len(faces))
        + " elements in the mesh",
        flush=True,
    )
    print("The minimum element quality is " + str(minimum_quality), flush=True)
    return points, faces
