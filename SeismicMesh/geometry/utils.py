import scipy.sparse as spsparse
import numpy as np

from . import signed_distance_functions as sdf

from .cpp import fast_geometry as gutils

from ..generation.cpp import c_cgal

"""
Routines to perform geometrical/topological operations on meshes
"""

# cpp implementation of 4x4 determinant calc
dete = gutils.calc_4x4determinant


def calc_re_ratios(points, cells):
    """
    Calculate radius edge ratios--mesh quality metric
    """
    # circumradius/shortest edge length
    bars = np.concatenate(
        [
            cells[:, [0, 1]],
            cells[:, [1, 2]],
            cells[:, [2, 0]],
            cells[:, [0, 3]],
            cells[:, [1, 3]],
            cells[:, [2, 3]],
        ]
    )
    barvec = points[bars[:, 0]] - points[bars[:, 1]]
    L = np.sqrt((barvec ** 2).sum(1))
    L = np.reshape(L, (6, len(cells)))
    # min edge length i every tetra
    minL = np.amin(L, axis=0)
    cc = c_cgal.circumballs3(points[cells.flatten()])
    r = cc[:, -1]
    return np.sqrt(r) / minL


def dump_mesh(points, cells, rank):
    import meshio

    meshio.write_points_cells(
        "test_" + str(rank) + ".vtk", points, [("tetra", cells)],
    )
    np.savetxt("points.txt", points, delimiter=",")
    np.savetxt("cells.txt", cells, delimiter=",")


def remove_duplicates(data):
    """
    removes duplicate rows for int data
    such as the face table
    """
    data = np.sort(data, axis=1)
    return np.unique(data, axis=1)


def remove_external_faces(points, faces, extent, new_idx, dim=2):
    """
    Remove entities with all dim+1 vertices outside block (external)
    and points that are "far" from local domain extents.
    """
    if dim == 2:
        signed_distance = sdf.drectangle(
            points[faces.flatten(), :],
            x1=extent[0],
            x2=extent[2],
            y1=extent[1],
            y2=extent[3],
        )
    elif dim == 3:
        signed_distance = sdf.dblock(
            points[faces.flatten(), :],
            x1=extent[0],
            x2=extent[3],
            y1=extent[1],
            y2=extent[4],
            z1=extent[2],
            z2=extent[5],
        )
    # keep faces that don't have all their nodes "out" of the local domain
    # and
    # faces that have all their nodes in the "close" to interior of the local block
    # isOut = np.reshape(signed_distance > 0, (-1, (dim + 1)))
    # determine minimum extent
    mIext = 99999999.0
    for ix in range(0, dim):
        mIext = np.amin(
            [mIext, extent[int(dim * 2) - int(ix) - 1] - extent[dim - int(ix) - 1]]
        )
    isOut = faces >= new_idx
    # if greater than x times the minimum lengthscale of the subdomain
    isFar = np.reshape(signed_distance > 0.50 * mIext, (-1, (dim + 1)))
    # high aspect ratio tetrahedrals sometimes occur on the boundary delete these
    faces_new = faces[
        (np.sum(isOut, axis=1) != (dim + 1)), :
    ]  # and (np.any(isFar, axis=1) != 1), :
    # ]
    points_new, faces_new, jx = fixmesh(points, faces_new, dim=dim)
    return points_new, faces_new, jx


def vertex_to_elements(points, faces, dim=2):
    """
    Determine which elements are connected to vertices
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


def fixmesh(p, t, ptol=2e-13, delunused=False, delslivers=False, dim=2, delete=None):
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
    if delete is not None:
        t = np.delete(t, delete, axis=0)

    # duplicate vertices
    snap = (p.max(0) - p.min(0)).max() * ptol
    _, ix, jx = unique_rows(np.round(p / snap) * snap, True, True)

    p = p[ix]
    t = jx[t]

    # duplicate elements
    t = np.sort(t, axis=1)
    t = unique_rows(t)

    # delete slivers
    if delslivers:
        dh_angles = np.rad2deg(gutils.calc_dihedral_angles(p, t))
        outOfBounds = np.argwhere((dh_angles[:, 0] < 5) | (dh_angles[:, 0] > 175))
        eleNums = np.floor(outOfBounds / 6).astype("int")
        eleNums, ix = np.unique(eleNums, return_index=True)
        print("Deleting " + str(len(eleNums)) + " slivers...", flush=True)
        t = np.delete(t, eleNums, axis=0)

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
    radius-to-edge ratio
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


def get_boundary_edges(faces, dim=2):
    """
    Get the boundary edges of the mesh.
    """
    edges = get_edges(faces, dim=dim)
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


def get_boundary_elements(points, faces, dim=2):
    """
    Determine the boundary elements of the mesh.
    """
    boundary_vertices = get_boundary_vertices(faces, dim=dim)
    vtoe, ptr = vertex_to_elements(points, faces, dim=dim)
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
    Get the facets shared by only 1 cell
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
    bele = get_boundary_elements(points, faces, dim=dim)
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


def sparse(Ix, J, S, shape=None, dtype=None):
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
        J = np.empty(Ix.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(Ix.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat
    II = Ix.flat
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


def ptInCell3(point, cell):
    """
    Does the 3D point lie in the face with vertices (x1,y1,z1,x2,y2,z2,x3,y3,z3)
    3D cell?
    """
    # utilize the fixed size determinant in cpp

    (x, y, z) = point
    (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4) = cell
    signs = np.zeros(5, dtype=int)
    A = np.array(
        [[x1, y1, z1, 1.0], [x2, y2, z2, 1.0], [x3, y3, z3, 1.0], [x4, y4, z4, 1.0]],
    )
    signs[0] = np.sign(dete(A))
    B = np.array(
        [[x, y, z, 1.0], [x2, y2, z2, 1.0], [x3, y3, z3, 1.0], [x4, y4, z4, 1.0]],
    )
    signs[1] = np.sign(dete(B))
    if signs[1] != signs[0]:
        return False
    C = np.array(
        [[x1, y1, z1, 1.0], [x, y, z, 1.0], [x3, y3, z3, 1.0], [x4, y4, z4, 1.0]],
    )
    signs[2] = np.sign(dete(C))
    if signs[2] != signs[1]:
        return False
    D = np.array(
        [[x1, y1, z1, 1.0], [x2, y2, z2, 1.0], [x, y, z, 1.0], [x4, y4, z4, 1.0]],
    )
    signs[3] = np.sign(dete(D))
    if signs[3] != signs[2]:
        return False
    E = np.array(
        [[x1, y1, z1, 1.0], [x2, y2, z2, 1.0], [x3, y3, z3, 1.0], [x, y, z, 1.0]],
    )
    signs[4] = np.sign(dete(E))
    if signs[4] != signs[3]:
        return False
    # point lies in cell if and only if all signs are the same
    return (signs == -1).all() or (signs == 1).all()


def get_centroids(points, entities, dim=2):
    """
    Calculate the centroids of all the faces
    """
    return points[entities].sum(1) / (dim + 1)


def doAnyOverlap(points, entities, dim=2):
    """
    Check if any entities of dim D connected to boundary of the mesh overlap
    with another entity D connected to the boundary ignoring self-intersections.
    Checks only the 1-ring around each entity for potential intersections
    using barycentric coordinates.
    """

    vtoe, ptr = vertex_to_elements(points, entities, dim=dim)
    # centroids of these elements beles
    cents = get_centroids(points, entities, dim=dim)
    # store intersection pairs
    intersections = []
    # for all elements
    for ie, cent in enumerate(cents):
        if ie % 1000 == 0:
            print("INFO: " + str(100.0 * (ie / len(cents))) + " % done...", flush=True)
        # collect all elements neis around element ie
        neis = np.array([], dtype=int)
        for vertex in entities[ie, :]:
            for ele in zip(vtoe[ptr[vertex] : ptr[vertex + 1]]):
                neis = np.append(neis, ele)
        # for all neighboring elements  to element ie
        for iee, ele in enumerate(neis):
            # centroid ie should be in face ie by definition!
            if ie == ele:
                continue
            if dim == 2:
                # does this centroid live inside another neighboring face?
                x1, x2, x3 = points[entities[ele, :], 0]
                y1, y2, y3 = points[entities[ele, :], 1]
                if ptInFace2((cent[0], cent[1]), (x1, y1, x2, y2, x3, y3)):
                    # centroid ie is inside face iee
                    print(
                        "Alert: face "
                        + str(ie)
                        + " intersects with face "
                        + str(ele)
                        + ". These will be adjusted.",
                        flush=True,
                    )
                    # record all intersection pairs
                    intersections.append((ie, ele))
            elif dim == 3:
                # does this centroid live inside another nei cell?
                x1, x2, x3, x4 = points[entities[ele], 0]
                y1, y2, y3, y4 = points[entities[ele], 1]
                z1, z2, z3, z4 = points[entities[ele], 2]
                if ptInCell3(
                    (cent[0], cent[1], cent[2]),
                    (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4),
                ):
                    # centroid ie is inside cell iee
                    print(
                        "Alert: cell "
                        + str(ie)
                        + " intersects with cell "
                        + str(ele)
                        + ". These will be adjusted.",
                        flush=True,
                    )
                    # record all intersection pairs
                    intersections.append((ie, ele))

    return intersections


def linter(points, faces, minqual=0.10, dim=2):
    """
    Remove and check mesh for defects
    """
    print("Performing mesh linting...", flush=True)
    qual = simpqual(points, faces)
    # determine if there's degenerate overlapping elements
    import time

    t1 = time.time()
    intersections = doAnyOverlap(points, faces, dim=dim)
    # delete the lower quality in the pair
    delete = []
    for intersect in intersections:
        ix = [i for i in intersect]
        sel = np.argmin(qual[ix])
        delete = np.append(delete, intersect[sel])
    delete = np.unique(delete)
    print("Deleting " + str(len(delete)) + " overlapped faces", flush=True)
    faces = np.delete(faces, delete, axis=0)
    print(time.time() - t1)

    # clean up
    points, faces, _ = fixmesh(points, faces, delunused=True)
    # delete remaining low quality boundary elements
    if dim == 2:
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
