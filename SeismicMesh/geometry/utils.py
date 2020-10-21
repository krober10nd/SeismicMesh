import numpy as np
import scipy.sparse as spsparse

from ..generation.cpp import c_cgal
from . import signed_distance_functions as sdf
from .cpp import fast_geometry as gutils

# cpp implementation of 4x4 determinant calc
dete = gutils.calc_4x4determinant


def calc_re_ratios(vertices, entities, dim=2):
    """Calculate radius edge ratios--mesh quality metric

    :param vertices: point coordinates of the mesh vertices'
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity table
    :type entities: numpy.ndarray[`int` x (dim + 1)]

    :return: ce_ratios: array of radius-to-edge ratios
    :rtype: ce_ratios: numpy.ndarray[`float` x 1]
    """
    # circumradius/shortest edge length
    if dim == 2:
        bars = np.concatenate(
            [entities[:, [0, 1]], entities[:, [1, 2]], entities[:, [2, 0]]]
        )
        cc = c_cgal.circumballs2(vertices[entities.flatten()])
    elif dim == 3:
        bars = np.concatenate(
            [
                entities[:, [0, 1]],
                entities[:, [1, 2]],
                entities[:, [2, 0]],
                entities[:, [0, 3]],
                entities[:, [1, 3]],
                entities[:, [2, 3]],
            ]
        )
        cc = c_cgal.circumballs3(vertices[entities.flatten()])
    else:
        raise ValueError("Dimension invalid")
    r = cc[:, -1]
    bar_vec = vertices[bars[:, 0]] - vertices[bars[:, 1]]
    L = np.sqrt((bar_vec ** 2).sum(1))
    L = np.reshape(L, (3 * (dim - 1), len(entities)))
    # min edge length i every tetra
    minL = np.amin(L, axis=0)
    return np.sqrt(r) / minL


def remove_external_entities(vertices, entities, extent, dim=2):
    """Remove entities with all dim+1 vertices outside block.

    :param vertices: point coordinates of mesh
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim + 1)]
    :param extent: coords. of the local block extents
    :type extent: numpy.ndarray[tuple(`float` x (2*dim))]
    :param dim: dimension of mesh
    :type dim: `int`, optional

    :return: vertices_new: point coordinates of mesh w/ removed entities
    :rtype: numpy.ndarray[`float` x dim]
    :return: entities_new: mesh connectivity w/ removed entities
    :rtype: numpy.ndarray[`int` x (dim +1)]
    :return: jx: mapping from old point indexing to new point indexing
    :rtype: numpy.ndarray[`int` x 1]
    """

    if dim == 2:
        signed_distance = sdf.drectangle(
            vertices[entities.flatten(), :],
            x1=extent[0],
            x2=extent[2],
            y1=extent[1],
            y2=extent[3],
        )
    elif dim == 3:
        signed_distance = sdf.dblock(
            vertices[entities.flatten(), :],
            x1=extent[0],
            x2=extent[3],
            y1=extent[1],
            y2=extent[4],
            z1=extent[2],
            z2=extent[5],
        )
    isOut = np.reshape(signed_distance > 0.0, (-1, (dim + 1)))
    entities_new = entities[(np.sum(isOut, axis=1) != (dim + 1))]
    vertices_new, entities_new, jx = fix_mesh(vertices, entities_new, dim=dim)
    return vertices_new, entities_new, jx


def vertex_to_entities(vertices, entities, dim=2):
    """Determine which elements are connected to which vertices.

    :param vertices: point coordinates of mesh vertices
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim + 1)]
    :param dim: dimension of mesh
    :type dim: `int`, optional

    :return: vtoe: indices of entities connected to each vertex
    :rtype: numpy.ndarray[`int` x 1]

    :return: vtoe_pointer: indices into `vtoe` such that vertex `v` is connected to
                          `vtoe[vtoe_pointer[v]:vtoe_pointer[v+1]]` entities
    :rtype: numpy.ndarray[`int` x 1]
    """
    num_entities = len(entities)

    ext = np.tile(np.arange(0, num_entities), (dim + 1, 1)).reshape(-1, order="F")
    ve = np.reshape(entities, (-1,))
    ve = np.vstack((ve, ext)).T
    ve = ve[ve[:, 0].argsort(), :]

    idx = np.insert(np.diff(ve[:, 0]), 0, 0)
    vtoe_pointer = np.argwhere(idx)
    vtoe_pointer = np.insert(vtoe_pointer, 0, 0)
    vtoe_pointer = np.append(vtoe_pointer, num_entities * (dim + 1))

    vtoe = ve[:, 1]

    return vtoe, vtoe_pointer


def unique_rows(A, return_index=False, return_inverse=False):
    """Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    :param  A: array of data
    :type A: numpy.ndarray[`int`/`float` x N]
    :param return_index: whether to return the indices of unique data
    :type return_index: `boolean`, optional
    :param return_inverse: whether to return the inverse mapping back to A from B.
    :type return_inverse: `boolean`, optional

    :return: B: array of data with duplicates removed
    :rtype: numpy.ndarray[`int`/`float` x N]
    :return: I: array of indices to unique data B.
    :rtype: numpy.ndarray[`int` x 1]
    :return: J: array of indices to A from B.
    :rtype: numpy.ndarray[`int` x 1]
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


def simp_vol(p, t):
    """Signed volumes of the simplex elements in the mesh.

    :param p: point coordinates of mesh
    :type p: numpy.ndarray[`float` x dim]
    :param t: mesh connectivity
    :type t: numpy.ndarray[`int` x (dim + 1)]

    :return: volume: signed volume of entity/simplex.
    :rtype: numpy.ndarray[`float` x 1]
    """

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


def fix_mesh(p, t, ptol=2e-13, dim=2, delete_unused=False):
    """Remove duplicated/unused vertices and entities and
       ensure orientation of entities is CCW.

    :param p: point coordinates of mesh
    :type p: numpy.ndarray[`float` x dim]
    :param t: mesh connectivity
    :type t: numpy.ndarray[`int` x (dim + 1)]
    :param ptol: point tolerance to detect duplicates
    :type ptol: `float`, optional
    :param dim: dimension of mesh
    :type dim: `int`, optional
    :param delete_unused: flag to delete disjoint vertices.
    :type delete_unused: `boolean`, optional

    :return: p: updated point coordinates of mesh
    :rtype: numpy.ndarray[`float` x dim]
    :return: t: updated mesh connectivity
    :rtype: numpy.ndarray[`int` x (dim+1)]
    """

    # duplicate vertices
    snap = (p.max(0) - p.min(0)).max() * ptol
    _, ix, jx = unique_rows(np.round(p / snap) * snap, True, True)

    p = p[ix]
    t = jx[t]

    # duplicate entities
    t = np.sort(t, axis=1)
    t = unique_rows(t)

    # delete disjoint vertices
    if delete_unused:
        pix, _, jx = np.unique(t, return_index=True, return_inverse=True)
        t = np.reshape(jx, (t.shape))
        p = p[pix, :]

    # entity orientation is CCW
    flip = simp_vol(p, t) < 0
    t[flip, :2] = t[flip, 1::-1]

    return p, t, jx


def simp_qual(p, t):
    """Simplex quality radius-to-edge ratio

    :param p: vertex coordinates of mesh
    :type p: numpy.ndarray[`float` x dim]
    :param t: mesh connectivity
    :type t: numpy.ndarray[`int` x (dim + 1)]

    :return: signed mesh quality: signed mesh quality (1.0 is perfect)
    :rtype: numpy.ndarray[`float` x 1]

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


def get_edges(entities, dim=2):
    """Get the undirected edges of mesh in no order (NB: are repeated)

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional

    :return: edges: the edges that make up the mesh
    :rtype: numpy.ndarray[`int`x 2]
    """

    num_entities = len(entities)
    entities = np.array(entities)
    if dim == 2:
        edges = entities[:, [[0, 1], [0, 2], [1, 2]]]
        edges = edges.reshape((num_entities * 3, 2))
    elif dim == 3:
        edges = entities[:, [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]]]
        edges = edges.reshape((num_entities * 6, 2))
    return edges


def get_boundary_edges(entities, dim=2):
    """Get the boundary edges of the mesh. Boundary edges only appear (dim-1) times

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional

    :return: boundary_edges: the edges that make up the boundary of the mesh
    :rtype: numpy.ndarray[`int` x 2]
    """
    edges = get_edges(entities, dim=dim)
    edges = np.sort(edges, axis=1)
    unq, cnt = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = np.array([e for e, c in zip(unq, cnt) if c == (dim - 1)])
    return boundary_edges


def get_winded_boundary_edges(entities):
    """Order the boundary edges of the mesh in a winding fashiono

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]

    :return: boundary_edges: the edges that make up the boundary of the mesh in a winding order
    :rtype: numpy.ndarray[`int` x 2]
    """

    boundary_edges = get_boundary_edges(entities)
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


def get_boundary_vertices(entities, dim=2):
    """Get the indices of the mesh representing boundary vertices.

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int `x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional

    :return: indices: indices into the vertex array that are on the boundary.
    :rtype: numpy.ndarray[`float` x dim]
    """
    if dim == 2:
        b = get_boundary_edges(entities)
    elif dim == 3:
        b = get_boundary_facets(entities)
    else:
        raise ValueError("Dimension not supported.")
    indices = np.unique(b.reshape(-1))
    return indices


def get_boundary_entities(vertices, entities, dim=2):
    """Determine the entities that lie on the boundary of the mesh.

    :param vertices: vertex coordinates of mesh
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional

    :return: bele: indices of entities on the boundary of the mesh.
    :rtype: numpy.ndarray[`int` x 1]
    """
    boundary_vertices = get_boundary_vertices(entities, dim=dim)
    vtoe, ptr = vertex_to_entities(vertices, entities, dim=dim)
    bele = np.array([], dtype=int)
    for vertex in boundary_vertices:
        for ele in zip(vtoe[ptr[vertex] : ptr[vertex + 1]]):
            bele = np.append(bele, ele)
    bele = np.unique(bele)
    return bele


def get_facets(entities):
    """Gets the four facets of each tetrahedral.

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]

    :return: facets: facets of a tetrahedral entity.
    :rtype: numpy.ndarray[`int` x 4]
    """
    ix = [[0, 1, 3], [1, 2, 3], [2, 0, 3], [1, 2, 0]]
    return np.array(entities[:, ix]).reshape((len(entities) * 4, 3))


def get_boundary_facets(entities):
    """Get the facets that represent the boundary of a 3D mesh.

    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]

    :return: boundary_facets: facets on the boundary of a 3D mesh.
    :rtype: numpy.ndarray[`int` x 4]
    """
    if entities.shape[1] < 4:
        raise ValueError("Only works for triangles")
    facets = get_facets(entities)
    facets = np.sort(facets, axis=1)
    unq, cnt = np.unique(facets, axis=0, return_counts=True)
    boundary_facets = np.array([e for e, c in zip(unq, cnt) if c == 1])
    return boundary_facets


def delete_boundary_entities(vertices, entities, dim=2, min_qual=0.10):
    """Delete boundary entities with poor geometric quality (i.e., < min. quality)

    :param vertices: vertex coordinates of mesh
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional
    :param min_qual: minimum geometric quality to consider "poor" quality
    :type min_qual: `float`, optional

    :return: vertices: updated vertex array of mesh
    :rtype: numpy.ndarray[`int` x dim]
    :return: entities: update mesh connectivity
    :rtype: numpy.ndarray[`int` x (dim+1)]
    """
    qual = simp_qual(vertices, entities)
    bele = get_boundary_entities(vertices, entities, dim=dim)
    qualBou = qual[bele]
    delete = qualBou < min_qual
    print(
        "Deleting " + str(np.sum(delete)) + " poor quality boundary entities...",
        flush=True,
    )
    delete = np.argwhere(delete == 1)
    entities = np.delete(entities, bele[delete], axis=0)
    vertices, entities, _ = fix_mesh(vertices, entities, delete_unused=True, dim=dim)
    return vertices, entities


def _sparse(Ix, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...)
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


def laplacian2(vertices, entities, max_iter=20, tol=0.01):
    """Move vertices to the average position of their connected neighbors
    with the goal to hopefully improve geometric entity quality.

    :param vertices: vertex coordinates of mesh
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param max_iter: maximum number of iterations to perform
    :type max_iter: `int`, optional
    :param tol: iterations will cease when movement < tol
    :type tol: `float`, optional

    :return vertices: updated vertices of mesh
    :rtype: numpy.ndarray[`float` x dim]
    :return: entities: updated mesh connectivity
    :rtype: numpy.ndarray[`int` x (dim+1)]
    """
    if vertices.ndim != 2:
        raise NotImplementedError("Laplacian smoothing is only works in 2D for now")

    eps = np.finfo(float).eps

    n = len(vertices)

    S = _sparse(
        entities[:, [0, 0, 1, 1, 2, 2]],
        entities[:, [1, 2, 0, 2, 0, 1]],
        1,
        shape=(n, n),
    )
    bnd = get_boundary_vertices(entities)
    edge = get_edges(entities)
    W = np.sum(S, 1)
    if np.any(W == 0):
        print("Invalid mesh. Disjoint vertices found. Returning", flush=True)
        print(np.argwhere(W == 0), flush=True)
        return vertices, entities

    L = np.sqrt(
        np.sum(np.square(vertices[edge[:, 0], :] - vertices[edge[:, 1], :]), axis=1)
    )
    L[L < eps] = eps
    L = L[:, None]
    for it in range(max_iter):
        pnew = np.divide(S * np.matrix(vertices), np.hstack((W, W)))
        pnew[bnd, :] = vertices[bnd, :]
        vertices = pnew
        Lnew = np.sqrt(
            np.sum(np.square(vertices[edge[:, 0], :] - vertices[edge[:, 1], :]), axis=1)
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
    vertices = np.array(vertices)
    return vertices, entities


def is_manifold(vertices, entities, dim=2):
    """Determine if mesh is manifold by checking for the following:
    1. A boundary edge should be a member of one entity
    2. A non-boundary edge should be a member of two entities
    3. The number of boundary vertices == number of boundary edges

    :param vertices: vertex coordinates of mesh
    :type vertices: numpy.ndarray[`float` x dim]
    :param entities: the mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of the mesh
    :type dim: `int`, optional

    :return: is_manifold: flag to indicate if the mesh has a manifold boundary.
    :rtype: `boolean`.
    """
    bedges = get_boundary_edges(entities, dim=dim)
    if bedges.size != vertices[np.unique(bedges), :].size:
        print("Mesh has a non-manifold boundary...", flush=True)
        return False
    return True


def vertex_in_entity2(vertex, entity):
    """
    Does the 2D vertex lie in the entity defined by vertices (x1,y1,x2,y2,x3,y3)?

    :param vertex: vertex coordinates of mesh
    :type vertex: numpy.ndarray[`float` x dim]
    :param entity: connectivity of an entity
    :type entity: numpy.ndarray[`int` x (dim+1)]

    :return: vertex_in_entity2: logical flag indicating if it is or isn't.
    :rtype: `boolean`
    """
    (x, y) = vertex
    (x1, y1, x2, y2, x3, y3) = entity
    a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / (
        (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    )
    b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / (
        (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    )
    c = 1 - a - b
    # vt lies in entity if and only if 0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1
    return 0 <= a and a <= 1 and 0 <= b and b <= 1 and 0 <= c and c <= 1


def vertex_in_entity3(vertex, entity):
    """
    Does the 3D vertex lie in the entity defined by vertices (x1,y1,z1,x2,y2,z2,x3,y3,z3)?

    :param vertex: vertex coordinates of mesh
    :type vertex: numpy.ndarray[`float` x dim]
    :param entity: connectivity of an entity
    :type entity: numpy.ndarray[`int` x (dim+1)]

    :return: vertex_in_entity3: logical flag indicating if it is or isn't.
    :rtype: `boolean`
    """
    (x, y, z) = vertex
    (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4) = entity
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


def get_centroids(vertices, entities, dim=2):
    """Calculate the centroids of all the entities.

    :param vertex: vertex coordinates of mesh
    :type vertex: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of mesh
    :type dim: `int`, optional

    :return: centroids of entities
    :rtype: numpy.ndarray[`float` x dim]
    """
    return vertices[entities].sum(1) / (dim + 1)


def do_any_overlap(vertices, entities, dim=2):
    """
    Check if any entities connected to boundary of the mesh overlap
    ignoring self-intersections. This routine checks only the 1-ring around
    each entity for potential intersections using barycentric coordinates.

    :param vertex: vertex coordinates of mesh
    :type vertex: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of mesh
    :type dim: `int`, optional

    :return: intersections: a list of 2-tuple of entity indices that intersect
    :rtype: List[tuple(num_intersections x 2)]
    """
    vtoe, ptr = vertex_to_entities(vertices, entities, dim=dim)
    # centroids of these elements beles
    cents = get_centroids(vertices, entities, dim=dim)
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
                x1, x2, x3 = vertices[entities[ele, :], 0]
                y1, y2, y3 = vertices[entities[ele, :], 1]
                if vertex_in_entity2((cent[0], cent[1]), (x1, y1, x2, y2, x3, y3)):
                    # centroid ie is inside face iee
                    print(
                        "Alert: entity "
                        + str(ie)
                        + " intersects with entity "
                        + str(ele)
                        + ". These will be adjusted.",
                        flush=True,
                    )
                    # record all intersection pairs
                    intersections.append((ie, ele))
            elif dim == 3:
                # does this centroid live inside another nei cell?
                x1, x2, x3, x4 = vertices[entities[ele], 0]
                y1, y2, y3, y4 = vertices[entities[ele], 1]
                z1, z2, z3, z4 = vertices[entities[ele], 2]
                if vertex_in_entity3(
                    (cent[0], cent[1], cent[2]),
                    (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4),
                ):
                    # centroid ie is inside cell iee
                    print(
                        "Alert: entity "
                        + str(ie)
                        + " intersects with entity "
                        + str(ele)
                        + ". These will be adjusted.",
                        flush=True,
                    )
                    # record all intersection pairs
                    intersections.append((ie, ele))

    return intersections


def linter(vertices, entities, dim=2, min_qual=0.10):
    """Remove and check mesh for geometric and toplogical defects.

    :param vertex: vertex coordinates of mesh
    :type vertex: numpy.ndarray[`float` x dim]
    :param entities: mesh connectivity
    :type entities: numpy.ndarray[`int` x (dim+1)]
    :param dim: dimension of mesh
    :type dim: `int`, optional
    :param min_qual: minimum geometric quality to consider "poor" quality
    :type min_qual: `float`, optional

    :return vertices: updated mesh vertices
    :rtype: numpy.ndarray[`float` x dim]
    :return: entities: updated mesh connectivity
    :rtype: numpy.ndarray[`int` x (dim+1)]
    """
    print("Performing mesh linting...", flush=True)
    qual = simp_qual(vertices, entities)
    # determine if there's degenerate overlapping elements
    intersections = do_any_overlap(vertices, entities, dim=dim)
    # delete the lower quality in the pair
    delete = np.array([], dtype=int)
    for intersect in intersections:
        ix = [i for i in intersect]
        sel = np.argmin(qual[ix])
        delete = np.append(delete, intersect[sel])
    delete = np.unique(delete)
    print("Deleting " + str(len(delete)) + " overlapped entities", flush=True)
    entities = np.delete(entities, delete, axis=0)

    # clean up
    vertices, entities, _ = fix_mesh(vertices, entities, delete_unused=True)
    # delete remaining low quality boundary elements
    if dim == 2:
        vertices, entities = delete_boundary_entities(
            vertices, entities, min_qual=min_qual
        )
        # check for non-manifold boundaries and alert
        _ = is_manifold(vertices, entities)
    # calculate final minimum simplex quality
    qual = simp_qual(vertices, entities)
    minimum_quality = np.amin(qual)
    print(
        "There are "
        + str(len(vertices))
        + " vertices and "
        + str(len(entities))
        + " elements in the mesh",
        flush=True,
    )
    print("The minimum element quality is " + str(minimum_quality), flush=True)
    return vertices, entities
