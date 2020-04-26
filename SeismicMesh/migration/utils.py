import numpy as np

from ..geometry import signed_distance_functions as sdf


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


def fixmesh(p, t, ptol=2e-13):
    """Remove duplicated/unused nodes
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

    return p, t, jx
