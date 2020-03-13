# -----------------------------------------------------------------------------
#  Copyright (C) 2020 Keith Roberts

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
import numpy as np
import scipy.spatial as spspatial
import scipy.sparse as spsparse


def dense(I, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...), but instead returning a
    dense array.

    Usage
    -----
    >>> shape = (m, n)
    >>> A = dense(I, J, S, shape, dtype)
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
    return spsparse.coo_matrix((S, (II, J)), shape, dtype).toarray()


def setdiff_rows(A, B, return_index=False):
    """
    Similar to MATLAB's setdiff(A, B, 'rows'), this returns C, I
    where C are the row of A that are not in B and I satisfies
    C = A[I,:].

    Returns I if return_index is True.
    """
    A = np.require(A, requirements="C")
    B = np.require(B, requirements="C")

    assert A.ndim == 2, "array must be 2-dim'l"
    assert B.ndim == 2, "array must be 2-dim'l"
    assert A.shape[1] == B.shape[1], "arrays must have the same number of columns"
    assert A.dtype == B.dtype, "arrays must have the same data type"

    # NumPy provides setdiff1d, which operates only on one dimensional
    # arrays. To make the array one-dimensional, we interpret each row
    # as being a string of characters of the appropriate length.
    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize * ncolumns))
    C = (
        np.setdiff1d(A.view(dtype), B.view(dtype))
        .view(A.dtype)
        .reshape((-1, ncolumns), order="C")
    )
    if return_index:
        raise NotImplementedError
    else:
        return C


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


class MeshGenerator:
    """
    MeshGenerator: using sizng function and signed distance field build a mesh

    Usage
    -------
    >>>> obj = MeshGenerator(MeshSizeFunction_obj,method=DistMesh)


    Parameters
    -------
        MeshSizeFunction_obj: self-explantatory
                        **kwargs
        method: verbose name of mesh generation method to use  (default=DistMesh)


    Returns
    -------
        A mesh object


    Example
    ------
    msh = MeshSizeFunction(ef,fd)

    """

    def __init__(self, SizingFunction, method="DistMesh"):
        self.SizingFunction = SizingFunction

    # SETTERS AND GETTERS
    @property
    def SizingFunction(self):
        return self.__SizingFunction

    @SizingFunction.setter
    def SizingFunction(self, value):
        self.__SizingFunction = value

    ### PUBLIC METHODS ###
    def build(self, pfix=None, max_iter=10, nscreen=5, plot=False):
        """
        distmeshnd: 2D/3D mesh generator using distance functions.

        Usage
        -----
        >>> p, t = build(self, pfix=None, max_iter=20, plot=True)

        Parameters
        ----------
        pfix: points that you wish you constrain
        max_iter: maximum number of iterations
        nscreen: output to screen nscreen
        plot: Visualize incremental meshes

        Returns
        -------
        p:         Node positions (np, dim)
        t:         Triangle indices (nt, dim+1)
        """

        _ef = self.SizingFunction
        fd = _ef.fd
        fh = _ef.fh
        h0 = _ef.hmin
        bbox = _ef.bbox

        if plot:
            import matplotlib.pyplot as plt

        dim = int(len(bbox) / 2)
        bbox = np.array(bbox).reshape(2, -1)

        ptol = 0.001
        ttol = 0.1
        L0mult = 1 + 0.4 / 2 ** (dim - 1)
        deltat = 0.1
        geps = 1e-1 * h0
        deps = np.sqrt(np.finfo(np.double).eps) * h0

        if pfix is not None:
            pfix = np.array(pfix, dtype="d")
            nfix = len(pfix)
        else:
            pfix = np.empty((0, dim))
            nfix = 0

        # 1. Create initial distribution in bounding box (equilateral triangles)
        p = np.mgrid[tuple(slice(min, max + h0, h0) for min, max in bbox.T)]
        p = p.reshape(dim, -1).T

        # 2. Remove points outside the region, apply the rejection method
        p = p[fd(p) < geps]  # Keep only d<0 points
        r0 = fh(p)
        p = np.vstack(
            (pfix, p[np.random.rand(p.shape[0]) < r0.min() ** dim / r0 ** dim])
        )
        N = p.shape[0]

        count = 0
        pold = float("inf")  # For first iteration

        while True:

            # 3. Retriangulation by the Delaunay algorithm
            def dist(p1, p2):
                return np.sqrt(((p1 - p2) ** 2).sum(1))

            if (dist(p, pold) / h0).max() > ttol:  # Any large movement?
                pold = p.copy()  # Save current positions
                t = spspatial.Delaunay(p).vertices  # List of triangles
                pmid = p[t].sum(1) / (dim + 1)  # Compute centroids
                t = t[fd(pmid) < -geps]  # Keep interior triangles
                # 4. Describe each bar by a unique pair of nodes
                if dim == 2:
                    bars = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
                    bars = np.sort(bars, axis=1)
                elif dim == 3:
                    # TODO PROVIDE SUPPORT FOR 3D
                    bars = []
                bars = unique_rows(bars)  # Bars as node pairs
                # 5. Graphical output of the current mesh
                if plot:
                    if dim == 2:
                        plt.triplot(p[:, 0], p[:, 1], t)
                        plt.title("Retriangulation %d" % count)
                        plt.axis("equal")
                        plt.show()
                    elif dim == 3:
                        if count % 5 == 0:
                            # TODO ALL 3D VIZ
                            plt.title("Retriangulation %d" % count)
                            plt.axis("equal")
            # 6. Move mesh points based on bar lengths L and forces F
            barvec = p[bars[:, 0]] - p[bars[:, 1]]  # List of bar vectors
            L = np.sqrt((barvec ** 2).sum(1))  # L = Bar lengths
            hbars = fh(p[bars].sum(1) / 2)
            L0 = (
                hbars
                * L0mult
                * ((L ** dim).sum() / (hbars ** dim).sum()) ** (1.0 / dim)
            )
            F = L0 - L
            F[F < 0] = 0  # Bar forces (scalars)
            Fvec = (
                F[:, None] / L[:, None].dot(np.ones((1, dim))) * barvec
            )  # Bar forces (x,y components)
            Ftot = dense(
                bars[:, [0] * dim + [1] * dim],
                np.repeat([list(range(dim)) * 2], len(F), axis=0),
                np.hstack((Fvec, -Fvec)),
                shape=(N, dim),
            )
            Ftot[:nfix] = 0  # Force = 0 at fixed points
            p += deltat * Ftot  # Update node positions

            # 7. Bring outside points back to the boundary
            d = fd(p)
            ix = d > 0  # Find points outside (d>0)
            if ix.any():

                def deps_vec(i):
                    a = [0] * dim
                    a[i] = deps
                    return a

                dgrads = [(fd(p[ix] + deps_vec(i)) - d[ix]) / deps for i in range(dim)]
                dgrad2 = sum(dgrad ** 2 for dgrad in dgrads)
                p[ix] -= (d[ix] * np.vstack(dgrads) / dgrad2).T  # Project

            # 8a. Termination criterion: All interior nodes move less than dptol (scaled)
            maxdp = deltat * np.sqrt((Ftot[d < -geps] ** 2).sum(1)).max()
            if count % nscreen == 0:
                print(
                    "Iteration #%d, max movement is %f, there are %d vertices and %d cells"
                    % (count + 1, maxdp, len(p), len(t))
                )
            if maxdp < ptol * h0:
                print("Termination reached...all interior nodes move less than dptol.")
                break
            # 8b. Number of iterations reached.
            if count == max_iter - 1:
                print("Termination reached...maximum number of iterations reached.")
                break
            count += 1
        return p, t
