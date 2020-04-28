# -----------------------------------------------------------------------------
#  Copyright (C) 2020 Keith Roberts

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
import numpy as np
import scipy.spatial as spspatial
import matplotlib.pyplot as plt
import time

from . import utils as mutils
from .. import decomp
from .. import migration

from .cpp import c_cgal


class MeshGenerator:  # noqa: C901
    """
    MeshGenerator: using sizng function and signed distance field build a mesh

    Usage
    -------
    >>> obj = MeshGenerator(MeshSizeFunction_obj,method='qhull')

    Parameters
    -------
        MeshSizeFunction_obj: self-explantatory
                        **kwargs
        method: verbose name of mesh generation method to use  (default=qhull) or cgal
    -------


    Returns
    -------
        A mesh object
    -------


    Example
    ------
    >>> msh = MeshSizeFunction(ef)
    -------

    """

    def __init__(self, SizingFunction, method="qhull"):
        self.SizingFunction = SizingFunction
        self.method = method

    # SETTERS AND GETTERS
    @property
    def SizingFunction(self):
        return self.__SizingFunction

    @SizingFunction.setter
    def SizingFunction(self, value):
        self.__SizingFunction = value

    @property
    def method(self):
        return self.__method

    @method.setter
    def method(self, value):
        self.__method = value

    ### PUBLIC METHODS ###
    def build(  # noqa: ignore=C901
        self,
        pfix=None,
        max_iter=10,
        nscreen=5,
        plot=False,
        seed=None,
        COMM=None,
        axis=0,
        points=None,
    ):
        """
        Interface to either DistMesh2D/3D mesh generator using signed distance functions.
        User has option to use either qhull or cgal for Del. retriangulation.

        Usage
        -----
        >>> p, t = build(self, max_iter=20)

        Parameters
        ----------
        pfix: points that you wish you constrain (default==None)
        max_iter: maximum number of iterations (default==20)
        nscreen: output to screen nscreen (default==5)
        plot: Visualize incremental meshes (default==False)
        seed: Random seed to ensure results are deterministic (default==None)
        COMM: MPI4py communicator (default==None)
        axis: axis to decomp the domain wrt (default==0)
        points: initial point distribution (default==None)

        Returns
        -------
        p:         Point positions (np, dim)
        t:         Triangle indices (nt, dim+1)
        """
        _ef = self.SizingFunction
        fd = _ef.fd
        fh = _ef.fh
        h0 = _ef.hmin
        bbox = _ef.bbox
        _method = self.method
        comm = COMM
        _axis = axis
        _points = points

        if comm is not None:
            PARALLEL = True
            rank = comm.Get_rank()
            size = comm.Get_size()
        else:
            PARALLEL = False
            rank = 0
            size = 1

        # set random seed to ensure deterministic results for mesh generator
        if seed is not None:
            if rank == 0:
                print("Setting psuedo-random number seed to " + str(seed), flush=True)
            np.random.seed(seed)

        dim = int(len(bbox) / 2)
        bbox = np.array(bbox).reshape(-1, 2)

        ptol = 0.001
        ttol = 0.1
        L0mult = 1 + 0.4 / 2 ** (dim - 1)
        deltat = 0.1
        geps = 1e-1 * h0
        deps = np.sqrt(np.finfo(np.double).eps) * h0

        if pfix is not None and not PARALLEL:
            pfix = np.array(pfix, dtype="d")
            nfix = len(pfix)
        else:
            pfix = np.empty((0, dim))
            nfix = 0

        if points is None:
            # 1. Create initial distribution in bounding box (equilateral triangles)
            p = np.mgrid[tuple(slice(min, max + h0, h0) for min, max in bbox)]
            p = p.reshape(dim, -1).T

            # 2. Remove points outside the region, apply the rejection method
            p = p[fd(p) < geps]  # Keep only d<0 points
            r0 = fh(p)
            p = np.vstack(
                (pfix, p[np.random.rand(p.shape[0]) < r0.min() ** dim / r0 ** dim])
            )
        else:
            # user has supplied initial points
            p = _points

        # we add jitter to avoid co-spherical points
        if PARALLEL:
            jitter = np.random.uniform(size=(len(p), dim), low=-h0 / 10, high=h0 / 10)
            p += jitter

        # call domain decomposition
        if PARALLEL:

            p, extents = decomp.blocker(points=p, rank=rank, nblocks=size, axis=_axis)

            N = p.shape[0]

            count = 0
            pold = float("inf")  # For first iteration

            print(
                "Commencing mesh generation with %d vertices on rank %d." % (N, rank),
                flush=True,
            )
        else:
            N = p.shape[0]

            count = 0
            pold = float("inf")  # For first iteration

            print("Commencing mesh generation with %d vertices." % N, flush=True)

        while True:

            # 3. Retriangulation by the Delaunay algorithm
            def dist(p1, p2):
                return np.sqrt(((p1 - p2) ** 2).sum(1))

            start = time.time()
            if (dist(p, pold) / h0).max() > ttol:  # Any large movement?

                # Make sure all points are unique
                p = np.unique(p, axis=0)
                pold = p.copy()  # Save current positions
                if _method == "qhull":
                    if PARALLEL:
                        tria = spspatial.Delaunay(p, incremental=True)
                        exports = migration.enqueue(
                            extents, p, tria.vertices, rank, size
                        )

                        recv = migration.exchange(comm, rank, size, exports)
                        tria.add_points(recv, restart=True)
                        p, t, inv = migration.utils.remove_external_faces(
                            tria.points, tria.vertices, extents[rank]
                        )
                        N = p.shape[0]
                        recv_ix = len(recv)  # we do not allow new points to move
                    else:
                        t = spspatial.Delaunay(p).vertices  # List of triangles
                elif _method == "cgal":
                    if dim == 2:
                        t = c_cgal.delaunay2(p[:, 0], p[:, 1])  # List of triangles
                    elif dim == 3:
                        t = c_cgal.delaunay3(
                            p[:, 0], p[:, 1], p[:, 2]
                        )  # List of triangles

                pmid = p[t].sum(1) / (dim + 1)  # Compute centroids
                t = t[fd(pmid) < -geps]  # Keep interior triangles

                # 4. Describe each bar by a unique pair of nodes
                if dim == 2:
                    bars = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
                elif dim == 3:
                    bars = np.concatenate(
                        [
                            t[:, [0, 1]],
                            t[:, [1, 2]],
                            t[:, [2, 0]],
                            t[:, [0, 3]],
                            t[:, [1, 3]],
                            t[:, [2, 3]],
                        ]
                    )
                bars = np.sort(bars, axis=1)
                bars = mutils.unique_rows(bars)  # Bars as node pairs

                # 5. Graphical output of the current mesh
                if plot and not PARALLEL:
                    if count % nscreen == 0:
                        if dim == 2:
                            plt.triplot(p[:, 0], p[:, 1], t)
                            plt.title("Retriangulation %d" % count)
                            plt.axis("equal")
                            plt.show()
                        elif dim == 3:
                            plt.title("Retriangulation %d" % count)

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

            Ftot = mutils.dense(
                bars[:, [0] * dim + [1] * dim],
                np.repeat([list(range(dim)) * 2], len(F), axis=0),
                np.hstack((Fvec, -Fvec)),
                shape=(N, dim),
            )

            Ftot[:nfix] = 0  # Force = 0 at fixed points

            if PARALLEL:
                if count < max_iter - 1:
                    p += deltat * Ftot
            else:
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
                dgrad2 = np.where(dgrad2 < deps, deps, dgrad2)
                p[ix] -= (d[ix] * np.vstack(dgrads) / dgrad2).T  # Project

            maxdp = deltat * np.sqrt((Ftot[d < -geps] ** 2).sum(1)).max()
            if count % nscreen == 0 and rank == 0:
                if PARALLEL:
                    print("On rank 0: ", flush=True)
                print(
                    "Iteration #%d, max movement is %f, there are %d vertices and %d cells"
                    % (count + 1, maxdp, len(p), len(t)),
                    flush=True,
                )

            # 8a. Termination criterion: All interior nodes move less than dptol (scaled)
            if maxdp < ptol * h0 and not PARALLEL:
                print(
                    "Termination reached...all interior nodes move less than dptol.",
                    flush=True,
                )

                if PARALLEL:
                    p, t = migration.aggregate(p, t, comm, size, rank)

                break

            # 8b. Number of iterations reached (only PARALLEL CAN EXIT THROUGH HERE)
            if count == max_iter - 1:
                if rank == 0:
                    print(
                        "Termination reached...maximum number of iterations reached.",
                        flush=True,
                    )

                if PARALLEL:
                    p, t = migration.aggregate(p, t, comm, size, rank)

                break

            # Delete new points
            if PARALLEL:
                p = np.delete(p, inv[-recv_ix::], axis=0)

                comm.barrier()

            count += 1

            end = time.time()
            if rank == 0:
                print("     Elapsed wall-clock time %f : " % (end - start), flush=True)
        return p, t
