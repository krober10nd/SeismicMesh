# -----------------------------------------------------------------------------
#  Copyright (C) 2020 Keith Roberts

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
import numpy as np
from mpi4py import MPI
import scipy.spatial as spspatial
import matplotlib.pyplot as plt
import time

from . import utils as mutils
from .. import decomp
from .. import migration
from .. import geometry

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
        perform_checks=False,
        mesh_improvement=False,
        max_dh_bound=170,
        min_dh_bound=10,
        improvement_method="circumsphere",  # also volume or random
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
        nscreen: output to screen every nscreen iterations (default==5)
        plot: Visualize incremental meshes (default==False)
        seed: Random seed to ensure results are deterministic (default==None)
        COMM: MPI4py communicator (default==None)
        axis: axis to decomp the domain wrt (default==0)
        points: initial point distribution (default==None)
        perform_checks: run serial linting (slow)
        mesh_improvement: run mesh improvement (default=False)
        min_dh_bound: minimum dihedral angle allowed (default=5)
        improvement_method: method to perturb slivers (default=circumsphere)

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

        L0mult = 1 + 0.4 / 2 ** (dim - 1)
        deltat = 0.1
        geps = 1e-1 * h0
        deps = np.sqrt(np.finfo(np.double).eps) * h0

        if pfix is not None:
            pfix = np.array(pfix, dtype="d")
            nfix = len(pfix)
            if rank == 0:
                print(
                    "INFO: Constraining " + str(nfix) + " fixed points..", flush=True,
                )
        else:
            pfix = np.empty((0, dim))
            nfix = 0

        if _points is None:
            p = None
            if rank == 0:
                # 1. Create initial distribution in bounding box (equilateral triangles)
                p = np.mgrid[
                    tuple(slice(min, max + h0, h0) for min, max in bbox)
                ].astype(float)
                p = p.reshape(dim, -1).T

                # 2. Remove points outside the region, apply the rejection method
                p = p[fd(p) < geps]  # Keep only d<0 points
                r0 = fh(p)
                p = np.vstack(
                    (pfix, p[np.random.rand(p.shape[0]) < r0.min() ** dim / r0 ** dim],)
                )
            USER_DEFINED_POINTS = False
        else:
            # If the user has supplied initial points
            if PARALLEL:
                p = None
                if rank == 0:
                    p = _points
                USER_DEFINED_POINTS = True
            else:
                USER_DEFINED_POINTS = True
                p = _points

        # 2b. Call domain decomposition and localize points
        if PARALLEL:
            if rank == 0:
                blocks, extents = decomp.blocker(
                    points=p, rank=rank, nblocks=size, axis=_axis
                )
            else:
                blocks = None
                extents = None
            # send points to each subdomain
            p, extents = migration.localize(blocks, extents, comm, dim)

        N = len(p)

        count = 0
        print(
            "Commencing mesh generation with %d vertices on rank %d." % (N, rank),
            flush=True,
        )
        while True:

            # 3. Retriangulation by the Delaunay algorithm
            start = time.time()

            # Using the SciPy qhull wrapper for Delaunay triangulation.
            if _method == "qhull":
                if PARALLEL:
                    tria = spspatial.Delaunay(p, incremental=True)
                    # This greatly avoids coplanar and colinear points (just done once)
                    if count == 0 and not USER_DEFINED_POINTS:
                        jitter = np.random.uniform(
                            size=(len(p), dim), low=-h0 / 10, high=h0 / 10
                        )
                        p += jitter
                    exports = migration.enqueue(
                        extents, p, tria.simplices, rank, size, dim=dim
                    )
                    recv = migration.exchange(comm, rank, size, exports, dim=dim)
                    tria.add_points(recv, restart=True)
                    p, t, inv = geometry.remove_external_faces(
                        tria.points, tria.simplices, extents[rank], dim=dim
                    )
                    N = p.shape[0]
                    recv_ix = len(recv)  # we do not allow new points to move
                else:
                    # SERIAL
                    t = spspatial.Delaunay(p).vertices  # List of triangles
            # Using CGAL's Delaunay triangulation algorithm
            elif _method == "cgal":
                if PARALLEL:
                    print(
                        "Parallel mesh generation with CGAL not yet supported!",
                        flush=True,
                    )
                    quit()
                else:
                    # SERIAL
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
            bars = mutils.unique_rows(bars)  # Bars as node pairs
            bars = bars[0]

            # 5a. Graphical output of the current mesh
            if plot and not PARALLEL:
                if count % nscreen == 0:
                    if dim == 2:
                        plt.triplot(p[:, 0], p[:, 1], t)
                        plt.title("Retriangulation %d" % count)
                        plt.axis("equal")
                        plt.show()
                    elif dim == 3:
                        plt.title("Retriangulation %d" % count)

            # Slow the movement of points periodically as things converge
            if mesh_improvement:
                deltat /= 2.0

            # Periodic sliver removal
            num_move = 0
            if mesh_improvement and count != (max_iter - 1) and dim == 3:
                dh_angles = np.rad2deg(geometry.calc_dihedral_angles(p, t))
                outOfBounds = np.argwhere(
                    (dh_angles[:, 0] < min_dh_bound) | (dh_angles[:, 0] > max_dh_bound)
                )
                eleNums = np.floor(outOfBounds / 6).astype("int")
                eleNums, ix = np.unique(eleNums, return_index=True)
                if count % nscreen == 0:
                    print(
                        "On rank: "
                        + str(rank)
                        + " There are "
                        + str(len(eleNums))
                        + " slivers...",
                        flush=True,
                    )
                move = t[eleNums, 0]
                num_move = move.size
                if PARALLEL:
                    g_num_move = comm.allreduce(num_move, op=MPI.SUM)
                    if num_move == 0 and g_num_move != 1:
                        if count % nscreen == 0:
                            print("Rank " + str(rank) + " is locked...", flush=True)
                        nfix = N
                    if g_num_move == 0:
                        if rank == 0:
                            print(
                                "Terimation reached...No slivers detected!", flush=True
                            )
                        p, t = migration.aggregate(p, t, comm, size, rank, dim=dim)
                        return p, t
                else:
                    if num_move == 0:
                        print("Terimation reached...No slivers detected!", flush=True)
                        return p, t

                p0, p1, p2, p3 = (
                    p[t[eleNums, 0], :],
                    p[t[eleNums, 1], :],
                    p[t[eleNums, 2], :],
                    p[t[eleNums, 3], :],
                )

                # perturb the points associated with the out-of-bound dihedral angle
                if improvement_method == "random":
                    # pertubation vector is random
                    perturb = np.random.uniform(size=(num_move, dim), low=-1, high=1)

                if improvement_method == "volume":
                    # perturb vector is based on REDUCING slivers volume
                    perturb = geometry.calc_volume_grad(p1, p2, p3)

                if improvement_method == "circumsphere":
                    # perturb vector is based on INCREASING circumsphere's radius
                    perturb = geometry.calc_circumsphere_grad(p0, p1, p2, p3)
                    perturb[np.isinf(perturb)] = 1.0

                # normalize
                perturb_norm = np.sum(np.abs(perturb) ** 2, axis=-1) ** (1.0 / 2)
                perturb /= perturb_norm[:, None]

                # perturb % of minimum mesh size
                mesh_size = fh(p[move])
                if PARALLEL:
                    if count < max_iter - 1:
                        if improvement_method == "volume":
                            p[move] -= 0.10 * mesh_size[:, None] * perturb
                        else:
                            p[move] += 0.10 * mesh_size[:, None] * perturb

                else:
                    if improvement_method == "volume":
                        p[move] -= 0.10 * mesh_size[:, None] * perturb
                    else:
                        p[move] += 0.10 * mesh_size[:, None] * perturb

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

            # lock the moved point
            if num_move > 0:
                Ftot[move] = 0.0

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

            # 8. Number of iterations reached
            if count == max_iter - 1:
                if rank == 0:
                    print(
                        "Termination reached...maximum number of iterations reached.",
                        flush=True,
                    )
                if PARALLEL:
                    p, t = migration.aggregate(p, t, comm, size, rank, dim=dim)
                    if rank == 0 and perform_checks:
                        p, t = geometry.linter(p, t, dim=dim)
                break

            # Delete new points
            if PARALLEL:
                p = np.delete(p, inv[-recv_ix::], axis=0)

                comm.barrier()

            if count % nscreen == 0 and rank == 0:
                if PARALLEL:
                    print("On rank 0: ", flush=True)
                print(
                    "Iteration #%d, max movement is %f, there are %d vertices and %d cells"
                    % (count + 1, maxdp, len(p), len(t)),
                    flush=True,
                )

            count += 1

            end = time.time()
            if rank == 0 and count % nscreen == 0:
                print("     Elapsed wall-clock time %f : " % (end - start), flush=True)
        return p, t

    def parallel_build(  # noqa: ignore=C901
        self,
        pfix=None,
        max_iter=10,
        nscreen=5,
        plot=False,
        seed=None,
        COMM=None,
        axis=0,
        points=None,
        perform_checks=False,
    ):
        """
        Thin wrapper for build to simplify user interaction when building in parallel
        See MeshGenerator.build for inputs
        """
        p, t = self.build(
            pfix=pfix,
            max_iter=max_iter,
            nscreen=nscreen,
            plot=plot,
            seed=seed,
            COMM=COMM,
            axis=axis,
            points=None,
            perform_checks=False,
        )
        # alternate the decomposed axis
        if axis == 1:
            axis == 0
        elif axis == 0:
            axis = 1
        # finalize mesh (switch decomposition axis and perform serial linting)
        # improves mesh quality near decomp boundaries
        p, t = self.build(
            pfix=pfix,
            max_iter=30,
            nscreen=nscreen,
            plot=plot,
            seed=seed,
            COMM=COMM,
            axis=axis,
            points=p,
            perform_checks=True,
        )
        return p, t
