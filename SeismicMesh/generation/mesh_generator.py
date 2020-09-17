#   Copyright (C) 2020 Keith Roberts
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import math
import time
import warnings

import numpy as np
from mpi4py import MPI

from .. import decomp, geometry, migration
from .. import sizing
from . import utils as mutils
from .cpp.delaunay_class import DelaunayTriangulation as DT2
from .cpp.delaunay_class3 import DelaunayTriangulation3 as DT3

__all__ = ["sliver_removal", "generate_mesh"]

opts = {
    "nscreen": 1,
    "max_iter": 50,
    "seed": 0,
    "perform_checks": False,
    "pfix": None,
    "axis": 1,
    "min_dh_angle_bound": 10.0,
    "max_dh_angle_bound": 170.0,
    "points": None,
}


def sliver_removal(points, domain, cell_size, h0, comm=None, **kwargs):
    r"""Improve an existing 3D mesh by removing degenerate elements
    commonly referred to as `slivers`.

    :param points:
        An array of points that describe the vertices of an existing (higher-quality) mesh.
    :type filename: ``string``
    :param domain:
        A function that takes a point and returns the signed nearest distance to the domain boundary Ω
    :type domain: A :class:`geometry.Rectangle/Cube/Circle` object or a function object.
    :param cell_size:
        A function that can evalulate a point and return a mesh size.
    :type cell_size: A :class:`SizeFunction` object or a function object.
    :param h0:
        The minimum element size in the domain
    :type h0: `float`
    :param comm:
        MPI4py communicator
    :type comm: MPI4py communicator object, optional

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *nscreen* (``float``) --
            Output to the screen `nscreen` timestep. (default==1)
        * *max_iter* (``float``) --
            Maximum number of meshing iterations. (default==50)
        * *perform_checks* (`boolean`) --
            Whether or not to perform mesh linting/mesh cleanup. (default==False)
        * *pfix* (`array-like`) --
            An array of points to constrain in the mesh. (default==None)
        * *axis* (`int`) --
            The axis to decompose the mesh (1,2, or 3). (default==1)
        * *min_dh_angle_bound* (`float`) --
            The minimum allowable dihedral angle bound. (default==10 degrees)
        * *max_dh_angle_bound* (`float`) --
            The maximum allowable dihedral angle bound. (default==170 degrees)

    """
    comm = comm or MPI.COMM_WORLD
    if comm.rank > 0:
        if comm.rank == 1:
            warnings.warn("Sliver removal only works in serial for now")
        return True, True

    opts.update(kwargs)
    _parse_kwargs(kwargs)

    dim = points.shape[1]
    if dim == 2:
        if comm.rank == 0:
            raise Exception("Mesh improvement currently on works in 3D")

    # unpack domain
    fd, bbox0 = _unpack_domain(domain)

    fh, bbox1 = _unpack_sizing(cell_size)

    # take maxmin of boxes
    bbox = _minmax(bbox0, bbox1)

    # rebuild the Rectangle or cube if domain padding
    if bbox0 != bbox1:
        fd = geometry.Cube(bbox).eval

    if not isinstance(bbox, tuple):
        raise ValueError("`bbox` must be a tuple")

    if opts["max_iter"] < 0:
        raise ValueError("`max_iter` must be > 0")
    max_iter = opts["max_iter"]
    print("Will attempt " + str(max_iter) + " to bound the dihedral angles...")

    geps = 1e-1 * h0
    deps = np.sqrt(np.finfo(np.double).eps) * h0
    min_dh_bound = opts["min_dh_angle_bound"] * math.pi / 180
    max_dh_bound = opts["max_dh_angle_bound"] * math.pi / 180

    print(
        "Enforcing a min. dihedral bound of: "
        + str(min_dh_bound * 180 / math.pi)
        + " degrees..."
    )
    print(
        "Enforcing a max. dihedral bound of: "
        + str(max_dh_bound * 180 / math.pi)
        + " degrees..."
    )

    DT = _select_cgal_dim(dim)

    p = points

    N = len(p)

    print(
        "Commencing sliver removal with %d vertices on rank %d." % (N, comm.rank),
        flush=True,
    )

    count = 0
    pold = None
    nscreen = opts["nscreen"]

    while True:

        start = time.time()

        # Using CGAL's incremental Delaunay triangulation capabilities.
        if count == 0:
            dt = DT()
            dt.insert(p.flatten().tolist())
        else:
            to_move = np.where(_dist(p, pold) > 0)[0]
            dt.move(to_move.flatten().tolist(), p[to_move].flatten().tolist())

        # Get the current topology of the triangulation
        p, t = _get_topology(dt)

        pold = p.copy()

        # Remove points outside the domain
        t = _remove_triangles_outside(p, t, fd, geps)

        # Sliver removal
        if count != (max_iter - 1):
            num_move = 0
            # calculate dihedral angles in mesh
            dh_angles = geometry.calc_dihedral_angles(p, t)
            out_of_bounds = np.argwhere(
                (dh_angles[:, 0] < min_dh_bound) | (dh_angles[:, 0] > max_dh_bound)
            )
            ele_nums = np.floor(out_of_bounds / 6).astype("int")
            ele_nums, ix = np.unique(ele_nums, return_index=True)

            if count % nscreen == 0:
                print(
                    "On rank: "
                    + str(comm.rank)
                    + " There are "
                    + str(len(ele_nums))
                    + " slivers...",
                    flush=True,
                )

            move = t[ele_nums, 0]
            num_move = move.size
            if num_move == 0:
                print("Termination reached...no slivers detected!", flush=True)
                p, t, _ = geometry.fix_mesh(p, t, dim=dim, delete_unused=True)
                return p, t

            p0, p1, p2, p3 = (
                p[t[ele_nums, 0], :],
                p[t[ele_nums, 1], :],
                p[t[ele_nums, 2], :],
                p[t[ele_nums, 3], :],
            )

            # perturb vector is based on INCREASING circumsphere's radius
            perturb = geometry.calc_circumsphere_grad(p0, p1, p2, p3)
            perturb[np.isinf(perturb)] = 1.0

            # normalize perturbation vector
            perturb_norm = np.sum(np.abs(perturb) ** 2, axis=-1) ** (1.0 / 2)
            perturb /= perturb_norm[:, None]

            # perturb % of local mesh size
            p[move] += 0.10 * h0 * perturb

        # Bring outside points back to the boundary
        p = _project_points_back(p, fd, deps)

        # Number of iterations reached, stop.
        if count == (max_iter - 1):
            p, t = _termination(p, t, opts, comm)
            break

        count += 1

        end = time.time()
        if comm.rank == 0 and count % nscreen == 0:
            print("     Elapsed wall-clock time %f : " % (end - start), flush=True)

    return p, t


def generate_mesh(domain, cell_size, h0, comm=None, **kwargs):
    r"""Generate a 2D/3D triangulation using callbacks to a sizing function `cell_size` and signed distance function :class:`domain`

    :param domain:
        A function that takes a point and returns the signed nearest distance to the domain boundary Ω.
    :type domain: A :class:`geometry.Rectangle/Cube/Circle` object or a function object.
    :param cell_size:
        Either a :class:`size_function` object or a function that can evalulate a point and return a mesh size.
    :type cell_size: A :class:`cell_size` object or a function object.
    :param h0:
        The minimum element size in the domain.
    :type h0: `float`
    :param comm:
        MPI4py communicator
    :type comm: MPI4py communicator object, optional

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *bbox* (``tuple``) --
            Bounding box containing domain extents. REQUIRED IF NOT USING :class:`cell_size`
        * *nscreen* (``float``) --
            Output to the screen `nscreen` timestep. (default==1)
        * *max_iter* (``float``) --
            Maximum number of meshing iterations. (default==50)
        * *seed* (``float`` or ``int``) --
            Psuedo-random seed to initialize meshing points. (default==0)
        * *perform_checks* (`boolean`) --
            Whether or not to perform mesh linting/mesh cleanup. (default==False)
        * *pfix* (`array-like`) --
            An array of points to constrain in the mesh. (default==None)
        * *axis* (`int`) --
            The axis to decompose the mesh (1,2, or 3). (default==1)

    :return: points: vertex coordinates of mesh
    :rtype: points: (numpy.ndarray[`float` x dim])
    :return: t: mesh connectivity table.
    :rtype: t: numpy.ndarray[`int` x (dim + 1)]

    """
    comm = comm or MPI.COMM_WORLD

    # check call was correct
    opts.update(kwargs)
    _parse_kwargs(kwargs)

    # unpack domain
    fd, bbox0 = _unpack_domain(domain)

    fh, bbox1 = _unpack_sizing(cell_size)

    # take maxmin of boxes
    bbox = _minmax(bbox0, bbox1)

    if not isinstance(bbox, tuple):
        raise ValueError("`bbox` must be a tuple")

    # check bbox shape
    dim = int(len(bbox) / 2)

    # rebuild the Rectangle or Cube if domain padding
    if bbox0 != bbox1:
        if dim == 2:
            tmp = geometry.Rectangle(bbox)
        elif dim == 3:
            tmp = geometry.Cube(bbox)

        fd = tmp.eval

    bbox = np.array(bbox).reshape(-1, 2)

    # check h0
    if h0 < 0:
        raise ValueError("`h0` must be > 0")

    if opts["max_iter"] < 0:
        raise ValueError("`max_iter` must be > 0")
    max_iter = opts["max_iter"]

    np.random.seed(opts["seed"])

    """These parameters originate from the original DistMesh"""
    L0mult = 1 + 0.4 / 2 ** (dim - 1)
    delta_t = 0.1
    geps = 1e-1 * h0
    deps = np.sqrt(np.finfo(np.double).eps) * h0

    DT = _select_cgal_dim(dim)

    pfix, nfix = _unpack_pfix(dim, opts, comm)

    fh, p, extents = _initialize_points(dim, geps, bbox, fh, fd, h0, opts, pfix, comm)

    N = p.shape[0]

    assert N > 0, "No vertices to mesh with!"

    count = 0
    print(
        "Commencing mesh generation with %d vertices on rank %d." % (N, comm.rank),
        flush=True,
    )

    nscreen = opts["nscreen"]
    while True:

        start = time.time()

        # (Re)-triangulation by the Delaunay algorithm
        dt = DT()
        dt.insert(p.flatten().tolist())

        # Get the current topology of the triangulation
        p, t = _get_topology(dt)

        # Add ghost points to perform Delaunay in parallel.
        if comm.size > 1:
            p, t, inv, recv_ix = _add_ghost_vertices(p, t, dt, extents, comm)

        # Remove points outside the domain
        t = _remove_triangles_outside(p, t, fd, geps)

        # Compute the forces on the bars
        Ftot = _compute_forces(p, t, fh, h0, L0mult)

        Ftot[:nfix] = 0  # Force = 0 at fixed points

        # Last timestep in parallel, we don't move points
        if comm.size > 1:
            if count < max_iter - 1:
                p += delta_t * Ftot
        else:
            p += delta_t * Ftot  # Update positions

        # Bring outside points back to the boundary
        p = _project_points_back(p, fd, deps)

        # Number of iterations reached, stop.
        if count == (max_iter - 1):
            p, t = _termination(p, t, opts, comm)
            break

        if comm.size > 1:
            # If continuing on, delete ghost points
            p = np.delete(p, inv[-recv_ix::], axis=0)
            comm.barrier()

        # Show the user some progress so they know something is happening
        if count % nscreen == 0 and comm.rank == 0:
            maxdp = delta_t * np.sqrt((Ftot ** 2).sum(1)).max()
            _display_progress(p, t, count, nscreen, maxdp, comm)

        count += 1

        end = time.time()
        if comm.rank == 0 and count % nscreen == 0:
            print("     Elapsed wall-clock time %f : " % (end - start), flush=True)

    return p, t


def _minmax(bbox0, bbox1):
    d = []
    for i, (a, b) in enumerate(zip(bbox0, bbox1)):
        a, b = (a, b)
        if i % 2:
            c = max(a, b)
        else:
            c = min(a, b)
        d.append(c)
    return tuple(d)


def _unpack_sizing(cell_size):
    if isinstance(cell_size, sizing.SizeFunction):
        bbox = cell_size.bbox
        fh = cell_size.eval
    elif callable(cell_size):
        bbox = opts["bbox"]
        fh = cell_size
    else:
        raise ValueError(
            "`cell_size` must either be a function or a `cell_size` object"
        )
    return fh, bbox


def _unpack_domain(domain):
    if isinstance(domain, geometry.Rectangle):
        bbox = (domain.x1, domain.x2, domain.y1, domain.y2)
        fd = domain.eval
    elif isinstance(domain, geometry.Cube):
        bbox = (domain.x1, domain.x2, domain.y1, domain.y2, domain.z1, domain.z2)
        fd = domain.eval
    elif callable(domain):
        # get the bbox from the name value pairs or quit
        bbox = opts["bbox"]
        fd = domain
    else:
        raise ValueError("`domain` must be a function or a :class:`geometry` object")
    return fd, bbox


def _parse_kwargs(kwargs):
    for key in kwargs:
        if key in {
            "nscreen",
            "max_iter",
            "seed",
            "perform_checks",
            "pfix",
            "axis",
            "points",
            "domain",
            "cell_size",
            "bbox",
            "min_dh_angle_bound",
            "max_dh_angle_bound",
        }:
            pass
        else:
            raise ValueError(
                "Option %s with parameter %s not recognized " % (key, kwargs[key])
            )


def _display_progress(p, t, count, nscreen, maxdp, comm):
    """print progress"""
    print(
        "Iteration #%d, max movement is %f, there are %d vertices and %d cells"
        % (count + 1, maxdp, len(p), len(t)),
        flush=True,
    )


def _termination(p, t, opts, comm):
    """Shut it down when reacing `max_iter`"""
    dim = p.shape[1]
    if comm.rank == 0:
        print("Termination reached...maximum number of iterations reached.", flush=True)
    if comm.size > 1:
        # gather onto rank 0
        p, t = migration.aggregate(p, t, comm, comm.size, comm.rank, dim=dim)
    # perform linting if asked
    if comm.rank == 0 and opts["perform_checks"]:
        p, t = geometry.linter(p, t, dim=dim)
    elif comm.rank == 0:
        p, t, _ = geometry.fix_mesh(p, t, dim=dim, delete_unused=True)
    return p, t


def _get_bars(t):
    """Describe each bar by a unique pair of nodes"""
    dim = t.shape[1] - 1
    bars = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    if dim == 3:
        bars = np.concatenate((bars, t[:, [0, 3]], t[:, [1, 3]], t[:, [2, 3]]), axis=0)
    bars = mutils.unique_rows(
        np.ascontiguousarray(bars, dtype=np.uint32)
    )  # Bars as node pairs
    return bars


def _compute_forces(p, t, fh, h0, L0mult):
    """Compute the forces on each edge based on the sizing function"""
    dim = p.shape[1]
    N = p.shape[0]
    bars = _get_bars(t)
    barvec = p[bars[:, 0]] - p[bars[:, 1]]  # List of bar vectors
    L = np.sqrt((barvec ** 2).sum(1))  # L = Bar lengths
    L[L == 0] = np.finfo(float).eps
    hbars = fh(p[bars].sum(1) / 2)
    L0 = hbars * L0mult * ((L ** dim).sum() / (hbars ** dim).sum()) ** (1.0 / dim)
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
    return Ftot


def _add_ghost_vertices(p, t, dt, extents, comm):
    """Parallel Delauany triangulation requires ghost vertices
    to be added each meshing iteration to maintain Delaunay-hood
    """
    dim = p.shape[1]
    exports = migration.enqueue(extents, p, t, comm.rank, comm.size, dim=dim)
    recv = migration.exchange(comm, comm.rank, comm.size, exports, dim=dim)
    recv_ix = len(recv)
    dt.insert(recv.flatten().tolist())
    p, t = _get_topology(dt)
    p, t, inv = geometry.remove_external_entities(
        p,
        t,
        extents[comm.rank],
        dim=dim,
    )
    return p, t, inv, recv_ix


def _remove_triangles_outside(p, t, fd, geps):
    """Remove vertices outside the domain"""
    dim = p.shape[1]
    pmid = p[t].sum(1) / (dim + 1)  # Compute centroids
    return t[fd(pmid) < -geps]  # Keep interior triangles


def _project_points_back(p, fd, deps):
    """Project points outsidt the domain back within"""
    dim = p.shape[1]

    d = fd(p)
    ix = d > 0  # Find points outside (d>0)
    if ix.any():

        def _deps_vec(i):
            a = [0] * dim
            a[i] = deps
            return a

        dgrads = [(fd(p[ix] + _deps_vec(i)) - d[ix]) / deps for i in range(dim)]
        dgrad2 = sum(dgrad ** 2 for dgrad in dgrads)
        dgrad2 = np.where(dgrad2 < deps, deps, dgrad2)
        p[ix] -= (d[ix] * np.vstack(dgrads) / dgrad2).T  # Project
    return p


def _user_defined_points(dim, fh, h0, bbox, points, comm, opts):
    """If the user has supplied initial points"""
    if comm.size > 1:
        # Domain decompose and localize points
        p = None
        if comm.rank == 0:
            blocks, extents = decomp.blocker(
                points=points, rank=comm.rank, num_blocks=comm.size, axis=opts["axis"]
            )
        else:
            blocks = None
            extents = None
        fh = migration.localize_sizing_function(fh, h0, bbox, dim, opts["axis"], comm)
        # send points to each subdomain
        p, extents = migration.localize_points(blocks, extents, comm, dim)
    else:
        p = points
        extents = None
    return fh, p, extents


def _generate_initial_points(h0, geps, dim, bbox, fh, fd, pfix, comm, opts):
    """User did not specify initial points"""
    if comm.size > 1:
        # Localize mesh size function grid.
        fh = migration.localize_sizing_function(fh, h0, bbox, dim, opts["axis"], comm)
        # Create initial points in parallel in local box owned by rank
        p = mutils.make_init_points(bbox, comm.rank, comm.size, opts["axis"], h0, dim)
    else:
        # Create initial distribution in bounding box (equilateral triangles)
        p = np.mgrid[tuple(slice(min, max + h0, h0) for min, max in bbox)].astype(float)
        p = p.reshape(dim, -1).T

    # Remove points outside the region, apply the rejection method
    p = p[fd(p) < geps]  # Keep only d<0 points
    r0 = fh(p)
    r0m = r0.min()
    # Make sure decimation occurs uniformly accross ranks
    if comm.size > 1:
        r0m = comm.allreduce(r0m, op=MPI.MIN)
    p = np.vstack(
        (
            pfix,
            p[np.random.rand(p.shape[0]) < r0m ** dim / r0 ** dim],
        )
    )
    extents = _form_extents(p, h0, comm)
    return fh, p, extents


def _initialize_points(dim, geps, bbox, fh, fd, h0, opts, pfix, comm):
    """Form initial point set to mesh with"""
    points = opts["points"]
    if points is None:
        # def _generate_initial_points(h0, geps, dim, bbox, fh, fd, pfix, comm, opts):
        fh, p, extents = _generate_initial_points(
            h0, geps, dim, bbox, fh, fd, pfix, comm, opts
        )
    else:
        fh, p, extents = _user_defined_points(dim, fh, h0, bbox, points, comm, opts)
    return fh, p, extents


def _form_extents(p, h0, comm):
    dim = p.shape[1]
    _axis = opts["axis"]
    if comm.size > 1:
        # min x min y min z max x max y max z
        extent = [*np.amin(p, axis=0), *np.amax(p, axis=0)]
        extent[_axis] -= h0
        extent[_axis + dim] += h0
        return [comm.bcast(extent, r) for r in range(comm.size)]
    else:
        return []


def _dist(p1, p2):
    """Euclidean distance between two sets of points"""
    return np.sqrt(((p1 - p2) ** 2).sum(1))


def _unpack_pfix(dim, opts, comm):
    """Unpack fixed points"""
    if opts["pfix"] is not None:
        if comm.size > 1:
            raise Exception("Fixed points are not yet supported in parallel.")
        else:
            pfix = np.array(opts["pfix"], dtype="d")
            nfix = len(pfix)
        if comm.rank == 0:
            print(
                "Constraining " + str(nfix) + " fixed points..",
                flush=True,
            )
    else:
        pfix = np.empty((0, dim))
        nfix = 0
    return pfix, nfix


def _select_cgal_dim(dim):
    """Select back-end CGAL Delaunay call"""
    if dim == 2:
        return DT2
    elif dim == 3:
        return DT3


def _get_topology(dt):
    """ Get points and entities from :clas:`CGAL:DelaunayTriangulation2/3` object"""
    p = dt.get_finite_vertices()
    t = dt.get_finite_cells()
    return p, t
