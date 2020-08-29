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

import warnings

import h5py
import matplotlib.pyplot as plt
import numpy as np
import segyio
from mpi4py import MPI
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

from ..geometry import signed_distance_functions as sdf
from .cpp import limgrad


def ReadSegy(fname):
    """Read a velocity model from a SEG-y file"""
    with segyio.open(fname, ignore_geometry=True) as f:
        nz, nx = len(f.samples), len(f.trace)
        vp = np.zeros(shape=(nz, nx))
        for index, trace in enumerate(f.trace):
            vp[:, index] = trace
        return np.flipud(vp)


class MeshSizeFunction:
    """The :class:`MeshSizeFunction` is used to build a rectangular or cubic isotropic mesh size function :math:`f(h)`.
    and assumes the domain is a rectangle (2D) or cube (3D).
    """

    def __init__(
        self,
        bbox,
        hmin,
        velocity_grid,
        units="m-s",
        wl=0.0,
        freq=5.0,
        grad=0.0,
        space_order=1,
        hmax=np.inf,
        dt=0.0,
        cr_max=1.0,
        grade=0.0,
        nx=None,
        ny=None,
        nz=None,
        domain_ext=0.0,
        padstyle="edge",
    ):
        """Class constructor for :class:`MeshSizeFunction`

        :param bbox: bounding box containing domain extents.
        :type bbox: tuple with size (2*dim). For example, in 2D `(zmin, zmax, xmin, xmax)`
        :param hmin: minimum triangular edgelength populating the domain in meters.
        :type hmin: float64
        :param velocity_grid: a grid of velocity values.
        :type velocity_grid: array-like, numpy.ndarray()
        :param units: units of the velocity model (either `m-s` or `km-s`)
        :type units: str, optional, default=`m-s`
        :param wl: number of vertices per wavelength for a given :math:`f_{max}`
        :type wl: int, optional
        :param grad: the resolution in m nearby sharp gradients in velociy.
        :type grad: float64, optional
        :param freq: :math:`f_{max}` in hertz for which to estimate `wl`
        :type freq: float64, optional
        :param space_order: the polynomial order of the basis functions.
        :type space_order: int, optional
        :param hmax: maximum mesh size in meters allowed in the domain
        :type hmax: float64, optional
        :param dt: theoretical maximum stable timestep in seconds given Courant number `Cr`
        :type dt: float64, optional
        :param cr_max: `dt` is stable with this Courant number.
        :type cr_max: float64, optional
        :param grade: maximum allowable variation in mesh size in decimal percent.
        :type grade: float64, optional
        :param domain_ext: width of domain extension in `-z`, `+x`, `-x`, `+y`, `-y` directions
        :type domain_ext: float64, optional
        :param padstyle: method to pad velocity in the domain extension
        :type padstyle: str, optional, `edge`, `linear_ramp`, `constant`

        :return: object populated with meta-data, :math:`f(h)`, and a :math:`f(d)`.
        :rtype: :class:`MeshSizeFunction` object
        """

        self.bbox = bbox
        self.velocity_grid = velocity_grid
        self.dim = int(len(self.bbox) / 2)
        self.width = bbox[3] - bbox[2]
        self.depth = bbox[1] - bbox[0]
        if self.dim == 3:
            self.length = bbox[5] - bbox[4]
        self.spacingZ = None
        self.spacingX = None
        self.units = units
        self.hmin = hmin
        self.hmax = hmax
        self.wl = wl
        self.freq = freq
        self.grad = grad
        self.space_order = space_order
        self.dt = dt
        self.cr_max = cr_max
        self.grade = grade
        self.fh = None
        self.fd = None
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.domain_ext = domain_ext
        self.padstyle = padstyle
        self.interpolant = None

    ### SETTERS AND GETTERS ###

    @property
    def interpolant(self):
        return self.__interpolant

    @interpolant.setter
    def interpolant(self, value):
        self.__interpolant = value

    @property
    def fh(self):
        return self.__fh

    @fh.setter
    def fh(self, value):
        self.__fh = value

    @property
    def fd(self):
        return self.__fd

    @fd.setter
    def fd(self, value):
        self.__fd = value

    @property
    def bbox(self):
        return self.__bbox

    @bbox.setter
    def bbox(self, value):
        assert (
            len(value) >= 4 and len(value) <= 6
        ), "bbox has wrong number of values. either 4 or 6."
        self.__bbox = value

    @property
    def hmin(self):
        return self.__hmin

    @hmin.setter
    def hmin(self, value):
        assert value > 0.0, "hmin must be non-zero"
        self.__hmin = value

    @property
    def dim(self):
        return self.__dim

    @dim.setter
    def dim(self, value):
        assert value == 2 or value == 3, "dim must be either 2 or 3"
        self.__dim = value

    @property
    def velocity_grid(self):
        return self.__velocity_grid

    @velocity_grid.setter
    def velocity_grid(self, value):
        if value is None:
            raise ValueError("Velocity grid must be a numpy array")
        if np.amin(value) < 1000:
            warnings.warn("Min. velocity < 1000 m-s. Units may be incorrect.")
        if np.amax(value) > 10000:
            warnings.warn("Max. velocity > 10,000 m-s. Units may be incorrect.")
        self.__velocity_grid = value

    @property
    def units(self):
        return self.__units

    @units.setter
    def units(self, value):
        assert value == "m-s" or value == "km-s", "units are not compatible"
        self.__units = value

    @property
    def wl(self):
        return self.__wl

    @wl.setter
    def wl(self, value):
        self.__wl = value

    @property
    def grad(self):
        return self.__grad

    @grad.setter
    def grad(self, value):
        self.__grad = value

    @property
    def freq(self):
        return self.__freq

    @freq.setter
    def freq(self, value):
        self.__freq = value

    @property
    def space_order(self):
        return self.__space_order

    @space_order.setter
    def space_order(self, value):
        self.__space_order = value

    @property
    def hmax(self):
        return self.__hmax

    @hmax.setter
    def hmax(self, value):
        self.__hmax = value

    @property
    def dt(self):
        return self.__dt

    @dt.setter
    def dt(self, value):
        assert value >= 0, "dt must be > 0"
        self.__dt = value

    @property
    def cr_max(self):
        return self.__cr_max

    @cr_max.setter
    def cr_max(self, value):
        assert value >= 0, "Cr_max must be > 0"
        self.__cr_max = value

    @property
    def grade(self):
        return self.__grade

    @grade.setter
    def grade(self, value):
        assert value >= 0, "grade must be > 0"
        self.__grade = value

    @property
    def domain_ext(self):
        return self.__domain_ext

    @domain_ext.setter
    def domain_ext(self, value):
        assert value >= 0, "domain extent must be > 0"
        self.__domain_ext = value

    @property
    def padstyle(self):
        return self.__padstyle

    @padstyle.setter
    def padstyle(self, value):
        assert value == "edge" or value == "constant" or value == "linear_ramp"
        self.__padstyle = value

    # ---PUBLIC METHODS---#

    def build(self, comm=None):  # noqa: C901

        """Builds the isotropic mesh size function according
        to the user arguments that were passed.
        """
        comm = comm or MPI.COMM_WORLD
        if comm is not None:
            rank = comm.Get_rank()
            size = comm.Get_size()
        else:
            rank = 0
            size = 1

        if rank == 0:
            self.__ReadVelocityModel()
            _vp = self.velocity_grid

            _bbox = self.bbox
            _dim = self.dim
            _width = self.width
            _nz = self.nz
            _nx = self.nx
            if _dim == 3:
                _ny = self.ny
            _domain_ext = self.domain_ext

            _hmax = self.hmax
            _hmin = self.hmin
            _grade = self.grade

            _wl = self.wl
            _freq = self.freq
            _space_order = self.space_order

            _grad = self.grad

            _dt = self.dt
            _cr_max = self.cr_max

            if _dim == 2:
                hh_m = np.zeros(shape=(_nz, _nx), dtype=np.float32) + _hmin
            if _dim == 3:
                hh_m = np.zeros(shape=(_nz, _nx, _ny), dtype=np.float32) + _hmin

            if _wl > 0:
                if rank == 0:
                    print(
                        "Mesh sizes will be built to resolve an estimate of wavelength with "
                        + str(_wl)
                        + " vertices...",
                        flush=True,
                    )
                hh_m = _vp / (_freq * _wl)

            if _grad > 0:
                if rank == 0:
                    print(
                        "Refining mesh sizes near sharp velocity gradients...",
                        flush=True,
                    )
                if _dim == 2:
                    win_rows, win_cols = 100, 100
                    win_mean = ndimage.uniform_filter(_vp, (win_rows, win_cols))
                    win_sqr_mean = ndimage.uniform_filter(
                        _vp ** 2, (win_rows, win_cols)
                    )
                elif _dim == 3:
                    win_rows, win_cols, win_cols2 = 100, 100, 100
                    win_mean = ndimage.uniform_filter(
                        _vp, (win_rows, win_cols, win_cols2)
                    )
                    win_sqr_mean = ndimage.uniform_filter(
                        _vp ** 2, (win_rows, win_cols, win_cols2)
                    )
                win_var = win_sqr_mean - win_mean ** 2
                # normalize variance to [0,1]
                win_var /= np.amax(win_var)
                win_var -= np.amin(win_var)
                tmp = _grad / (win_var + 0.10)
                hh_m = np.minimum(tmp, hh_m)

            # enforce min (and optionally max) sizes
            hh_m = np.where(hh_m < _hmin, _hmin, hh_m)
            if _hmax < np.inf:
                if rank == 0:
                    print("Enforcing maximum mesh resolution...", flush=True)
                hh_m = np.where(hh_m > _hmax, _hmax, hh_m)

            # adjust mesh res. based on the CFL limit so cr < cr_max
            # NB: considering the space order p
            if _dt > 0:
                if rank == 0:
                    print(
                        "Enforcing timestep of " + str(_dt) + " seconds...", flush=True
                    )
                # assume isotropic mesh resolution
                cr_old = (_vp * _dt) / (_dim * hh_m)
                # divide maximum Courant by desired space order
                _cr_max = _cr_max / (_dim * _space_order)
                # determine resolution that statisfies desired _cr_max
                dxn = (_vp * _dt) / (_dim * _cr_max)
                # edit mesh size function to avoid violations of cr
                hh_m = np.where(cr_old > _cr_max, dxn, hh_m)

            # grade the mesh sizes
            if _grade > 0:
                if rank == 0:
                    # TODO: grader needs to consider varying constant resolution by dimension
                    print("Enforcing mesh gradation in sizing function...", flush=True)
                    hh_m = self.hj(hh_m, _width / _nx, 10000)

            # Domain extension
            # edit the bbox to reflect the new domain size
            if _domain_ext > 0:
                self = self.__CreateDomainExtension()  # edit the bbox
                hh_m = self.__EditMeshSizeFunction(
                    hh_m, rank
                )  # edit the sizing function

            # construct a interpolator object to be queried during mesh generation
            if rank == 0:
                print("Building a gridded interpolant...", flush=True)

            if _dim == 2:
                z_vec, x_vec = self.__CreateDomainVectors()

            if _dim == 3:
                z_vec, x_vec, y_vec = self.__CreateDomainVectors()

            assert np.all(hh_m > 0.0), "Mesh size_function must be strictly positive."

            if _dim == 2:
                self.interpolant = RegularGridInterpolator(
                    (z_vec, x_vec), hh_m, bounds_error=False, fill_value=None
                )
            if _dim == 3:
                self.interpolant = RegularGridInterpolator(
                    (z_vec, x_vec, y_vec), hh_m, bounds_error=False, fill_value=None
                )
            # Python can't bcast pickles so this is done after in parallel
            if size == 1:
                # create a mesh size function interpolant
                self.fh = lambda p: self.interpolant(p)

                _bbox = self.bbox

                # signed distance function for a rectangle
                def fdd(p):
                    return sdf.drectangle(p, _bbox[0], _bbox[1], _bbox[2], _bbox[3])

                # signed distance function for a cube
                def fdd2(p):
                    return sdf.dblock(
                        p, _bbox[0], _bbox[1], _bbox[2], _bbox[3], _bbox[4], _bbox[5]
                    )

                if _dim == 2:
                    self.fd = lambda p: fdd(p)
                if _dim == 3:
                    self.fd = lambda p: fdd2(p)

        # if parallel
        if size > 1:
            self._construct_lambdas(comm)
        return self

    def _construct_lambdas(self, comm):
        """
        Build lambda fields (for parallel only) for
        SDF and mesh size function
        """
        self.bbox = comm.bcast(self.bbox, 0)

        _dim = self.dim
        _bbox = self.bbox

        # create a mesh size function interpolant
        self.fh = lambda p: self.interpolant(p)

        def fdd(p):
            return sdf.drectangle(p, _bbox[0], _bbox[1], _bbox[2], _bbox[3])

        def fdd2(p):
            return sdf.dblock(
                p, _bbox[0], _bbox[1], _bbox[2], _bbox[3], _bbox[4], _bbox[5]
            )

        # create a signed distance function
        if _dim == 2:
            self.fd = lambda p: fdd(p)
        if _dim == 3:
            self.fd = lambda p: fdd2(p)
        return self

    def plot(self, stride=5, comm=None):
        """ Plot the isotropic mesh size function"""
        comm = comm or MPI.COMM_WORLD
        if comm.rank == 0:
            _dim = self.dim
            _fh = self.fh
            _domain_ext = self.domain_ext
            _width = self.width
            _depth = self.depth

            if _dim == 2:
                zg, xg = self.__CreateDomainMatrices()
                sz1z, sz1x = zg.shape
                _zg = np.reshape(zg, (-1, 1))
                _xg = np.reshape(xg, (-1, 1))
                hh = _fh((_zg, _xg))
                hh = np.reshape(hh, (sz1z, sz1x))

                fig, ax = plt.subplots()
                plt.pcolormesh(
                    xg[0::stride], zg[0::stride], hh[0::stride], edgecolors="none"
                )
                if _domain_ext > 0:
                    rect = plt.Rectangle(
                        (0, -_depth + _domain_ext),
                        _width - 2 * _domain_ext,
                        _depth - _domain_ext,
                        fill=False,
                        edgecolor="black",
                    )
                    ax.add_patch(rect)

                plt.title("Isotropic mesh sizes")
                plt.colorbar(label="mesh size (m)")
                plt.xlabel("x-direction (m)")
                plt.ylabel("z-direction (m)")
                ax.axis("equal")
                # ax.set_xlim(0 - _domain_ext, _width)
                # ax.set_ylim(-_depth, 0)
                plt.show()
            elif _dim == 3:
                print("visualization in 3D not yet supported!")
        return None

    def GetDomainMatrices(self):
        """ Accessor to private method"""
        if self.dim == 2:
            zg, xg = self.__CreateDomainMatrices()
            return zg, xg
        elif self.dim == 3:
            zg, xg, yg = self.__CreateDomainMatrices()
            return zg, xg, yg

    def hj(self, fun, elen, imax):
        """Call-back to the cpp gradient limiter code """
        _dim = self.dim
        _nz = self.nz
        _nx = self.nx
        _ny = 1
        if _dim == 3:
            _ny = self.ny
        _grade = self.grade

        sz = (_nz, _nx, _ny)

        fun = fun.flatten("F")
        tmp = limgrad([*sz], elen, _grade, imax, fun)
        if _dim == 2:
            return np.reshape(tmp, (_nz, _nx), "F")
        if _dim == 3:
            return np.reshape(tmp, (_nz, _nx, _ny), "F")

    def WriteVelocityModel(self, ofname, comm=None):
        """Writes a velocity model as a hdf5 file for use in Spyro """
        comm = comm or MPI.COMM_WORLD
        if comm.rank == 0:
            _dim = self.dim
            _vp = self.velocity_grid
            _nz = self.nz
            _nx = self.nx
            if _dim == 3:
                _ny = self.ny
            _domain_ext = self.domain_ext
            _spacingZ = self.spacingZ
            _spacingX = self.spacingX
            _padstyle = self.padstyle
            nnz = int(_domain_ext / _spacingZ)
            nnx = int(_domain_ext / _spacingX)
            # create domain extension in velocity model
            mx = np.amax(_vp)
            if _dim == 2:
                if _padstyle == "edge":
                    _vp = np.pad(_vp, ((nnz, 0), (nnx, nnx)), "edge")
                elif _padstyle == "constant":
                    # set to maximum value in domain
                    _vp = np.pad(
                        _vp,
                        ((nnz, 0), (nnx, nnx)),
                        "constant",
                        constant_values=(mx, mx),
                    )
                elif _padstyle == "linear_ramp":
                    # linearly ramp to maximum value in domain
                    _vp = np.pad(
                        _vp, ((nnz, 0), (nnx, nnx)), "linear_ramp", end_values=(mx, mx)
                    )
            if _dim == 3:
                if _padstyle == "edge":
                    _vp = np.pad(_vp, ((nnz, 0), (nnx, nnx), (nnx, nnx)), "edge")
                elif _padstyle == "linear_ramp":
                    _vp = np.pad(
                        _vp,
                        ((nnz, 0), (nnx, nnx), (nnx, nnx)),
                        "linear_ramp",
                        end_values=(mx, mx),
                    )
                elif _padstyle == "constant":
                    # set to maximum value in domain
                    _vp = np.pad(
                        _vp,
                        ((nnz, 0), (nnx, nnx), (nnx, nnx)),
                        "constant",
                        constant_values=(mx, mx),
                    )

            print(_vp.shape)
            _nz += nnz  # only bottom
            _nx += nnx * 2  # left and right
            if _dim == 3:
                _ny += nnx * 2  # behind and in front

            ofname += ".hdf5"
            print("Writing velocity model " + ofname, flush=True)
            with h5py.File(ofname, "w") as f:
                f.create_dataset("velocity_model", data=_vp, dtype="f")
                if _dim == 2:
                    f.attrs["shape"] = (_nz, _nx)
                elif _dim == 3:
                    f.attrs["shape"] = (_nz, _nx, _ny)
                f.attrs["units"] = "m/s"

    def SaveMeshSizeFunctionOptions(self, ofname):
        " Save your mesh size function options as a hdf5 file" ""
        ofname += ".hdf5"
        print("Saving mesh size function options " + ofname)
        with h5py.File(ofname, "a") as f:
            f.attrs["bbox"] = self.bbox
            f.attrs["units"] = self.units
            f.attrs["hmin"] = self.hmin
            f.attrs["hmax"] = self.hmax
            f.attrs["wl"] = self.wl
            f.attrs["freq"] = self.freq
            f.attrs["dt"] = self.dt
            f.attrs["cr_max"] = self.cr_max
            f.attrs["grade"] = self.grade
            f.attrs["nx"] = self.nx
            if self.dim == 3:
                f.attrs["ny"] = self.ny
            f.attrs["nz"] = self.nz
            f.attrs["domain_ext"] = self.domain_ext
            f.attrs["padstyle"] = self.padstyle

    # ---PRIVATE METHODS---#

    def __ReadVelocityModel(self):
        """Parses a numpy.ndarray velocity model"""
        # nz, nx, ny and spacingZ, and spacingX
        sd = self.dim
        if sd == 2:
            self.nz, self.nx = self.velocity_grid.shape
        elif sd == 3:
            self.nz, self.nx, self.ny = self.velocity_grid.shape
        self.spacingZ = self.depth / self.nz
        self.spacingX = self.width / self.nx

        if self.__units == "km-s":
            print("Converting model velocities to m-s...")
            self.velocity_grid *= 1e3

    def __CreateDomainVectors(self):
        """ it is what it is """
        _dim = self.dim
        _nz = self.nz
        _nx = self.nx
        if _dim == 3:
            _ny = self.ny
        _domain_ext = self.domain_ext
        _spacingZ = self.spacingZ
        _spacingX = self.spacingX
        _bbox = self.bbox

        # if domain extension is enabled, we augment the domain vectors
        nnz = int(_domain_ext / _spacingZ)
        nnx = int(_domain_ext / _spacingX)
        if _domain_ext > 0:
            _nz += nnz  # only bottom
            _nx += nnx * 2  # left and right
            if _dim == 3:
                _ny += nnx * 2  # behind and in front

        zvec = np.linspace(_bbox[0], _bbox[1], _nz, dtype=np.float32)
        xvec = np.linspace(
            _bbox[2] - _domain_ext, _bbox[3] + _domain_ext, _nx, dtype=np.float32
        )

        if _dim == 2:
            return zvec, xvec
        elif _dim == 3:
            yvec = np.linspace(
                _bbox[4] - _domain_ext, _bbox[5] + _domain_ext, _ny, dtype=np.float32
            )
            return zvec, xvec, yvec

    def __CreateDomainMatrices(self):
        _dim = self.dim
        if _dim == 2:
            zvec, xvec = self.__CreateDomainVectors()
            zg, xg = np.meshgrid(zvec, xvec, indexing="ij")
            return zg, xg
        elif _dim == 3:
            zvec, xvec, yvec = self.__CreateDomainVectors()
            xg, yg, zg = np.meshgrid(
                xvec, yvec, zvec, indexing="ij", sparse=True, copy=False
            )
            return zg, xg, yg

    def __EditMeshSizeFunction(self, hh_m, rank):
        """ Edits sizing function to support domain extension of variable width """
        _domain_ext = self.domain_ext
        _dim = self.dim
        _hmax = self.hmax
        _spacingZ = self.spacingZ
        _spacingX = self.spacingX
        _padstyle = self.padstyle

        nnz = int(_domain_ext / _spacingZ)
        nnx = int(_domain_ext / _spacingX)

        if rank == 0:
            print(
                "Including a " + str(_domain_ext) + " meter domain extension...",
                flush=True,
            )
            print(
                "Using the padstyle " + _padstyle + " to extend mesh resolution.",
                flush=True,
            )
        if _dim == 2:
            mx = np.amax(hh_m)
            if _padstyle == "edge":
                hh_m = np.pad(hh_m, ((nnz, 0), (nnx, nnx)), "edge")
            elif _padstyle == "constant":
                # set to maximum value in domain
                hh_m = np.pad(
                    hh_m, ((nnz, 0), (nnx, nnx)), "constant", constant_values=(mx, mx)
                )
            elif _padstyle == "linear_ramp":
                # linearly ramp to maximum value in domain
                hh_m = np.pad(
                    hh_m, ((nnz, 0), (nnx, nnx)), "linear_ramp", end_values=(mx, mx)
                )
            hh_m = np.where(hh_m > _hmax, _hmax, hh_m)
            return hh_m
        if _dim == 3:
            mx = np.amax(hh_m)
            if _padstyle == "edge":
                hh_m = np.pad(hh_m, ((nnz, 0), (nnx, nnx), (nnx, nnx)), "edge")
            elif _padstyle == "linear_ramp":
                # linearly ramp to maximum value in domain
                hh_m = np.pad(
                    hh_m,
                    ((nnz, 0), (nnx, nnx), (nnx, nnx)),
                    "linear_ramp",
                    end_values=(mx, mx),
                )
            else:
                raise ValueError("3D pad style currently not supported yet")
            hh_m = np.where(hh_m > _hmax, _hmax, hh_m)
            return hh_m

    def __CreateDomainExtension(self):
        """Edit bbox to reflect domain extension"""
        _dim = self.dim
        _bbox = self.bbox
        _domain_ext = self.domain_ext
        if _domain_ext > 0:
            if _dim == 2:
                bbox_new = (
                    _bbox[0] - _domain_ext,
                    _bbox[1],
                    _bbox[2] - _domain_ext,
                    _bbox[3] + _domain_ext,
                )
            if _dim == 3:
                bbox_new = (
                    _bbox[0] - _domain_ext,
                    _bbox[1],
                    _bbox[2] - _domain_ext,
                    _bbox[3] + _domain_ext,
                    _bbox[4] - _domain_ext,
                    _bbox[5] + _domain_ext,
                )

            self.bbox = bbox_new
            self.width = bbox_new[3] - bbox_new[2]
            self.depth = bbox_new[1] - bbox_new[0]
            if _dim == 3:
                self.length = bbox_new[5] - bbox_new[4]
        return self
