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

from .size_function import SizeFunction

from .cpp import limgrad

__all__ = [
    "get_sizing_function_from_segy",
    "write_velocity_model",
    "plot_sizing_function",
]

# reasonable sizing function options go here
opts = {
    "hmin": 150.0,
    "hmax": 10000.0,
    "wl": 0,
    "freq": 2.0,
    "grad": 0.0,
    "grade": 0.0,
    "space_order": 1,
    "dt": 0.0,
    "cr_max": 1.0,
    "pad_style": "edge",
    "domain_pad": 0.0,
    "units": "m-s",
    "nz": None,
    "nx": None,
    "ny": None,
    "byte_order": "byte_order",
}


def get_sizing_function_from_segy(filename, bbox, comm=None, **kwargs):
    r"""Build a mesh size function from a seismic velocity model.

    :param filename:
        The name of a SEG-y or binary file containing a seismic velocity model
    :type filename: ``string``
    :param bbox:
        Bounding box containing domain extents of the velocity model contained in `filename`.
    :type bbox: `tuple` with size (2*dim). For example, in 2D `(zmin, zmax, xmin, xmax)`
    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *hmin* (``float``) --
            Minimum element size in the domain (default==150 m)
        * *hmax* (``float``) --
            Maximum element size in the domain (default==10,000 m)
        * *wl* (``int``) --
            Number of vertices per wavelength for a given 𝑓𝑚𝑎𝑥 (default==0 vertices)
        * *freq* (``float``) --
            𝑓𝑚𝑎𝑥 in hertz for which to estimate `wl` (default==2 Hertz)
        * *grad* (``float``) --
            Resolution in meters nearby sharp gradients in velociy (default==0 m)
        * *grade* (``float``) --
            Maximum allowable variation in mesh size in decimal percent (default==0.0)
        * *space_order* (``int``) --
            Simulation will be attempted with a mesh using the polynomial order `space_order` of the basis functions (default==1)
        * *dt* (``float``) --
            Theoretical maximum stable timestep in seconds given Courant number Cr (default==0.0 s)
        * *cr_max* (``float``) --
            The mesh simulated with this `dt` has this maximum Courant number (default==1.0)
        * *pad_style* (``string``) --
             The method (`edge`, `linear_ramp`, `constant`) to pad velocity in the domain pad region (default==None)
        * *domain_pad* (``float``) --
             The width of the domain pad in -z, +x, -x, +y, -y directions (default==0.0 m).
        * *units* (``string``) --
             The units of the seismic velocity model (default='m-s')
        * *nz* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the z-direction in the velocity model.
        * *ny* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the y-direction in the velocity model.
        * *nx* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the x-direction in the velocity model.
        * *byte_order* (``string``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The order of bytes in a 3D sesimic velocity model (`little` or `big`).

    :return: a :class:`SizeFunction` object with a `obj.bbox` field and an `obj.eval` method.
    :rtype: a :class:`SizeFunction` object

    """
    comm = comm or MPI.COMM_WORLD
    cell_size = None

    if comm.rank == 0:
        opts.update(kwargs)

        vp, nz, nx, ny = _read_velocity_model(
            filename=filename,
            nz=opts["nz"],
            nx=opts["nx"],
            ny=opts["ny"],
            byte_order=opts["byte_order"],
        )

        if opts["units"] == "km-s":
            vp *= 1000.0

        # check the bbox
        if len(bbox) == 4:
            dim = 2
        elif len(bbox) == 6:
            dim = 3
        else:
            raise ValueError("Dimension not supported")

        cell_size = _initialize_sizing_function(dim, opts["hmin"], nz, nx, ny)

        for key in kwargs:
            if key in {
                "hmin",
                "hmax",
                "wl",
                "freq",
                "cr_max",
                "dt",
                "space_order",
                "grad",
                "grade",
                "pad_style",
                "domain_pad",
                "units",
                "nz",
                "nx",
                "ny",
                "byte_order",
            }:
                pass
            else:
                raise ValueError(
                    "Option %s with parameter %s not recognized " % (key, kwargs[key])
                )
        if np.any([opts["wl"] > 0, opts["grad"] > 0]):
            cell_size = np.minimum(
                _wavelength_sizing(vp, opts["wl"], opts["freq"]),
                _gradient_sizing(vp, opts["grad"]),
            )

        print("Enforcing minimum element size of " + str(opts["hmin"]))
        cell_size[cell_size < opts["hmin"]] = opts["hmin"]

        print("Enforcing maximum element size of " + str(opts["hmax"]))
        cell_size[cell_size > opts["hmax"]] = opts["hmax"]

        cell_size = _enforce_courant_sizing(
            vp, cell_size, opts["cr_max"], opts["dt"], opts["space_order"]
        )

        cell_size = _enforce_gradation_sizing(
            cell_size, opts["grade"], (bbox[3] - bbox[2]) / nx
        )

        cell_size, vp, bbox = _build_domain_pad(cell_size, vp, bbox, opts)

        fh = _build_sizing_function(cell_size, vp, bbox)
    else:

        def fh(p):
            return 1

    # agreement re the bbox
    if comm.size > 1:
        bbox = comm.bcast(bbox, 0)

    return SizeFunction(bbox, fh)


def write_velocity_model(filename, ofname=None, comm=None, **kwargs):
    r"""Reads and then writes a velocity model as a hdf5 file

    :param filename: filename of a velocity model (either .segy or .bin)
    :type filename: `string`
    :param ofname: filename of output hdf5 file
    :type ofname: `string`, optional
    :param comm: MPI communicator
    :type comm: MPI4py communicator, optional

    :param \**kwargs:
        See below

    :Keyword Arguments:
        * *nz* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the z-direction in the velocity model.
        * *ny* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the y-direction in the velocity model.
        * *nx* (``int``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The number of grid points in the x-direction in the velocity model.
        * *byte_order* (``string``) --
             REQUIRED FOR BINARY VELOCITY MODEL. The order of bytes in a 3D sesimic velocity model (`little` or `big`).
        * *bbox* (``tuple``) --
             Coordinates of the velocity model's domain extents. Only required if padding the domain.
        * *domain_pad* (``float``) --
             Width of the domain pad in meters.
        * *pad_style* (``string``) --
             Type of padding.


    """
    comm = comm or MPI.COMM_WORLD
    if comm.rank == 0:
        opts.update(kwargs)

        if ofname is None:
            warnings.warn("No output filename specified, name will be `filename`")
            ofname = filename

        vp, nz, nx, ny = _read_velocity_model(
            filename=filename,
            nz=opts["nz"],
            nx=opts["nx"],
            ny=opts["ny"],
            byte_order=opts["byte_order"],
        )

        if opts["domain_pad"] > 0.0:
            print("Adding a domain pad to the velocity model...")
            tmp = vp
            _, vp, _ = _build_domain_pad(tmp, vp, opts["bbox"], opts)

        ofname += ".hdf5"
        print("Writing velocity model: " + ofname, flush=True)
        with h5py.File(ofname, "w") as f:
            f.create_dataset("velocity_model", data=vp, dtype="f")
            f.attrs["shape"] = vp.shape
            f.attrs["units"] = "m/s"


def plot_sizing_function(cell_size, stride=1, comm=None):
    """Plot the mesh size function in 2D

    :param cell_size: a callable function that takes a point and gives a size
    :type cell_size: a callable function object
    :param stride: skip `stride` points to save on memory when plotting
    :type stride: `int`, optional
    :param comm: MPI communicator
    :type comm: MPI4py communicator, optional

    """
    comm = comm or MPI.COMM_WORLD
    if comm.rank == 0:

        if not isinstance(cell_size, SizeFunction):
            raise ValueError("Can only plot a :class:`SizeFunction`")

        bbox = cell_size.bbox
        fh = cell_size.eval

        if len(bbox) != 4:
            raise ValueError("Visualization in 3D not supported")

        zg, xg = np.meshgrid(
            np.arange(bbox[0], bbox[1], 50.0),
            np.arange(bbox[2], bbox[3], 50.0),
            indexing="ij",
        )
        cell_size = fh((zg, xg))

        fig, ax = plt.subplots()
        plt.pcolormesh(
            xg[0::stride], zg[0::stride], cell_size[0::stride], shading="auto"
        )
        plt.title("Isotropic mesh sizes")
        plt.colorbar(label="mesh size (m)")
        plt.xlabel("x-direction (m)")
        plt.ylabel("z-direction (m)")
        ax.axis("equal")
        plt.show()
    return ax


def _build_sizing_function(cell_size, vp, bbox):
    """Builds a regular gridded interpolant to query during mesh generation"""
    dim = cell_size.ndim
    if dim == 2:
        nz, nx, dz, dx = _get_dimensions(vp, bbox)
        z_vec, x_vec = _get_vectors(dim, bbox, nz, nx)
        interpolant = RegularGridInterpolator(
            (z_vec, x_vec), cell_size, bounds_error=False, fill_value=None
        )
    elif dim == 3:
        nz, nx, ny, dz, dx, dy = _get_dimensions(vp, bbox)
        z_vec, x_vec, y_vec = _get_vectors(dim, bbox, nz, nx, ny)
        interpolant = RegularGridInterpolator(
            (z_vec, x_vec, y_vec), cell_size, bounds_error=False, fill_value=None
        )
    else:
        raise ValueError("Dimension not supported")
    return interpolant


def _wavelength_sizing(vp, wl=5, freq=2.0):
    """Mesh sizes are distributed according to an estimate of the wavelength
    of an acoustic/elastic wave.
    """
    if wl == 0.0:
        return 99999
    if wl < 0:
        raise ValueError("Parameter `wl` must be set > 0")
    if freq < 0.0:
        raise ValueError("Parameter `freq` must be set > 0.0")
    print(
        "Mesh sizes will be built to resolve an estimate of wavelength of a "
        + str(freq)
        + " hz wavelet with "
        + str(wl)
        + " vertices...",
        flush=True,
    )
    return vp / (freq * wl)


def _gradient_sizing(vp, grad):
    """Refine the mesh near sharp gradients in seismic velocity."""
    if grad == 0.0:
        return 99999
    print("Refining mesh sizes near sharp velocity gradients...")
    if grad < 0:
        raise ValueError("Parameter grad must be > 0")

    window = [100] * vp.ndim
    win_mean = ndimage.uniform_filter(vp, tuple(window))
    win_sqr_mean = ndimage.uniform_filter(vp ** 2, tuple(window))
    win_var = win_sqr_mean - win_mean ** 2
    # normalize variance to [0,1]
    win_var /= np.amax(win_var)
    win_var -= np.amin(win_var)
    return grad / (win_var + 0.10)


def _enforce_courant_sizing(vp, cell_size, cr_max, dt, space_order):
    """Ensure mesh resolution distribution doesn't violate CFL (cr_max > 1.0)"""
    if (cr_max == 0.0) or (dt == 0.0) or (space_order == 0.0):
        return cell_size
    print("Enforcing timestep of " + str(dt) + " seconds...", flush=True)
    if cr_max < 0:
        raise ValueError("Parameter `cr_max` must be > 0.0")
    if dt < 0:
        raise ValueError("Parameter `dt` must be > 0.0")
    if space_order < 1:
        raise ValueError("Parameter `space_order` must be >= 1 ")
    dim = vp.ndim
    cr_old = (vp * dt) / (dim * cell_size)
    cr_max = cr_max / (dim * space_order)
    dxn = (vp * dt) / (dim * cr_max)
    return np.where(cr_old > cr_max, dxn, cell_size)


def _enforce_gradation_sizing(cell_size, grade, elen):
    """Call-back to the cpp gradient limiter code """
    if grade == 0.0:
        warnings.warn(
            "Mesh size gradient is deactiavted. This may compromise mesh quality"
        )
        return cell_size
    if grade < 0:
        raise ValueError("Parameter `grade` must be > 0.0")
    if grade > 1.0:
        warnings.warn("Parameter `grade` is set pretty high (> 1.0)!")
    print("Enforcing mesh size gradation of " + str(grade) + " decimal percent...")

    dim = cell_size.ndim
    sz = cell_size.shape
    if len(sz) == 2:
        sz = (sz[0], sz[1], 1)

    cell_size = cell_size.flatten("F")
    tmp = limgrad([*sz], elen, grade, 10000, cell_size)
    if dim == 2:
        return np.reshape(tmp, (sz[0], sz[1]), "F")
    elif dim == 3:
        return np.reshape(tmp, (sz[0], sz[1], sz[2]), "F")
    else:
        raise ValueError("Dimension not supported")


def _get_dimensions(vp, bbox):
    dim = vp.ndim
    if dim == 2:
        nz, nx = vp.shape
        dz = int((bbox[1] - bbox[0]) / nz)
        dx = int((bbox[3] - bbox[2]) / nx)
        return nz, nx, dz, dx
    elif dim == 3:
        nz, nx, ny = vp.shape
        dz = int((bbox[1] - bbox[0]) / nz)
        dx = int((bbox[3] - bbox[2]) / nx)
        dy = int((bbox[5] - bbox[4]) / ny)
        return nz, nx, ny, dz, dx, dy
    else:
        raise ValueError("Dimension not supported")


def _get_vectors(dim, bbox, nz, nx, ny=None):
    zvec = np.linspace(bbox[0], bbox[1], nz, dtype=np.float32)
    xvec = np.linspace(bbox[2], bbox[3], nx, dtype=np.float32)
    if dim == 2:
        return zvec, xvec
    elif dim == 3:
        yvec = np.linspace(bbox[4], bbox[5], ny, dtype=np.float32)
        return zvec, xvec, yvec
    else:
        raise ValueError("Dimension not supported")


def _build_domain_pad(cell_size, vp, bbox, opts):
    """Building a domain extension"""
    dim = vp.ndim
    domain_pad = opts["domain_pad"]
    pad_style = opts["pad_style"]
    if domain_pad < 0:
        raise ValueError("Domain extension must be >= 0")

    if domain_pad > 0:
        print("Including a " + str(domain_pad) + " meter domain extension...")
        if dim == 2:
            nz, nx, dz, dx = _get_dimensions(vp, bbox)
            nnz = int(domain_pad / dz)
            nnx = int(domain_pad / dx)
            bbox = (
                bbox[0] - domain_pad,
                bbox[1],
                bbox[2] - domain_pad,
                bbox[3] + domain_pad,
            )
        elif dim == 3:
            nz, nx, ny, dz, dx, dy = _get_dimensions(vp, bbox)
            nnz = int(domain_pad / dz)
            nnx = int(domain_pad / dx)
            nny = int(domain_pad / dy)
            bbox = (
                bbox[0] - domain_pad,
                bbox[1],
                bbox[2] - domain_pad,
                bbox[3] + domain_pad,
                bbox[4] - domain_pad,
                bbox[5] + domain_pad,
            )

        print("Using the pad_style: " + pad_style)
        if dim == 2:
            padding = ((nnz, 0), (nnx, nnx))
        elif dim == 3:
            padding = ((nnz, 0), (nnx, nnx), (nny, nny))

        max_cell_size = np.amax(cell_size)
        max_vp = np.amax(vp)

        cell_size = _pad_it(cell_size, padding, pad_style, [max_cell_size] * 2)
        vp = _pad_it(vp, padding, pad_style, [max_vp] * 2)

    return cell_size, vp, bbox


def _pad_it(array, padding, style, extra):
    """Add a domain extension to `cell_size` and `vp`"""
    if style == "edge":
        array = np.pad(array, padding, "edge")
    elif style == "constant":
        array = np.pad(array, padding, "constant", constant_values=tuple(extra))
    elif style == "linear_ramp":
        array = np.pad(array, padding, "linear_ramp", end_values=tuple(extra))
    else:
        raise ValueError(
            "pad style currently not supported. Try `linear_ramp`, `edge`, or `constant`"
        )
    return array


def _read_velocity_model(filename, nz=None, nx=None, ny=None, byte_order=None):
    """Read a velocity model"""
    if filename.endswith(".segy"):
        return _read_segy(filename)
    else:
        return _read_bin(filename, nz, nx, ny, byte_order)


def _read_bin(filename, nz, nx, ny, byte_order):
    """Read a velocity model from a binary"""
    if (nz is None) or (nx is None) or (ny is None):
        raise ValueError(
            "Please specify the number of grid points in each dimension (e.g., `nz`, `nx`, `ny`)..."
        )
    with open(filename, "r") as file:
        print("Reading binary file: " + filename)
        if byte_order == "big":
            vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder(">"))
        elif byte_order == "little":
            vp = np.fromfile(file, dtype=np.dtype("float32").newbyteorder("<"))
        else:
            raise ValueError("Please specify byte_order as either: little or big.")
        vp = vp.reshape(nx, ny, nz, order="F")
        return np.flipud(vp.transpose((2, 0, 1))), nz, nx, ny  # z, x and then y


def _read_segy(filename):
    """Read a velocity model from a SEG-y file"""
    with segyio.open(filename, ignore_geometry=True) as f:
        nz, nx = len(f.samples), len(f.trace)
        vp = np.zeros(shape=(nz, nx))
        for index, trace in enumerate(f.trace):
            vp[:, index] = trace
        if np.amin(vp) < 1000.0:
            raise Warning(
                "Velocity appear to be in km/s. Maybe pass `units` km-s key pair?"
            )
        return np.flipud(vp), nz, nx, 0


def _initialize_sizing_function(dim, hmin, nz, nx, ny=None):
    """initialize a sizing function grid"""
    if dim == 2:
        cell_size = np.full((nz, nx), hmin, dtype=float)
    elif dim == 3:
        cell_size = np.full((nz, nx, ny), hmin, dtype=float)
    else:
        raise ValueError("Dimension not supported")
    return cell_size
