Modules
==============

Here we document the public API

*SeimsicMesh.Geometry*
-------------------------------

Routines to perform geometrical/topological operations and calculate things on meshes.

.. automodule:: SeismicMesh.geometry
    :members:

*SeimsicMesh.get_sizing_function_from_segy*
---------------------------------------------

Function to build a :math:`f(h)` mesh sizing function from a seismic velocity model.
Assumes the domain can be represented by a rectange (2D) or cube (3D) and thus builds a :math:`f(d)` accordingly.
The sizing function can be passed to `generate_mesh` to build a simplical mesh with varying mesh resolution.

.. automodule:: SeismicMesh.get_sizing_function_from_segy
    :members:

*SeimsicMesh.generate_mesh*
-------------------------------

Function to build a simplical mesh that conforms to the signed distance function :math:`f(d)`
and :math:`f(h)`.

.. automodule:: SeismicMesh.generate_mesh
    :members:
