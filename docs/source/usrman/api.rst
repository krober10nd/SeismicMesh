Modules
==============

Here we document the public API

*SeimsicMesh.Geometry*
-------------------------------

Routines to perform geometrical/topological operations and calculate things on meshes.

.. automodule:: SeismicMesh.geometry
    :members:

*SeimsicMesh.sizing*
---------------------------------------------

Function to build a :math:`f(h)` mesh sizing function from a seismic velocity model.
Assumes the domain can be represented by a rectange (2D) or cube (3D) and thus builds a :math:`f(d)` accordingly.

.. automodule:: SeismicMesh.sizing
    :members:

*SeimsicMesh.generation*
-------------------------------

Functions to build and improve a simplical mesh that conforms to the signed distance function :math:`f(d)`
and :math:`f(h)`.

.. automodule:: SeismicMesh.generation
    :members:
