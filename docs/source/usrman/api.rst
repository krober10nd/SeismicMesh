Modules
==============

Here we document the public API (e.g., functions, classes and their methods)

*SeimsicMesh.Geometry*
-------------------------------

Routines to perform geometrical/topological operations and calculate things on meshes

.. automodule:: SeismicMesh.geometry
    :members:

*SeimsicMesh.MeshSizeFunction*
-------------------------------

Convenience class to build a :math:`f(h)` mesh sizing function from a velocity model.
This class can be passed to :class:`MeshGenerator` to build a simplical mesh.

.. autoclass:: SeismicMesh.MeshSizeFunction
    :members:

*SeimsicMesh.MeshGenerator*
-------------------------------

Class to build a simplical mesh that conforms to the signed distance function :math:`f(d)`
and :math:`f(h)`.

.. autoclass:: SeismicMesh.MeshGenerator
    :members:
