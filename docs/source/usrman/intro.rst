Introduction
============

Generating a high-quality simplicial mesh for a geophysical
domain represents a modern challenge for sophiscated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (triangles/tetrahedrals)
that adapt in size to features of interest, while considering numerical
aspects. These meshes are commonly used with Finite Element and Finite volume methods to solve
Partial Differential Equations (PDEs) that represent phyiscal processes.

Depending on the nature of the geophysical problems, mesh elements must be well-shaped,
sufficiently sized to maintain numerical stability, minimize numerical dissipatation, and maximize accuracy.
Further, for geophysical domains in which explicit geometries of structures do not exist, mesh generation workflows
tend to become very manual requiring the use of hand-tracing and graphical user interfaces. Unfortuantely,
this largely compromises the reproducibility of mesh connectivity and mesh generation procedure.




Why *another* mesh generator?
------------





