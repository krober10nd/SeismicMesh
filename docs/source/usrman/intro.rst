Introduction
============

Generating a high-quality graded mesh for a geophysical domain represents a modern challenge for sophiscated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element and Finite volume methods to solve
Partial Differential Equations (PDEs) that model phyiscal processes such as the acoustic or elastic wave equation. These kind of simulations are
often used in geophysical exploration studies to solve inverse problems.

There are many aspects to consider when building a mesh for a geophysical problem.
Depending on the nature of the geophysical problem, mesh elements must be sufficiently well-shaped,
sized to maintain numerical stability, minimize numerical dissipatation, and maximize physical and numerical accuracy.
Further, for geophysical domains in which explicit geometries of structures do not exist, mesh generation workflows
can often become manual and require the use of graphical user interfaces, which I believe, hinders automation and reproducibility
of simulation results.

Why *another* mesh generator?
-------------------------------

While there are a vast array of mesh generation tools in existince, most mesh generation tools are stand alone and general purpose simply producing a set of points and triangles. However, the mesh generation step represents only one stage of a typical geophysical modeling workflow. For example, geophysical data is often required to guide the appliation of resolution in an objective and reproducible manner and to ensure accurate numerical simulations.
These pre-processing steps should ideally be script-able and included in the software stack for automation and reproduction of the mesh.
Geometry creation should also be as objective as possible and not rely on ad hoc manual specification of segments/surfaces.

As a result, we implemented a "batteries included" type approach that contains all the tools necessary to build isotropic mesh sizing functions that lead to high-quality, graded meshes with large variations in element sizes that are numerically stable. We rely on signed distance functions to define geometry, which avoids the need to supply explicit locations of segments/surfaces and simplifies workflows. Additionally methods are provided to convert existing meshes into signed distance functions to adapt existing meshes with the software.
