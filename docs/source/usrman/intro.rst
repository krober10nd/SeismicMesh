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

While there are a vast array of mesh generation tools in existince, most tools are general purpose and simply produce a set of points and triangles while requriing the user to define their own sizing distribution. Mesh generation represents only one stage of a typical geophysical modeling workflow. First, the integration of geophysical data is often required to guide the appliation of mesh resolution to ensure accurate and efficient numerical simulations.
These critical pre-processing steps should ideally be script-able and included in the software stack for automation and reproduction of the model; however, this then requires the user write their own ad hoc scripts. Additionally, geometry creation becomes non-trivial when irregular geometry structures are present. We argue this needs to be as objective as possible and not rely on manual tracing/drawing of segments/surfaces.

As a result, we've implemented a "batteries included" type approach to mesh generation for seismology. It  contains all the tools necessary to go from a raw geophsyical seismic velocity model to a unstructured triangular mesh. We include an approach to build isotropic mesh sizing functions that lead to high-quality, graded meshes which are numerically stable. For geometry creation, we rely on signed distance functions to define features implictly, which avoids the need to supply explicit locations of segments/surfaces and thus simplifies workflows. Additionally, methods are provided to convert existing meshes into signed distance functions to perform adapation. 
