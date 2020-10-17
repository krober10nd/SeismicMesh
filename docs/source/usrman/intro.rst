Introduction
============

Generating a high-quality graded mesh for a geophysical domain represents a modern challenge for sophisticated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element and Finite Volume methods to solve
Partial Differential Equations (PDEs) that model physical processes such as the acoustic or elastic wave equation. These kind of simulations are
often used in geophysical exploration studies to solve inverse problems such as Full Waveform Inversion.

There are many aspects to consider when building a mesh for a geophysical inverse problem. Mesh elements must be sufficiently well-shaped,
sized to maintain numerical stability, minimize numerical dissipation, and maximize physical and numerical accuracy. For applications such as Travel Time Tomography and Reverse Time Migration, material discontinuities need to be well represented to ensure reflections are as accurate as possible. In cases with complex and irregular rock structures, explicit geometries may not exist complicating mesh generation workflows and requiring the use of graphical user interfaces to create these geometries.

The purpose of this software is quickly assemble meshes for geophysical simulation using the Finite Element/Finite Volume methods. This is accomplished by giving users simple controls to design their own meshes without having to learn or write much code or more complex software. This program strives for automation and reproduction of the mesh.

This package contains the technology to build a 2D/3D mesh from a seismic velocity model in a scriptable manner. An approach to build mesh sizing functions that can lead to high-quality, graded meshes with the generator is detailed. For geometry creation, we rely on signed distance functions to define domain features implicitly, which avoids the need to supply explicit locations of segments/surfaces.
