Introduction
============

Generating a high-quality graded mesh for a geophysical domain represents a modern challenge for sophiscated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element and Finite volume methods to solve
Partial Differential Equations (PDEs) that model phyiscal processes such as the acoustic or elastic wave equation. These kind of simulations are
often used in geophysical exploration studies to solve inverse problems such as Full Waveform Inversion.

There are many aspects to consider when building a mesh for a geophysical inverse problem. Mesh elements must be sufficiently well-shaped,
sized to maintain numerical stability, minimize numerical dissipatation, and maximize physical and numerical accuracy. For applications such as Travel Time Tomography and Reverse Time Migration, material discontinuties need to be well represented to ensure reflections are as accurate as possible. In cases with complex and irregular rock structures, explicit geometries may not exist complicating mesh generation workflows and requiring the use of graphical user interfaces to create these geometries.

The purpose of this software is quickly assemble meshes for geophysical simulation using the Finite Element method, while giving users simple controls to design their own models without having to learn or write much code. This program strives for automation and reproduction of the model. As a result, I've implemented a "batteries included" type approach to mesh generation for seismological domains. It contains the technology necessary to go from a raw geophsyical seismic velocity model to a unstructured triangular mesh for both 2D and 3D domains in a scriptable manner. An approach to build isotropic mesh sizing functions that lead to high-quality, graded meshes with the generator is detailed. For geometry creation, I rely on signed distance functions to define domain features implictly, which avoids the need to supply explicit locations of segments/surfaces. Additionally, methods are provided to build signed distance functions directly from isocontours extracted from seismic velocity models.
