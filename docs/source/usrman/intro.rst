Introduction
============

Generating a high-quality graded mesh for a geophysical domain represents a modern challenge for sophiscated geophysical simulation workflows.
In these applications, a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element and Finite volume methods to solve
Partial Differential Equations (PDEs) that model phyiscal processes such as the acoustic or elastic wave equation. These kind of simulations are
often used in geophysical exploration studies to solve inverse problems such as Full Waveform Inversion.

There are many aspects to consider when building a mesh for a geophysical inverse problem. Mesh elements must be sufficiently well-shaped,
sized to maintain numerical stability, minimize numerical dissipatation, and maximize physical and numerical accuracy. For applications such as Travel Time Tomography and Reverse Time Migration, material discontinuties need to be well represented to ensure reflections are as accurate as possible. In cases with complex and irregular rock structures, explicit geometries may not exist complicating mesh generation workflows and requiring the use of graphical user interfaces.

Why *another* mesh generator?
-------------------------------

While there are an array of mesh generation tools in existince, most are very general purpose and simply produce a set of points and triangles while requriing the user to define their own element sizing distribution. It's also important to stress, mesh generation represents only one stage of a typical geophysical modeling workflow. First, the integration of geophysical data is often required to guide the workflow and design of mesh resolution to ensure accurate and efficient numerical simulations. These critical pre-processing steps should ideally be script-able and included in the software stack for automation and reproduction of the model. However, this requires the user write their own ad hoc scripts to this end. Additionally, geometry creation becomes non-trivial when irregular geometry structures are present in the domain. Irregular geometries common to seimsology such as faults, salt-bodies, and the like aren't necessarily easily to handle or constrain in 3D meshes explictly. Our contention is geometry creation needs to be as objective as possible and not rely on manual tracing/drawing of segments/surfaces, while integrating together with mesh sizing capabilities.

As a result, we've implemented a "batteries included" type approach to mesh generation for seismological domains. It contains all the tools necessary to go from a raw geophsyical seismic velocity model to a unstructured triangular mesh for both 2D and 3D domains. We include an approach to build isotropic mesh sizing functions that lead to high-quality, graded meshes which are shown to be numerically stable. For geometry creation, we rely on signed distance functions to define features implictly, which avoids the need to supply explicit locations of segments/surfaces and thus hopefully simplifies workflows. Additionally, methods are provided to build signed distance functions directly from seismic velocity models.
