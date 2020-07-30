---
title: 'SeismicMesh: Triangular meshing for seismology'
tags:
  - Python
  - Seismology
  - Full Waveform Inversion
  - Parallel mesh generation
authors:
  - name: Keith Jared Roberts
    affiliation: 1
affiliations:
 - name: Escola Politécnica da Universidade de São Paulo
   index: 1
date: 30 July 2020
bibliography: paper.bib

---

# Summary

`SeismicMesh` is a Python package for quality two and three dimensional triangular and tetrehedral
meshing directly from seismic velocity models using signed distance functions. Low-level C++ code
is wrapped with Python for high performance at scale without losing flexibility or ease-of-use in
the user-interface. At a low level, the program makes direct calls to the Computational Geometry
Algorithms Library [@cgal:hs-chdt3-20a] to robustly and efficiently produce multiscale meshes.
For large scale problems the computationally expensive operations
are parallelized and mesh generation is accomplished through combined modifications to the
DistMesh alogrithm [@doi:10.1137/S0036144503429121] and a parallel Delaunay triangulation algorithm [@peterka2014high].
Much like other popular mesh generators `gmsh` [doi:10.1002/nme.2579] and `tetgen` [@si2015tetgen], mesh generation can be made fully scriptable.
However, there are two important and unique aspects of this work: first we provide an ability to generate graded mesh sizing functions that
control distribution of mesh sizes directly from seismic P-wave velocity models.
In this way, mesh resolution respects material properties and numerical stability can be ensured with little effort on part of the user.
Secondly, mesh geometry is defined using set operations with signed distance functions thus requiring no
manual tracing or use of graphical user interfaces to mesh irregular 2D/3D geometrical features.


# Background

Generating a high-quality graded mesh for a geophysical domain represents a modern
challenge for sophiscated geophysical simulation workflows. In these applications,
a domain is discretized typically with simplicial elements (e.g., triangles/tetrahedrals)
that adapt in size to features of interest. These meshes are commonly used with Finite Element
methods to solve Partial Differential Equation that model phyiscal processes such as
the acoustic or elastic wave equation. These simulations are often used in geophysical exploration
studies to solve inverse problems such as Full Waveform Inversion (FWI) [doi:10.1190/1.1441754; @virieux2009overview]
and Reverse Time Migration (RTM) [@10.1093/gji/ggv380].

# Statement of need

There are many aspects to consider when building a mesh for a geophysical inverse problem
such as Full Waveform Inversion, which make mesh generation for geophysical domains quite complex.
For instacne, mesh elements must be sufficiently well-shaped, sized to maintain numerical stability, minimize numerical dissipatation,
and maximize physical and numerical accuracy. For applications such as FWI and RTM, material discontinuties
need to be well represented to ensure reflections are as accurate as possible. In cases
with complex and irregular rock structures, explicit geometries may not exist complicating
mesh generation workflows and requiring the use of graphical user interfaces to create these geometries.


# Geometry treatment

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Fenced code blocks are rendered with syntax highlighting:
```python
for n in range(10):
    yield f(n)
```

# Acknowledgements

Shell and Research Center for Gas Innovation

# References
