---
title: 'SeismicMesh: Triangular meshing for seismology'
tags:
  - Python
  - Seismology
  - Mesh generation
  - DistMesh
  - Parallel mesh generation
authors:
  - name: Keith J. Roberts
    affiliation: 1
  - name: Rafael dos Santos Gioria
    affiliation: 1
  - name: William J. Pringle
    affiliation: 2
affiliations:
 - name: Research Center for Gas Innovation, Escola Politécnica da Universidade de São Paulo, São Paulo, Brazil.
   index: 1
 - name:  Dept. of Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, 156 Fitzpatrick Hall, Notre Dame, IN, U.S.A.
   index: 2
date: 20 August 2020
bibliography: paper.bib

---
# Summary

`SeismicMesh` is a Python package for generating quality two and three dimensional triangular and tetrahedral meshes from seismic velocity models using signed distance functions to define geometry. C++ code is wrapped with Python for high performance without losing ease-of-use. Mesh generation is performed in either serial or using distributed memory parallelism through modifications to a mesh generation algorithm known as `DistMesh` [@doi:10.1137/S0036144503429121; @peterka2014high]. The program uses the Computational Geometry Algorithms Library [@cgal:hs-chdt3-20a] to perform many geometric operations efficiently.


# Background

Generating a high-geometric quality graded mesh for a geophysical domain represents a modern challenge for seismological modeling using the Finite Element Method (FEM). In these applications, a domain is discretized typically with triangles/tetrahedral elements that vary widely in size around features of interest. In seismology, these meshes are commonly used with the FEM to solve partial differential equations that model acoustic or elastic waves, which are used in seismic velocity model building algorithms such as Full Waveform Inversion (FWI) [@doi:10.1190/1.1441754; @virieux2009overview] and Reverse Time Migration [@10.1093/gji/ggv380].

# Statement of Need

Self-contained automatic meshing workflows that incorporate geophysical data to generate unstructured meshes have been successful in several domains such as coastal ocean modeling [e.g., `OceanMesh2D` @roberts2019oceanmesh2d; `GMT` and `Terreno` @gorman2008systematic] and reservoir modeling [e.g., `MeshIT` @cacace2015meshit]. Besides generating a high-quality mesh, the distribution of elements must accurately represent physical features in the domain such as faults and non-conformities. Element distributions must be placed such that the stability requirements of the numerical solver are satisfied and feature smooth transitions in element size. These aspects motivate the development of a self-contained mesh generation program.

Parallel capabilities can be useful as well as tetrahedral meshes can become computationally expensive to generate. In acoustic and elastic wave propagation, a minimum number between 5 to 10 vertices per wavelength (~250-2000 m) of the source wavelet is applied to ensure the subsequent simulation is numerically accurate. For example, in a 3D mesh of a benchmark FWI model EAGE Salt [@doi:10.1190/1.1437283] this requires 507,862 cells when resolving a 2 Hz source frequency using 5 vertices per wavelength. The number of cells increases by a factor of approximately 8 (approximately 4,082,978 cells) as the source frequency is doubled.


# Core functionality

  1. The creation of 2D/3D graded mesh size functions defined on axis-aligned regular Cartesian grids. These mesh sizing functions encode mesh resolution distributions that conform to the variations from inputted seismic velocity model data and are distributed according to several heuristics [see @SeismicMeshDocs for further details]. Mesh size function grading is accomplished using [@persson2006mesh].

  2. Distributed memory parallelism. The generation of potentially large (> 10 million cells) high-geometric quality triangular or tetrahedral meshes in either using distributed memory parallelism with mesh resolution following sizing functions.

  3. An implementation of a 3D degenerate (i.e., sliver) tetrahedral element removal technique [@tournois2009perturbing] to bound a mesh quality metric while preserving the domain structure. Note that 2D mesh generation does not suffer from the formation of degenerate elements.

 Similar to other meshing programs such as `gmsh` [@doi:10.1002/nme.2579], `SeismicMesh` [@SeismicMeshDocs] enables both generation and improvement of simplical meshes through a Python scripting-based application programming interface.

The mesh's domain geometry is defined as the 0-level set of a signed distance function (SDF), which avoids the need to have explicit geometry information defining the boundary and can be particularly useful in geophysical domains.


# Performance Comparison

The 2D/3D serial performance against `gmsh` and `cgal` [@cgal:rty-m3-20b; @cgal:r-ctm2-20b] in terms of cell quality and creation time where cell quality is defined as dimension (2 or 3) multiplied by the circumcircle radius divided by the incircle radius. This cell quality is between 0 and 1, where 1 is a perfectly symmetrical simplex. For the analytical geometries (e.g. disk and ball) `gmsh` is the fastest to generate a mesh and then performance is approximately similar for both `SeismicMesh` and `cgal` with `cgal` outperforming `SeismicMesh` for the disk and vice versa for the ball. `gmsh` and `cgal` produce higher minimum cell qualities overall than `SeismicMesh` but importantly all generators are capable of producing sliver-free meshes. With that said, `SeismicMesh` produces higher mean cell qualities by about 5-10\% as compared to `gmsh` and `cgal`, which leads to reduced finite element matrix condition numbers. For the two seismic domains (e.g., BP2004 and EAGE), `SeismicMesh` is faster than `gmsh` for the 2D BP2004 benchmark but slightly slower for the 3D EAGE benchmark. Mesh quality results similarly follows the results observed in the two analytical cases. Interpolant-based mesh sizing functions significantly slow the mesh generation time of `gmsh` increasing its mesh generation time by a factor of 3x. The `DistMesh` algorithm used in `SeismicMesh` requires far less calls to the sizing function. For example in the EAGE benchmark, `gmsh` calls the sizing function 95,756 times whereas `DistMesh` calls it just 26 times.

![The mesh creation time (left columns) and resulting cell quality (right columns) for the four benchmark problems studied over a range of problem sizes. Two analytical problems (disk and a ball) and two non-analytical problems with sizing functions defined via regular gridded interpolants (BP2004 and EAGE). For the panels that show cell quality, solid lines indicate the mean and dashed lines indicate the minimum cell quality in the mesh. See the `SeismicMesh` github repository for information regarding the benchmarks and version of programs used. \label{fig:benchmark}](Benchmarks.png)


# Parallelism

A simplified version of the parallel Delaunay algorithm proposed by [ @peterka2014high] is implemented inside the `DistMesh` algorithm, which does not consider sophisticated domain decomposition or load balancing yet. \autoref{fig:speedup} shows a peak speed-up of approximately 4 times using 11 cores when performing 50 meshing iterations to generate both the light and heavy meshes of the EAGE P-wave velocity model. While the parallel performance is not perfect at this stage of development, the capability reduces the generation time of this relatively large example (e.g., 33 M cells) from 92.0 minutes to approximately 21.7 minutes. Results indicate that Python array manipulations dominate parallel execution time and inhibit scalability. The machine used for this experiment was an Intel Xeon Gold 6148 machine clocked at 2.4 GHz  with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.

![The speedup (left-panel) as compared to the serial version of SeismicMesh V3.1.0 for a relatively light and heavy mesh each adapted to P-wave data from the EAGE Salt seismic velocity model. The total mesh generation wall-clock time is annotated in decimal minutes next to each point. The panel on the right hand side shows the mesh generation rate normalized by the number of total number of cells in the mesh. \label{fig:speedup}](Performance.png)


# Ongoing and future applications

 Some future applications for this software:

 * `SeismicMesh` is being used by a group of researchers to build 2D/3D meshes for a seismological FEM model that has been developed in the Firedrake computing environment [@10.1145/2998441].

 * The usage of SDF to implicitly define the meshing domain presents potential use cases in a topology-optimization framework [@laurain2018level] for modeling the sharp interface of salt-bodies in seismological domains. In these applications, the 0-level set of a SDF is used to demarcate the boundary of the feature. Each inversion iteration, updates to an objective functional produce modifications to the 0-level set. In this framework, `SeismicMesh` can be embedded within the inversion algorithm to generate and adapt meshes.

 * Much like how the original `DistMesh` program has been used, `SeismicMesh` can be adapted for other domain-specific applications besides seismology (e.g., fluid dynamics, astrophysics, and oceanography). An open source project project is already under way to use the same mesh generation technology for a Python version of `OceanMesh2D` to build industrial-grade meshes of coastal oceans [@roberts2019oceanmesh2d].

 We expect future extensions of the program to introduce better domain decomposition algorithms to improve parallel performance, and support for both immersed and periodic mesh generation.

# Acknowledgements

This research was carried out in association with the ongoing R&D project registered as ANP 20714-2, "Software technologies for modelling and inversion, with applications in seismic imaging"  (University of São Paulo / Shell Brasil / ANP).

# References
