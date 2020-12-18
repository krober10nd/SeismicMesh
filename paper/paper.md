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

SeismicMesh is a Python package for simplex mesh generation in two or three dimensions. As an implementation of DistMesh [@doi:10.1137/S0036144503429121], it produces high-quality meshes at the expense of speed. For increased efficiency, the core package is written in C++, works in parallel, and uses the Computational Geometry Algorithms Library [@cgal:hs-chdt3-20b]. SeismicMesh can also produce mesh-density functions from seismological data to be used in the mesh generator.

# Background

Generating a high-quality graded mesh for a geophysical domain represents a challenge for seismological modeling using the finite element method (FEM). In these applications, a domain is discretized typically with triangular/tetrahedral elements that vary widely in size around features of interest. These meshes are commonly used with the FEM to solve partial differential equations that model acoustic or elastic waves, which are used in seismic velocity model building algorithms such as full waveform inversion (FWI) [@doi:10.1190/1.1441754; @virieux2009overview] and reverse time migration [@10.1093/gji/ggv380].

# Statement of Need

Despite the fact that many mesh generation programs exist such as Gmsh [@doi:10.1002/nme.2579] and CGAL [@cgal:rty-m3-20b; @cgal:r-ctm2-20b], it is uncommon to find capabilities that incorporate geophysical data into the mesh generation process to appropriately size elements. This in part contributes to the reality that automatic mesh generation for geophysical domains is not user-friendly.

Some packages have been created to script mesh generation from geophysical datasets such as in coastal ocean modeling [@roberts2019oceanmesh2d; @gorman2008systematic] and reservoir modeling [@cacace2015meshit]. In a similar manner, the aim of this package is to provide a straightforward Python package to script mesh generation directly from seismic velocity models. This is accomplished first by building a mesh density function using seismic velocity data and then supplying these inputs to a mesh generator that can use these inputs and operate at scale.

The mesh density function can be used as input other mesh generators. However, the usage of a sizing function can have significant impact on the mesh generation performance. For example, Gmsh’s advancing front and Delaunay refinement methods construct the mesh incrementally and do not permit vectorization, which leads to reduced performance at scale in 2D/3D. In contrast, the DistMesh algorithm takes advantage of vectorization when querying a complex mesh density function making it efficient and competitive to Gmsh for this kind of meshing problem.

# Core functionality

  1. The creation of 2D/3D graded mesh size functions defined on axis-aligned regular Cartesian grids. These mesh sizing functions encode mesh resolution distributions that conform to the variations from inputted seismic velocity model data and are distributed according to several heuristics [see @SeismicMeshDocs for further details]. Mesh size function grading is accomplished using [@persson2006mesh].

  2. Distributed memory parallelism. The generation of potentially large (> 10 million cells) high-quality triangular or tetrahedral meshes in either using distributed memory parallelism with mesh resolution following sizing functions.

  3. An implementation of a 3D degenerate (i.e., sliver) tetrahedral element removal technique [@tournois2009perturbing] to bound a mesh quality metric. Note that 2D mesh generation does not suffer from the formation of degenerate elements.

Similar to other meshing programs such as Gmsh, SeismicMesh [@SeismicMeshDocs] enables generation of simplex meshes through a Python application programming interface.

The mesh's domain geometry is defined as the 0-level set of a signed distance function (SDF), which avoids the need to have explicit geometry information defining the boundary and can be particularly useful in geophysical domains.

# Performance Comparison

We compare the 2D/3D serial performance in terms of cell quality and mesh creation time between SeismicMesh, Gmsh [@doi:10.1002/nme.2579] and CGAL [@cgal:rty-m3-20b; @cgal:r-ctm2-20b]. The cell quality is defined as the product of the topological dimension of the mesh (2 or 3) and the incircle radius divided by the circumcircle radius and ranges between 0 and 1, where 1 is a perfectly symmetrical simplex. In mesh generation, there is always a trade-off between generation speed and mesh quality. We find that Gmsh produces high-quality meshes by far the fastest, SeismicMesh will produce meshes with the best quality, but much slower. Gmsh becomes comparatively slow when a user-defined mesh-density function is involved, which is SeismicMesh's primary use case.

For the two seismic domains (e.g., BP2004 and EAGE), SeismicMesh is faster than Gmsh for the 2D BP2004 benchmark but slightly slower for the 3D EAGE benchmark at scale. CGAL is not competitive for the 3D benchmark and is therefore not shown. Interpolant-based mesh sizing functions significantly slow the mesh generation time of Gmsh by a factor of $\sim 3$ as Gmsh calls the sizing function for each point individually (e.g., 95,756 times) whereas SeismicMesh does it for all points at once each meshing iteration (e.g., 26 times). 

![The mesh creation time (left columns) and resulting cell quality (right columns) for the four benchmarks studied over a range of problem sizes. Two analytical problems (disk and a ball) and two non-analytical problems with sizing functions defined via regular gridded interpolants (BP2004 and EAGE). For the panels that show cell quality, solid lines indicate the mean and dashed lines indicate the minimum cell quality in the mesh. See the SeismicMesh github repository for information regarding the benchmarks and version of programs used. \label{fig:benchmark}](Performance.svg)

# Parallelism

A simplified version of the parallel Delaunay algorithm proposed by [ @peterka2014high] is implemented inside the DistMesh algorithm, which does not consider sophisticated domain decomposition or load balancing yet. \autoref{fig:speedup} shows a peak speed-up of approximately 6 times using 11 cores when performing 50 meshing iterations to generate the 33M cell mesh of the EAGE P-wave velocity model. While the parallel performance is not perfect at this stage of development, the capability reduces the generation time of this relatively large example (e.g., 33 M cells) from 90.0 minutes to approximately 15.6 minutes. Results indicate that the simple domain decomposition approach inhibit perfect scalability. The machine used for this experiment was an Intel Xeon Gold 6148 machine clocked at 2.4 GHz  with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.

![The speedup (left-panel) as compared to the serial version of SeismicMesh V3.2.0 for a relatively light and heavy mesh each adapted to P-wave data from the EAGE Salt seismic velocity model. The total mesh generation wall-clock time is annotated in decimal minutes next to each point. The panel on the right hand side shows the mesh generation rate normalized by the number of total number of cells in the mesh. \label{fig:speedup}](Benchmarks.svg)

# Ongoing and future applications

 Some future applications for this software:

 * SeismicMesh is being used by a group of researchers to build 2D/3D meshes for a seismological FEM model that has been developed in the Firedrake computing environment [@10.1145/2998441].

 * The usage of SDF to implicitly define the meshing domain presents potential use cases in a topology-optimization framework [@laurain2018level] for modeling the sharp interface of salt-bodies in seismological domains. In these applications, the 0-level set of a SDF is used to demarcate the boundary of the feature. Each inversion iteration, an optimization problem is solved to produce modifications to the location of the 0-level set. In this framework, SeismicMesh can be used within the inversion algorithm to generate and adapt meshes.

 * Much like how the original DistMesh program has been used, SeismicMesh can be adapted for other domain-specific applications besides seismology (e.g., fluid dynamics, astrophysics, and oceanography). An open source project project is already under way to use the same mesh generation technology for a Python version of OceanMesh2D to build industrial-grade meshes of coastal oceans [@roberts2019oceanmesh2d].

We expect future extensions of the program to introduce better domain decomposition algorithms to improve parallel performance.

# Acknowledgements

This research was carried out in association with the ongoing R&D project registered as ANP 20714-2, "Software technologies for modelling and inversion, with applications in seismic imaging" (University of São Paulo / Shell Brasil / ANP). We would like to thank the two reviewers who helped improved the presentation and quality of this manuscript and package greatly.

# References
