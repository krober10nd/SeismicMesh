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
 - name: Reearch Center for Gas Innovation, Escola Politécnica da Universidade de São Paulo, São Paulo, Brazil.
   index: 1
 - name:  Dept. of Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, 156 Fitzpatrick Hall, Notre Dame, IN, U.S.A.
   index: 2
date: 20 August 2020
bibliography: paper.bib

---
# Summary

`SeismicMesh` is a Python package for generating quality two and three dimensional triangular and tetrahedral meshes from seismic velocity models using signed distance functions to define geometry. C++ code is wrapped with Python for high performance without losing flexibility or ease-of-use. Mesh generation is accomplished in either serial or using distributed memory parallelism through modifications to a smoothing-based mesh generation algorithm known as `DistMesh` [@doi:10.1137/S0036144503429121; @peterka2014high]. The program uses the Computational Geometry Algorithms Library [@cgal:hs-chdt3-20a] to perform many geometric operations efficiently while avoiding numerical round-off problems


  # Background

Generating a high-geometric quality graded mesh for a geophysical domain represents a modern challenge for sophisticated seismological workflows using the Finite Element Method (FEM). In these applications, a domain is discretized typically with triangles/tetrahedral elements that vary widely in size around features of interest. In seismology, these meshes are commonly used with the FEM to solve partial differential equations that model acoustic or elastic waves. Geophysical exploration studies solve inverse problems with these wave propagators for seismic velocity model building applications such as Full Waveform Inversion (FWI) [@doi:10.1190/1.1441754; @virieux2009overview] and Reverse Time Migration (RTM) [@10.1093/gji/ggv380]. In these inversion problems, many forward and adjoint simulations are computed, which motivates a mesh that needs to efficiently discretize the domain and be accurate to simulate  with.

# Statement of Need

Self-contained automatic meshing workflows that incorporate geophysical data to generate unstructured meshes have been successful in several geophysical domains such as coastal ocean modeling [e.g., `OceanMesh2D` @roberts2019oceanmesh2d, `GMT` and `Terreno` @gorman2008systematic] and reservoir modeling [e.g., `MeshIT` @cacace2015meshit] to name a few. The overall motivation behind these tools is that, besides providing a tool to generate a high-quality mesh, there are numerous application-specific aspects to consider in the design of the mesh. For example, the distribution of elements must accurately represent physical features in the domain such as faults and the location of wells. At the same time, the distribution of elements must be placed in such a way  that the stability requirements of the numerical solver are satisfied. For example, mesh resolution must respect the Courant-Friedrichs-Lewey (CFL) condition to ensure stability. Material discontinuities in the interior of the domain need to be resolved to ensure reflection and  refraction of waves are modeled accurately.  These domain-specific aspects motivate the development of a self-contained mesh generation program.

 # Core functionality

 A schematic of `SeismicMesh` is shown in \autoref{fig:workflow}. The core functionality is as follows:

  1. The creation of 2D/3D graded mesh size functions defined on regular Cartesian grids. These mesh sizing functions contain mesh resolution distributions that conform to the variations from inputted seismic velocity model data and are distributed according to several heuristics [see @SeismicMeshDocs for further details]. Note that mesh size function grading is accomplished using [@persson2006mesh].

  2. The generation of potentially large (> 10 million cells) high-geometric quality triangular or tetrahedral meshes in either serial or using distributed memory parallelism.

  3. An implementation of a 3D degenerate (sliver) tetrahedral element removal technique [@tournois2009perturbing] to bound a mesh quality metric while preserving the domain structure. Note that 2D Delaunay mesh generation does not suffer from the formation of sliver elements.

  4. The ability to represent irregular domains implicitly using signed distance functions.

 ![A schematic of `SeismicMesh`. On the right hand side, a P-wave seismic velocity in the Canadian Rockies is shown [@gray1995migration]. Note that mesh linting is a sequence of mesh improvement operations and checks that occur after mesh generation has ceased. \label{fig:workflow}](Workflow.jpg)

 Similar to other meshing programs such as `gmsh` [@doi:10.1002/nme.2579], `tetgen` [@si2015tetgen] and `mmg` [@mmg], `SeismicMesh` [@SeismicMeshDocs] provides both generation and improvement of meshes through a scripting-based application programming interface. A point of difference from the aforementioned software is a convenience class that can be used to generate graded mesh sizing functions from seismic velocity models. With this, mesh resolution is distributed to resolve material variations within the interior of the domain.

Domain geometry is defined as the 0-level set of a signed distance function (SDF), which avoids the need to have explicit geometry information defining the boundary. Geometries such as the free surface  boundary, seafloor, volcanoes, and salt-bodies are characterized by pronounced changes to seismic velocities making it possible to demarcate these regions by thresholding seismic velocity. A capability is provided in to create signed distance functions from isocontours of seismic velocity using the Fast-Marching method [@sethian1996fast].

 # Parallelism

Time-domain FWI and RTM require relatively high source frequencies (e.g., 5-7 Hz) to produce high-resolution seismic velocity images, which lead to large meshes. A minimum number of 5 to 10 vertices per wavelength (~250-2000 m) of the source wavelet is applied to ensure the subsequent simulation is numerically accurate, which makes the generation of tetrahedral meshes computationally expensive. For example, a 3D mesh of a benchmark FWI model EAGE Salt [@doi:10.1190/1.1437283] requires 507,862 cells when resolving a 2 Hz source frequency using 5 vertices per wavelength. The umber of cells by a factor of approximately 8 (approximately 4,082,978 cells) as the source frequency is doubled.

\autoref{fig:speedup} shows a peak speed-up of approximately 6.60 times using 11 cores when performing 50 meshing iterations to generate an approximately 4 million cell mesh. The usage of 11 cores reduces the generation time of this example from 20 minutes to approximately 2 wall-clock minutes.  The machine used was 2 Intel Xeon Gold 6148 clocked at 2.4 GHz  (40 cores in total, 27 MB cache, 10.4 GT/s) with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.


 ![Speed-up (left-axis) as compared to the sequential version of the program and wall-clock time in minutes to generate a 3D mesh (approximately 4.6 M cells) for the EAGE Salt seismic velocity model. The panel on the right hand side shows the a slice through the center of the generated mesh. \label{fig:speedup}](Performance.jpg)


 # Future applications

 Here are some future applications for this software:

 * `SeismicMesh` is being used by a group of researchers to build 2D/3D meshes for a seismological FEM model that has been developed in the Firedrake computing environment [@10.1145/2998441].

 * The usage of SDF to implicitly define the meshing domain presents potential use cases in a topology-optimization framework [@laurain2018level] for modeling the sharp interface of salt-bodies in seismological domains. In these applications, the 0-level set of a SDF is used to demarcate the boundary of the feature. Each inversion iteration, updates to an   objective functional produce modifications to the 0-level set. In this framework, `SeismicMesh` can be embedded within the inversion algorithm to generate and adapt meshes so that they conform accurately to the 0-level set.

 * Much like how the original `DistMesh` program has been used, `SeismicMesh` can be adapted for other domain-specific applications besides seismology (e.g., fluid dynamics, astrophysics, and oceanography). An open source project project is already under way to use the same mesh generation technology for a Python version of `OceanMesh2D` to build industrial-grade meshes of coastal oceans [@roberts2019oceanmesh2d].

# Acknowledgements

This research was carried out in association with the ongoing R&D project registered as ANP 20714-2, "Software technologies for modelling and inversion, with applications in seismic imaging"  (University of São Paulo / Shell Brasil / ANP).

# References
