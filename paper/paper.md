---
title: 'SeismicMesh: Triangular meshing for seismology'
tags:
  - Python
  - Seismology
  - Mesh generation
  - Full Waveform Inversion
  - Parallel mesh generation
authors:
  - name: Keith Jared Roberts
    affiliation: 1
affiliations:
 - name: Escola Politécnica da Universidade de São Paulo
   index: 1
date: 4 August 2020
bibliography: paper.bib

---

# Summary

`SeismicMesh` is a Python package for generating quality two and three dimensional triangular and tetrahedral
meshes directly from seismic velocity models using signed distance functions to define geometry. Low-level C++ code is wrapped with Python for high performance at scale without losing flexibility or ease-of-use. At a low level, the program makes direct calls to the Computational Geometry Algorithms Library [@cgal:hs-chdt3-20a] to robustly and efficiently produce triangulations. Mesh generation is accomplished through modifications to a smoothing-based mesh generation algorithm known as `DistMesh` [@doi:10.1137/S0036144503429121]. Computationally expensive operations are parallelized through an implementation of a parallel Delaunay triangulation algorithm [@peterka2014high].

# Background

Generating a high-quality graded mesh for a geophysical domain represents a modern
challenge for sophisticated geophysical simulation workflows. In these applications,
a geophysical domain is discretized typically with triangles/tetrahedral elements that vary widely in size around features of interest. In geophysics, these meshes are commonly used with Finite Element
Method (FEM) to solve Partial Differential Equation that model the acoustic or elastic wave equation. Geophysical exploration studies use these meshes as a discretization to solve inverse problems such as Full Waveform Inversion (FWI) [@doi:10.1190/1.1441754; @virieux2009overview]
and Reverse Time Migration (RTM) [@10.1093/gji/ggv380]. In these inverse problems, potentially hundreds of forward and adjoint simulations are conducted, which motivates a mesh that should be both efficient and accurate to simulate with.

# Rationale

Self-contained automatic modeling workflows that can incorporate geophysical datasets to produce unstructured meshes used in numerical simulations have shown success in several geophysical domains such as ocean modeling [@roberts2019oceanmesh2d] and reservoir modeling `MeshIT` [@cacace2015meshit]. Besides having an implementation to generate the triangulation, there are numerous application-specific aspects to consider for a geophysical inverse problem. The distribution of elements must accurately represent physical features in the domain such as faults and the location of wells. At the same time, mesh resolution must respect the Courant-Friedrichs-Lewey condition to ensure numerical stability. In applications such as FWI and RTM, material discontinuities in the interior of the domain need to be resolved to ensure reflection and refraction of waves are modeled accurately as compared to observed data. Domain extensions are also needed to built into the mesh to absorb outgoing waves and minimize artificial boundary reflections. In domains with irregular free surface boundaries, explicit geometry data of the boundary surface may not exist requiring the use of external programs to create it and then to mesh it.

# Software architecture

A high-level depiction of the execution of 'SeismicMesh' is shown in Figure \autoref{fig:workflow}. The core functionality is as follows:

 1. The creation of 2D/3D graded mesh size functions defined on regular Cartesian grids with mesh resolution variations that conform to the variations from inputted seismic velocity model data.

 2. The generation of potentially large (> :math`10` million cells) high-geometric quality triangular or tetrahedral meshes in either serial or using distributed parallelism that follow a mesh size function.

 3. An implementation of a 3D degenerate tetrahedral element removal technique [@tournois2009perturbing] to bound the minimum mesh quality while preserving the domain structure.


![A workflow to generate a mesh using `SeismicMesh`. On the right hand side, a BP2004 dataset of the P-wave seismic velocity in the Canadian Rockies [@gray1995migration]  \label{fig:workflow}](Workflow.jpg)

Similar to other meshing programs such as `gmsh` [@doi:10.1002/nme.2579], `tetgen` [@si2015tetgen] and `mmg` [@mmg], `SeismicMesh` provides both generation and improvement of meshes through a scripting-based approach. However, one point of difference from the aforementioned software programs is a convenience class that can be used to generate graded mesh sizing functions directly from geophysical datasets such as a seismic velocity model. Using this capability, mesh resolution is distributed to resolve material variations within the interior of the domain while at the same time statisfying numerical stability requirements and requiring little effort on part of the user.

In `SeismicMesh`, the domain geometry is defined as the 0-level set of a signed distance function (SDF). To create irregular geometries, set operations (e.g., union and intersection) can be performed with several SDFs. The usage of SDFs avoids the need to have explicit geometry information defining the boundary of the domain while supporting the meshing of complex and irregularly shaped domains. Geometries such as the free surface boundary, seafloor, and salt-bodies are characterized by pronounced changes to seismic velocities making it possible to demarcate these regions using seismic velocity ranges. To mesh these features, a capability is provided in `SesimicMesh` to create signed distance functions from isocontours of seismic velocity using the Fast-Marching method [@sethian1996fast].

# Parallelism

In applications such as time-domain FWI and RTM, relatively high source frequencies are required (e.g., 5-7 Hz) to produce high-resolution seismic velocity images. Depending on the spatial polynomial order of the finite elements, these models require a minimum number of vertices per wavelength (e.g., 5 to 10) of the source wavelet to ensure the subsequent numerical simulation can be considered numerically accurate. However, this can make the generation of tetrahedral meshes prohibitively computationally expensive for a realistic 3D inversion domain, which motivates the implemented distributed memory paralllelism. For example, a 3D mesh of a benchmark FWI model EAGE Salt [@doi:10.1190/1.1437283] requires 507,862 cells when resolving a 2 Hz source frequency using 5 vertices per wavelength. For the same model discretized for source wavelet with a peak frequency of 4 Hz, the tetrahedral mesh increases in number of cells by a factor of approximately 8 and becomes 4,082,978 cells.

The implemented distributed memory parallelism makes generating high-quality meshes for high-frequency FWI and RTM applications feasible on the order of minutes. Figure \autoref{fig:speedup} shows a peak speed-up of approximately 7 times using 11 cores when performing 50 meshing iterations with `SeismicMesh` to generate a 4 million cell mesh. The machine used was 2 Intel Xeon Gold 6148 clocked at 2.4 GHz (40 cores in total, 27 MB cache, 10.4 GT/s) with 192 GB of RAM connected together with a 100 Gb/s InfiniBand network.


![Speed-up obtained when generate a 3D mesh (approximately 4.6 M cells) for the EAGE Salt seismic velocity model as compared to the sequential version of the program. The panel on the right hand side shows the a slice through the center of the mesh. \label{fig:speedup}](Performance.jpg)


# Future applications

The usage of SDF to define the meshing domain present potential use cases in the topology-optimization framework [@laurain2018level] for modeling the sharp interface of salt-bodies in seismological domains. In these applications, the 0-level set of the SDF is used to demarcate the boundary of the feature. Each inversion iteration, updates to an objective functional produce a new 0-level set. In this scenario, `SeismicMesh` can be embedded within the algorithm to generate/adapt meshes that conform to the updated 0-level set.


# Acknowledgements

Shell and Research Center for Gas Innovation

# References
