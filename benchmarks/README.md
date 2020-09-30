Performance
------------
Here we compare `SeismicMesh` against well-established existing mesh generation tools such as [cgal](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html). Specifically:

    * a comparison in mesh creation speed in terms of wall-clock time and throughput.
    * a comparison in cell quality.
    * a comparison of dihedral angle distributions.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available:

    benchmark_cuboid: # a simple box-type geometry with linearly varying mesh resolution along the x-axis.
    benchmark_sphere: # a sphere with a ring of higher resolution near the center.

Run `python benchmark_cuboid.py` to run all benchmarks for a particular domain (e.g., cuboid). Run `python benchmark_cuboid.py --method METHODNAME` to select either `cgal` using [pygalmesh](https://github.com/nschloe/pygalmesh) or `sm` to use `SeismicMesh`.

Results
---------------

Using [termplotlib](https://github.com/nschloe/termplotlib) and [meshplex](https://github.com/nschloe/meshplex) to calculate some mesh statistics, the benchmarks produce histograms of [dihedral angles](https://en.wikipedia.org/wiki/Dihedral_angle#:~:text=A%20dihedral%20angle%20is%20the,line%20as%20a%20common%20edge) in the cells and histograms of cell quality.


![The computer used for benchmarking is a PC running MacOS with Dual-Core Intel Core i5 clocked at 2.00 GHz with 8GB of RAM. Both mesh generation programs have been compiled with g++ v8.3.0 with the -O3 option. These benchmarks have been done using CGAL v5.0 and SeismicMesh v3.0.3](https://user-images.githubusercontent.com/18619644/94751279-89e18700-035e-11eb-9ddf-b42995e4a041.jpg)


Notes
-----
* Mesh generation with `cgal` is accomplished via [pygalmesh](https://github.com/nschloe/pygalmesh)
* For CGAL's mesh generator in 3D, all default quality options are assumed (facet angle bound of 30 degrees and the radius edge bound 2--to their theoretical limit). A `cell_size` function is passed to create variable mesh resolution in a way that is equivalent to the mesh size function in `SeismicMesh`.
* Mesh generation with `gmsh` is accomplished via [pygmsh](https://github.com/nschloe/pygmsh). WIP
* All programs here are executed in seqeuntial mode. It's important to note however that a significant speed-up can be achieved for moderate to large problems using the [parallel capabilities](https://seismicmesh.readthedocs.io/en/par3d/tutorial.html#basics) in `SeismicMesh`. Threading based parallelism can be employed with `gmsh` and `cgal` but these benchmarks have not been explored here.
