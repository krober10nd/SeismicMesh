Benchmarking
------------
Here we compare *SeismicMesh* against well-established existing mesh generation approaches such as [CGAL](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html).
    * a comparison in mesh creation speed in terms of wall-clock time.
    * a comparison in cell quality.
    * a comparison of dihedral angle distributions.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available in `benchmark.py`:

    benchmark_uniform_cuboid: 1 # a simple box-type geometry with uniform mesh resolution

Run `python benchmark_uniform_cuboid.py --method METHODNAME` with `cgal` using [pygalmesh](https://github.com/nschloe/pygalmesh) to mesh and `sm` using SeismicMesh to mesh with.

