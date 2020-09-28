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

Expected output
---------------
```
$ python benchmarks/benchmark_uniform_cuboid.py --method sm
Commencing mesh generation with 9261 vertices on rank 0.
Iteration #1, max movement is 0.015040, there are 9261 vertices and 48000 cells
     Elapsed wall-clock time 0.233844 :
Iteration #11, max movement is 0.067096, there are 9261 vertices and 54648 cells
     Elapsed wall-clock time 0.226372 :
Iteration #21, max movement is 0.029739, there are 9261 vertices and 52903 cells
Termination reached...maximum number of iterations reached.
Will attempt 30 to bound the dihedral angles...
Enforcing a min. dihedral bound of: 10.0 degrees...
Enforcing a max. dihedral bound of: 170.0 degrees...
Commencing sliver removal with 9260 vertices on rank 0.
On rank: 0 There are 639 slivers...
Termination reached...no slivers detected!
8.290160179138184
Elapsed time is: 8.290189027786255
┌───────────────────────┬──────────────────────────┬────────────────────────────────────┬──────────────────────┐
│                       │                          │                                    │                      │
│  dihedral angles      │  min angle:      10.577  │  mesh quality metric               │  min quality: 0.179  │
│               ▃█▄     │  avg angle:      60.000  │                           ▄█▁      │  avg quality: 0.777  │
│              ▄███     │  max angle:      56.520  │                          ████      │  max quality: 0.996  │
│              ████▃    │  std dev angle:   8.940  │                         ▆█████     │                      │
│             ██████    │                          │                        ▆██████▁    │                      │
│            ▆██████    │                          │                       ▄████████    │                      │
│            ███████    │                          │                      ▄█████████    │                      │
│           ████████▅   │                          │                     ▅██████████▃   │                      │
│         ▃██████████   │                          │                   ▂▇████████████   │                      │
│       ▃▆███████████   │                          │               ▂▃▅███████████████   │                      │
│  ▃▅▆▇██████████████▃  │                          │  ▂▂▃▄▄▄▅▅▅▆▆▇███████████████████▇  │                      │
│                       │                          │                                    │                      │
└───────────────────────┴──────────────────────────┴────────────────────────────────────┴──────────────────────┘
```
