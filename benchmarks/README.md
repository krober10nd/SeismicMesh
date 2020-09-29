Benchmarking
------------
Here we compare *SeismicMesh* against well-established existing mesh generation approaches such as [CGAL](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html). Specifically:

    * a comparison in mesh creation speed in terms of wall-clock time.
    * a comparison in cell quality.
    * a comparison of dihedral angle distributions.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available in `benchmark.py`:

    benchmark_cuboid: # a simple box-type geometry with linearly varying mesh resolution
    benchmark_sphere: # a sphere with a ring of higher resolution in the center

Run `python benchmark_cuboid.py --method METHODNAME` with `cgal` using [pygalmesh](https://github.com/nschloe/pygalmesh) to mesh and `sm` using SeismicMesh to mesh with.

Expected output
---------------
Using  [termplotlib](https://github.com/nschloe/termplotlib), the benchmarks produce the following output:

```
┌───────────────────────┬──────────────────────────────────────────┬────────────────────────────────────┬──────────────────────┐
│                       │                                          │                                    │                      │
│  SeismicMesh          │  Mesh creation time (seconds):     8.27  │  Number of vertices:        11510  │                      │
│                       │                                          │  Number of cells:        62017     │                      │
│                       │                                          │                                    │                      │
├───────────────────────┼──────────────────────────────────────────┼────────────────────────────────────┼──────────────────────┤
│                       │                                          │                                    │                      │
│  dihedral angles      │  min angle:      10.584                  │  mesh quality metric               │  min quality: 0.137  │
│                █▇     │  avg angle:      42.594                  │                            ▇█▅     │  avg quality: 0.801  │
│               ███▃    │  max angle:      56.723                  │                           ████     │  max quality: 0.999  │
│              ▄████    │  std dev angle:   8.393                  │                          █████▇    │                      │
│              █████    │                                          │                         ▆██████    │                      │
│             ▇█████    │                                          │                        ▄███████    │                      │
│            ▂██████▂   │                                          │                       ▄████████▇   │                      │
│           ▂████████   │                                          │                      ▄██████████   │                      │
│          ▁█████████   │                                          │                    ▃▆███████████   │                      │
│        ▂▅██████████   │                                          │                ▁▃▄▆█████████████▃  │                      │
│  ▁▃▄▅▆█████████████▄  │                                          │  ▁▁▂▂▃▃▃▃▄▅▅▅▆▇██████████████████  │                      │
│                       │                                          │                                    │                      │
└───────────────────────┴──────────────────────────────────────────┴────────────────────────────────────┴──────────────────────┘
```
