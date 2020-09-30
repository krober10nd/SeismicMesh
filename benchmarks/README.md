Performance
------------
Here we compare `SeismicMesh` against well-established existing mesh generation approaches such as [CGAL](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html). Specifically:

    * a comparison in mesh creation speed in terms of wall-clock time.
    * a comparison in cell quality.
    * a comparison of dihedral angle distributions.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available in `benchmark.py`:

    benchmark_cuboid: # a simple box-type geometry with linearly varying mesh resolution
    benchmark_sphere: # a sphere with a ring of higher resolution in the center

Run `python benchmark_cuboid.py --method METHODNAME` with `cgal` using [pygalmesh](https://github.com/nschloe/pygalmesh) to mesh and `sm` using `SeismicMesh` to mesh with.

Expected output
---------------

Using [termplotlib](https://github.com/nschloe/termplotlib) and [meshplex](https://github.com/nschloe/meshplex) to calculate some mesh statistics, the benchmarks produce histograms of [dihedral angles](https://en.wikipedia.org/wiki/Dihedral_angle#:~:text=A%20dihedral%20angle%20is%20the,line%20as%20a%20common%20edge) in the cells and the cell quality.


For example running `benchmark_sphere.py` produces the following output:

┌───────────────────────────┬────────────────────────────────────────────────────┬────────────────────────────────────────────┬──────────────────────┐
│                           │                                                    │                                            │                      │
│  CGAL                     │  Mesh creation time (seconds):    15.84            │  Number of vertices:       108407          │                      │
│                           │  Mesh creation speed (vertices/seconds):  6843.10  │  Number of cells:           18118          │                      │
│                           │                                                    │                                            │                      │
├───────────────────────────┼────────────────────────────────────────────────────┼────────────────────────────────────────────┼──────────────────────┤
│                           │                                                    │                                            │                      │
│  dihedral angles          │  min angle:      12.764                            │  mesh quality metric                       │  min quality: 0.189  │
│                    ▇█     │  avg angle:      44.003                            │                                  ▁▆█▁      │  avg quality: 0.797  │
│                   ▇██     │  max angle:      56.870                            │                                 ▁████      │  max quality: 0.999  │
│                   ████    │  std dev angle:   6.886                            │                                 █████▆     │                      │
│                  █████    │                                                    │                                ▇██████     │                      │
│                  █████    │                                                    │                               ▅████████    │                      │
│                 ██████    │                                                    │                              ▅█████████    │                      │
│                ▃██████▇   │                                                    │                             ▅██████████▅   │                      │
│               ▃▉███████   │                                                    │                           ▃█████████████   │                      │
│             ▁▄█▉███████   │                                                    │                        ▁▄▇██████████████▁  │                      │
│        ▁▁▂▄▆███▉███████▃  │                                                    │              ▁▁▁▂▂▃▄▅▆███▉███████████████  │                      │
│                           │                                                    │                                            │                      │
└───────────────────────────┴────────────────────────────────────────────────────┴────────────────────────────────────────────┴──────────────────────┘
┌───────────────────────────┬────────────────────────────────────────────────────┬────────────────────────────────────────────┬──────────────────────┐
│                           │                                                    │                                            │                      │
│  SeismicMesh              │  Mesh creation time (seconds):    21.50            │  Number of vertices:       155364          │                      │
│                           │  Mesh creation speed (vertices/seconds):  7226.69  │  Number of cells:           25774          │                      │
│                           │                                                    │                                            │                      │
├───────────────────────────┼────────────────────────────────────────────────────┼────────────────────────────────────────────┼──────────────────────┤
│                           │                                                    │                                            │                      │
│  dihedral angles          │  min angle:      10.555                            │  mesh quality metric                       │  min quality: 0.091  │
│                     ▅█    │  avg angle:      46.413                            │                                       █▂   │  avg quality: 0.859  │
│                     ██    │  max angle:      56.935                            │                                      ▇██   │  max quality: 0.999  │
│                    ▅██▁   │  std dev angle:   7.112                            │                                     ▄███   │                      │
│                    ████   │                                                    │                                    ▁████   │                      │
│                   ▃████   │                                                    │                                    █████   │                      │
│                   █████   │                                                    │                                   ██████   │                      │
│                  ▄█████   │                                                    │                                  ▇██████▇  │                      │
│                 ▂██████   │                                                    │                                ▂▇████████  │                      │
│                ▃███████▁  │                                                    │                              ▂▅██████████  │                      │
│      ▁▁▁▂▂▂▃▄▅█▉████████  │                                                    │           ▁▁▁▁▁▁▁▁▂▂▂▂▃▃▄▄▅▆█████████████  │                      │
│                           │                                                    │                                            │                      │
└───────────────────────────┴────────────────────────────────────────────────────┴────────────────────────────────────────────┴──────────────────────┘

and for `benchmark_cuboid.py`....

┌───────────────────────────┬────────────────────────────────────────────────────┬────────────────────────────────────────────┬──────────────────────┐
│                           │                                                    │                                            │                      │
│  CGAL                     │  Mesh creation time (seconds):    11.47            │  Number of vertices:        32651          │                      │
│                           │  Mesh creation speed (vertices/seconds):  2846.85  │  Number of cells:          180996          │                      │
│                           │                                                    │                                            │                      │
├───────────────────────────┼────────────────────────────────────────────────────┼────────────────────────────────────────────┼──────────────────────┤
│                           │                                                    │                                            │                      │
│  dihedral angles          │  min angle:      12.778                            │  mesh quality metric                       │  min quality: 0.198  │
│                    ▅█     │  avg angle:      44.652                            │                                   ▃█▅      │  avg quality: 0.813  │
│                    ██▂    │  max angle:      56.977                            │                                  ▂███▃     │  max quality: 0.999  │
│                   ████    │  std dev angle:   6.610                            │                                  █████     │                      │
│                   ████    │                                                    │                                 ▇█████▂    │                      │
│                  █████    │                                                    │                                ▄███████    │                      │
│                 ▁█████▁   │                                                    │                               ▃████████    │                      │
│                 ███████   │                                                    │                              ▄█████████▅   │                      │
│                ▆███████   │                                                    │                            ▁▅███████████   │                      │
│              ▂▆▉███████   │                                                    │                         ▁▂▅█████████████   │                      │
│         ▁▂▃▅▇██▉███████▃  │                                                    │               ▁▁▁▂▂▃▃▄▆▆█▉███████████████  │                      │
│                           │                                                    │                                            │                      │
└───────────────────────────┴────────────────────────────────────────────────────┴────────────────────────────────────────────┴──────────────────────┘
┌───────────────────────────┬────────────────────────────────────────────────────┬────────────────────────────────────────────┬──────────────────────┐
│                           │                                                    │                                            │                      │
│  SeismicMesh              │  Mesh creation time (seconds):    28.50            │  Number of vertices:        45485          │                      │
│                           │  Mesh creation speed (vertices/seconds):  1596.01  │  Number of cells:          260120          │                      │
│                           │                                                    │                                            │                      │
├───────────────────────────┼────────────────────────────────────────────────────┼────────────────────────────────────────────┼──────────────────────┤
│                           │                                                    │                                            │                      │
│  dihedral angles          │  min angle:      10.556                            │  mesh quality metric                       │  min quality: 0.149  │
│                     █▅    │  avg angle:      45.937                            │                                      ▂█    │  avg quality: 0.859  │
│                    ▃██    │  max angle:      56.926                            │                                      ███   │  max quality: 1.000  │
│                    ███    │  std dev angle:   7.182                            │                                     ▅███   │                      │
│                   ▁███▂   │                                                    │                                    ▂████   │                      │
│                   █████   │                                                    │                                    █████   │                      │
│                   █████   │                                                    │                                   ▆█████   │                      │
│                  ██████   │                                                    │                                  ▅██████▄  │                      │
│                 ▄██████   │                                                    │                                 ▅████████  │                      │
│                ▄███████   │                                                    │                              ▁▄██████████  │                      │
│      ▁▁▁▂▂▃▃▄▆█▉████████  │                                                    │          ▁▁▁▁▁▁▁▁▂▂▂▂▂▃▃▄▄▅▆▇████████████  │                      │
│                           │                                                    │                                            │                      │
└───────────────────────────┴────────────────────────────────────────────────────┴────────────────────────────────────────────┴──────────────────────┘

![Benchmark meshes: cuboid and sphere](https://user-images.githubusercontent.com/18619644/94724275-69023d00-0330-11eb-8b47-d6ede11f46cd.jpg)

Notes
-----
* For CGAL's mesh generator in 3D, all default quality options are assumed. A `cell_size` function is passed to create variable resolution.
* SeismicMesh is run here in serial mode. It's important to note however that a significant speed-up can be achieved for moderate to large problems using its [parallel capabilities](https://seismicmesh.readthedocs.io/en/par3d/tutorial.html#basics).
