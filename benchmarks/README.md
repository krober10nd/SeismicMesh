Packages for benchmarks
------------------------
Some extra packages are required for benchmarking. These can be installed through the following command:
```
pip install SeismicMesh[benchmarking]
```

Performance
------------
Here we compare SeismicMesh against well-established existing mesh generation tools such as [cgal](https://doc.cgal.org/latest/Mesh_3/) and [gmsh](https://gmsh.info/doc/texinfo/gmsh.html). Specifically:

    * a comparison in mesh creation speed in terms of wall-clock time and throughput.
    * a comparison in cell quality.

where cell quality is defined as d * circumcircle_radius / incircle_radius (where d is 2 for triangles and 3 for tetrahedra). The value is between 0 and 1, where 1 is a perfectly symmetrical simplex.

Benchmarks
----------

The following set of benchmark problems are available:

    benchmark_BP2004: # meshing the 2D BP2004 seismic velocity model
    benchmark_EAGE: # meshing the 3D EAGE Salt velocity model
    benchmark_disk: # a uniform 2D mesh of a unit disk
    benchmark_ball: # a unit ball with a ring of higher resolution near the center.

Run `python benchmark_sphere.py` to run all benchmarks for a particular domain (e.g., sphere). Run `python benchmark_sphere.py --method METHODNAME` to select either CGAL using [pygalmesh](https://github.com/nschloe/pygalmesh), Gmsh using [pygmsh](https://github.com/nschloe/pygmsh) or `sm` to use SeismicMesh.

* Note: 3D CGAL results are not shown for the EAGE benchmark because they are several times slower than Gmsh and SeismicMesh. 2D CGAL results are not shown for BP2004 since they do not support user-defined variable mesh density functions.

Results
---------------

The computer used for benchmarking is a PC running MacOS with Dual-Core Intel Core i5 clocked at 2.00 GHz with 8GB of RAM. All mesh generation programs have been compiled similarly with gcc v8.3.0 with the -O3 option. These benchmarks have been done using CGAL v5.0, gmsh 4.7.0, and SeismicMesh v3.1.4. Each statistic is reported as the average of 5 executions.

Using [termplotlib](https://github.com/nschloe/termplotlib) and [meshplex](https://github.com/nschloe/meshplex) to calculate some mesh statistics, the benchmarks produce histograms of each cells' minimum [dihedral angles](https://en.wikipedia.org/wiki/Dihedral_angle) and histograms of cell quality (which was described above).

**NOTE: 2D mesh sizing functions are not supported by CGAL**

Average speed statistics can be computed via [pytest-benchmark](https://pypi.org/project/pytest-benchmark/) which is set up to run each domain 5 times. For example:

```python
py.test --benchmark-max-time=360 benchmarks/benchmark_ball.py
```

produces for example:

```
=========================================================================================== test session starts ============================================================================================
platform darwin -- Python 3.7.4, pytest-6.0.2, py-1.9.0, pluggy-0.13.1
benchmark: 3.2.3 (defaults: timer=time.perf_counter disable_gc=False min_rounds=5 min_time=0.000005 max_time=360 calibration_precision=10 warmup=False warmup_iterations=100000)
rootdir: /Users/Keith/junk/SeismicMesh, configfile: pytest.ini
plugins: xdist-2.1.0, cov-2.10.1, benchmark-3.2.3, forked-1.3.0
collected 3 items

benchmarks/benchmark_ball.py ...                                                                                                                                                                   [100%]

-------------------------------------------------------------------------------- benchmark: 3 tests --------------------------------------------------------------------------------
Name (time in s)          Min                Max               Mean            StdDev             Median               IQR            Outliers     OPS            Rounds  Iterations
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_gmsh              4.1087 (1.0)       4.6007 (1.0)       4.3494 (1.0)      0.1953 (1.0)       4.2846 (1.0)      0.2935 (1.0)           2;0  0.2299 (1.0)           5           1
test_seismic_mesh     14.3477 (3.49)     16.1393 (3.51)     15.2947 (3.52)     0.7355 (3.77)     15.3189 (3.58)     1.2345 (4.21)          2;0  0.0654 (0.28)          5           1
test_cgal             16.2922 (3.97)     19.2086 (4.18)     17.4405 (4.01)     1.2892 (6.60)     16.9541 (3.96)     2.2134 (7.54)          1;0  0.0573 (0.25)          5           1
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Legend:
  Outliers: 1 Standard Deviation from Mean; 1.5 IQR (InterQuartile Range) from 1st Quartile and 3rd Quartile.
  OPS: Operations Per Second, computed as 1 / Mean
================================================================================ 3 passed, 20 warnings in 186.87s (0:03:06) ===============================================================================
```


Details on experiments
-----------------------
* Mesh generation with CGAL is accomplished via [pygalmesh](https://github.com/nschloe/pygalmesh) version 0.8.2
* For CGAL's 3D mesh generator, all default quality options are assumed (e.g., facet angle bound of 30 degrees and the radius edge bound 2--to their theoretical limit). A `cell_size` function is passed to create variable mesh resolution in a way that is approximately equivalent to the mesh size function in SeismicMesh. In 2D we do not use Lloyd smoothing as it can significantly increase mesh generation time (but produce higher quality cells).
* Mesh generation with Gmsh is accomplished via [pygmsh](https://github.com/nschloe/pygmsh) 7.0.0 with all default options and, similar to CGAL, an approximately equivalent cell size function is passed.
* For SeismicMesh version 3.1.4, we perform all examples with 25 meshing iterations with a psuedo-timestep of 0.30 and then run the sliver removal implemention to bound the diheral angle to 10 degrees in 3D and, in 2D, we delete any lower quality elements on the boundary (with a cell quality less than 10 percent).
* All programs are executed in a seqeuntial mode. It's important to note however that a significant speed-up can be achieved for moderate to large problems using the [parallel capabilities](https://seismicmesh.readthedocs.io/en/master/tutorial.html#basics) provided in SeismicMesh. Threading based parallelism can be used with Gmsh and CGAL but these benchmarks have not been explored.
* The scripts with the prefix `run` iterate over a range of relevant problem sizes to produce the timining and quality scales at different scales.
