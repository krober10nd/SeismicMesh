# run EAGE benchmark and plot a figure
import matplotlib.pyplot as plt
import numpy
from benchmark_EAGE import run_gmsh, run_SeismicMesh, _build_sizing


colors1 = ["ko-", "ro-", "bo-"]
colors2 = ["ko--", "ro--", "bo--"]
labels = ["gmsh", "SeismicMesh", "cgal"]

plt.rcParams.update({"font.size": 18})

entries = []
# minimal mesh size
rg = numpy.linspace(300, 150.0, 5)
freqs = numpy.linspace(1, 2, 5)

# for i, func in enumerate([run_gmsh, run_SeismicMesh, run_cgal]):
for i, func in enumerate([run_gmsh, run_SeismicMesh]):
    q = []
    mq = []
    nv = []
    nc = []
    t = []

    for hmin, freq in zip(rg, freqs):
        ef = _build_sizing(HMIN=hmin, FREQ=freq)
        _, quality, elapsed, num_vertices, num_cells = func(
            ef,
            HMIN=hmin,
        )
        q.append(numpy.mean(quality))
        mq.append(numpy.min(quality))
        t.append(elapsed)
        nv.append(num_vertices)
        nc.append(num_cells)

    plt.subplot(1, 2, 1)
    h = plt.plot(nc, t, colors1[i], label=labels[i])

    plt.subplot(1, 2, 2)
    plt.plot(nc, q, colors1[i])
    plt.plot(nc, mq, colors2[i])
    entries.append(h)

plt.subplot(1, 2, 1)
plt.title("Number of cells vs. mesh generation time")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
plt.legend()
plt.xlabel("Number of cells")
plt.ylabel("Elapsed time (s)")
plt.grid()

plt.subplot(1, 2, 2)
plt.title("Number of cells. vs cell quality")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
plt.xlabel("Number of cells")
plt.ylabel("Cell quality")
plt.grid()

plt.show()
