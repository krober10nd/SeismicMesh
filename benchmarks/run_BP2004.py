# run BP2004 benchmark and plot a figure
import matplotlib.pyplot as plt
import numpy
from benchmark_BP2004 import run_gmsh, run_SeismicMesh, _build_sizing


colors1 = ["ko-", "ro-"]
colors2 = ["ko--", "ro--"]
labels = ["gmsh", "SeismicMesh"]

plt.rcParams.update({"font.size": 18})

entries = []
# minimize mesh size
rg = numpy.linspace(75.0, 25.0, 5)
freqs = [2, 3, 4, 5, 6]
for i, func in enumerate([run_gmsh, run_SeismicMesh]):
    q = []
    mq = []
    nv = []
    nc = []
    t = []

    for hmin, freq in zip(rg, freqs):
        ef = _build_sizing(HMIN=hmin, FREQ=freq)
        quality, elapsed, num_vertices, num_cells = func(ef, HMIN=hmin)
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
plt.legend()
plt.xlabel("Number of cells")
plt.ylabel("Elapsed time (s)")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
plt.grid()

plt.subplot(1, 2, 2)
plt.title("Number of cells vs. cell quality")
plt.xlabel("Number of cells")
plt.ylabel("Cell quality")
plt.xlabel("Number of cells")
plt.ylim(ymax=1, ymin=0)
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
plt.grid()

plt.show()
