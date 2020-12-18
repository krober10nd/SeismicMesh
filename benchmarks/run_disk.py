# run disk benchmark and plot a figure
import matplotlib
import matplotlib.pyplot as plt

import numpy
from benchmark_disk import run_gmsh, run_SeismicMesh, run_cgal

matplotlib.use("TkAgg")

# plt.rcParams.update({"font.size": 14})


colors1 = ["C7o-", "C3o-", "C0o-"]
colors2 = ["C7o--", "C3o--", "C0o--"]
labels = ["Gmsh", "SeismicMesh", "CGAL"]

entries = []
# minimize mesh size
rg = numpy.linspace(0.1, 0.01, 10)
for i, func in enumerate([run_gmsh, run_SeismicMesh, run_cgal]):
    q = []
    mq = []
    nv = []
    nc = []
    t = []

    for hmin in rg:
        quality, elapsed, num_vertices, num_cells = func(HMIN=hmin)
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
plt.title("# of cells vs. mesh generation time [s]")
plt.legend()
plt.xlabel("# of cells")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
axes = plt.gca()
axes.yaxis.grid()
axes.set_frame_on(False)

plt.subplot(1, 2, 2)
plt.title("# of cells. vs cell quality")
plt.xlabel("# of cells")
plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
plt.ylim(ymax=1.1, ymin=0)
axes = plt.gca()
axes.yaxis.grid()
axes.set_frame_on(False)


plt.tight_layout(pad=3.0)
plt.savefig("out.svg", transparent=True, bbox_inches="tight", pad_inches=0)
