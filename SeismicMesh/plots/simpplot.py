import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3

from .. import geometry


def plot_tets(points, cells, hold_on=False):
    """
    vizualize tetrahedral mesh
    """
    axes = a3.Axes3D(plt.figure())
    tri = geometry.get_facets(cells)
    vts = points[tri, :]
    tri = a3.art3d.Poly3DCollection(vts)
    tri.set_alpha(0.2)
    tri.set_color("grey")
    axes.add_collection3d(tri)
    axes.plot(points[:, 0], points[:, 1], points[:, 2], "ko")
    axes.set_axis_off()
    if hold_on is not False:
        plt.show()
    return None


def plot_facets(points, facets, color="red", marker="gx", hold_on=False):
    """
    visualize the facets
    """
    axes = a3.Axes3D(plt.figure())
    vts2 = points[facets, :]
    tri2 = a3.art3d.Poly3DCollection(vts2)
    tri2.set_alpha(0.2)
    tri2.set_color(color)
    axes.add_collection3d(tri2)
    axes.plot(points[:, 0], points[:, 1], points[:, 2], marker)
    axes.set_axis_off()
    if hold_on is not False:
        plt.show()
    return None
