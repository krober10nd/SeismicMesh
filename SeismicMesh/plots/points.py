import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3


def plot_points3(points, marker="gx", hold_on=False):
    """
    plot points in 3d
    """
    axes = a3.Axes3D(plt.figure())
    axes.plot(points[:, 0], points[:, 1], points[:, 2], marker)
    axes.set_axis_off()
    if hold_on is False:
        plt.show()
    return None
