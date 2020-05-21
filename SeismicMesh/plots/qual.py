import numpy as np

import matplotlib.pyplot as plt


def dh_angles_hist(dh_angles):
    """
    Histograms of dihedral angles in tetrahedrons
    """

    n, bins, patches = plt.hist(
        x=dh_angles, bins="auto", color="#0504aa", alpha=0.7, rwidth=0.85
    )
    maxN = np.amax(n)
    minDh = np.amin(dh_angles)
    maxDh = np.amax(dh_angles)
    plt.plot(np.linspace(minDh, minDh, 10), np.linspace(0, maxN / 2, 10), "r-")
    plt.plot(np.linspace(maxDh, maxDh, 10), np.linspace(0, maxN / 2, 10), "r-")
    plt.text(minDh, maxN, minDh)
    plt.text(maxDh, maxN, maxDh)
    frame1 = plt.gca()
    frame1.axes.yaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticks([])
    # Hide the right and top spines
    frame1.spines["right"].set_visible(False)
    frame1.spines["top"].set_visible(False)
    plt.show()


def re_ratios_hist(re_ratios):
    """
    Histograms of radius-edge metric
    """
    n, bins, patches = plt.hist(
        x=re_ratios, bins="auto", color="#0504aa", alpha=0.7, rwidth=0.85
    )
    maxN = np.amax(n)
    minRe = np.amin(re_ratios)
    maxRe = np.amax(re_ratios)
    plt.plot(np.linspace(minRe, minRe, 10), np.linspace(0, maxN / 2, 10), "r-")
    plt.plot(np.linspace(maxRe, maxRe, 10), np.linspace(0, maxN / 2, 10), "r-")
    plt.text(minRe, maxN, minRe)
    plt.text(maxRe, maxN, maxRe)
    frame1 = plt.gca()
    frame1.axes.yaxis.set_ticklabels([])
    frame1.axes.yaxis.set_ticks([])
    # Hide the right and top spines
    frame1.spines["right"].set_visible(False)
    frame1.spines["top"].set_visible(False)
    plt.show()
