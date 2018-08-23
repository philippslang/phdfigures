import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatches
from matplotlib import rcParams


def plot(
    theta,
    radii,
    cprop=[],
    label="",
    cbformat=None,
    save_as=None,
    tformat="{:d}",
    as_marker=False,
    alpha=1.0,
    nbars=None
):
    """Equal angle projection.

    See for example
    Statistical geological discrete fracture network model, A Fox et al., 2007"""
    fig = plt.figure(figsize=(7, 6))
    ax = plt.subplot(111, projection="polar")

    if len(cprop) > 0:
        cda = ax.scatter(
            theta,
            radii,
            s=140,
            c=cprop,
            marker=".",
            edgecolors="none",
            cmap=plt.get_cmap("hot_r"),
            zorder=3,
            alpha=alpha,
        )
        if label:
            cbr = fig.colorbar(cda, pad=0.075, shrink=0.75, format=cbformat)
            cbr.set_label(label)
    if nbars:
        width_bars = 2.0 * np.pi / nbars
        theta_bars = np.arange(0.0, 2*np.pi, 2*np.pi/nbars)
        radii_bars = []
        for i in range(nbars):
            theta_start = width_bars * i
            theta_end = width_bars * (i + 1)
            where = np.where((theta > theta_start) & (theta <= theta_end))
            count = len(where[0])
            radii_bars.append(count)
        radii_bars = np.array(radii_bars)
        print(np.sum(radii_bars))
        radii_bars = radii_bars / radii_bars.max()
        radii_bars = radii_bars * radii.max()
        ax.bar(theta_bars, radii_bars, width=width_bars, bottom=0.0, color='grey', alpha=0.2, zorder=2)
    if len(cprop) == 0 or as_marker:
        ax.scatter(theta, radii, s=60, c="k", marker="x", alpha=0.5, zorder=1)

    ax.grid(color="gray", linewidth=1, linestyle="-", alpha=0.4)
    ax.set_xticklabels(
        ["$\mathit{x}$", "", "$\mathit{y}$", "", "$\mathit{x}$", "", "$\mathit{y}$"]
    )
    ax.set_yticklabels([])
    ax.set_rmax(1.01)
    fig.text(0.02, 0.9, "{0} = {1:d}".format("$\mathit{N}$", len(theta)))
    plt.tight_layout()
    if save_as:
        plt.savefig(save_as)
    plt.show()

