import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatches
import json
from matplotlib import rcParams


def plot(
    theta,
    radii,
    cprop=[],
    label="",
    as_text=False,
    cbformat=None,
    save_as=None,
    tformat="{:d}",
    as_marker=False,
    alpha=1.0,
):
    """Equal angle projection.

    See for example
    Statistical geological discrete fracture network model, A Fox et al., 2007"""
    fig = plt.figure(figsize=(7, 6))
    ax = plt.subplot(111, projection="polar")
    if not as_text:
        if len(cprop) > 0:
            cda = ax.scatter(
                theta,
                radii,
                s=80,
                c=cprop,
                marker=".",
                edgecolors="none",
                cmap=plt.get_cmap("hot_r"),
                zorder=2,
                alpha=alpha,
            )
            if label:
                cbr = fig.colorbar(cda, pad=0.075, shrink=0.75, format=cbformat)
                cbr.set_label(label)
        if len(cprop) == 0 or as_marker:
            ax.scatter(theta, radii, s=60, c="k", marker="x", alpha=0.5, zorder=1)
    else:
        [
            ax.text(theta[i], radii[i], tformat.format(cprop[i]), fontsize=8)
            for i in range(len(theta))
        ]
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

