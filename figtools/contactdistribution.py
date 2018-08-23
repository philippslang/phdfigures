import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatches
from matplotlib import rcParams


def plot(
    allareas,
    show=True,
    save_as="",
    bins=4,
    label="",
    figsize=None,
    colors=None,
    alphas=None,
    hook=None,
):
    """
    Log-log plot of cluster siye vs count.
    """
    fig, ax = plt.subplots(figsize=figsize)
    legend = len(label)
    for i, areas in enumerate(allareas):
        b = np.logspace(np.log10(areas.min()), np.log10(areas.max()), num=bins + 1)
        hist, bin_edges = np.histogram(areas, bins=b)
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.
        color, alpha = None, 1.0
        if colors is not None:
            color = colors[i]
        if alphas is not None:
            alpha = alphas[i]
        lab = ""
        try:
            lab = L[i]
        except:
            pass
        ax.loglog(bin_centers, hist, label=lab, color=color, alpha=alpha)

    ax.set_xlabel("Area of contacts ($\mu$m$^2$)")
    ax.set_ylabel("Frequency")
    if legend:
        ax.legend()
    if hook is not None:
        hook(ax)
    plt.tight_layout()
    if len(save_as) != 0:
        plt.savefig(save_as)
    if show:
        plt.show()
    plt.close()
