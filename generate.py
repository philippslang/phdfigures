#!/usr/bin/env python3

import json
import os
import sys
import collections

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import figtools


HMRunsData = collections.namedtuple(
    "HMRunsData",
    "kmatrix cprops radii theta theta_kmax radii_kmax theta_kmed radii_kmed theta_kmin radii_kmin kmax kmed kmin",
)


def parse_array(data, key):
    return np.array(data[key])


def parse_hmruns_data(data):
    kmatrix = data["kmatrix"]
    cprops = data["cprops"]
    radii = parse_array(data, "radii")
    theta = parse_array(data, "theta")
    theta_kmax = parse_array(data, "theta_kmax")
    radii_kmax = parse_array(data, "radii_kmax")
    theta_kmed = parse_array(data, "theta_kmed")
    radii_kmed = parse_array(data, "radii_kmed")
    theta_kmin = parse_array(data, "theta_kmin")
    radii_kmin = parse_array(data, "radii_kmin")
    kmax = parse_array(data, "kmax")
    kmed = parse_array(data, "kmed")
    kmin = parse_array(data, "kmin")
    return HMRunsData(
        kmatrix,
        cprops,
        radii,
        theta,
        theta_kmax,
        radii_kmax,
        theta_kmed,
        radii_kmed,
        theta_kmin,
        radii_kmin,
        kmax,
        kmed,
        kmin,
    )


def parse_contactarea_data(data):
    areas = data["areas"]
    return [np.array(a) for a in areas]


def plot_hemispheres(data, nbars=None):
    fnameapp = ""
    if nbars is not None:
        fnameapp = "_bars"
    figtools.hemisphere.plot(
            data.theta_kmax,
            data.radii_kmax,
            np.log10(kmax_n),
            label=r"$\mathrm{log}_\mathrm{10}\mathit{k_{max}^\prime}$",
            cbformat="%.1f",
            save_as=os.path.join(".", "figures", "hemisphere_kmax{}.png".format(fnameapp)),
            nbars=nbars,
            alpha=0.8,
        )
    figtools.hemisphere.plot(
            data.theta,
            data.radii,
            ah_prime,
            label=r"$\mathit{a_{h}^\prime}$",
            cbformat="%.4f",
            save_as=os.path.join(".", "figures", "hemisphere_ahprime{}.png".format(fnameapp)),
            nbars=nbars,
            alpha=0.8,
        )


if __name__ == "__main__":

    # hemisphere plots
    if 1:
        with open(os.path.join(".", "data", "hmruns.json")) as f:
            data = parse_hmruns_data(json.load(f))

        kmax_n = data.kmax / data.kmatrix
        kmed_n = data.kmed / data.kmatrix
        kmin_n = data.kmin / data.kmatrix
        ah_prime = np.array(data.cprops["ah"]) / (2 * np.array(data.cprops["radius"]))

        from matplotlib import rc
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        ## for Palatino and other serif fonts use:
        #rc('font',**{'family':'serif','serif':['Palatino']})
        rc('text', usetex=True)

        plot_hemispheres(data, nbars=None)
        plot_hemispheres(data, nbars=8)

    # contact size plot
    if 0:
        with open(os.path.join(".", "data", "contactareas.json")) as f:
            allareas = parse_contactarea_data(json.load(f))

        lsim = len(allareas) - 1
        colors, alphas = ["blue"] + ["red"] * lsim, [1.0] + [0.3] * lsim

        def hook(ax):
            i, exponent = 120000, -0.6
            def powerlaw(x):
                return i * x**exponent
            slopeAx = np.array([1.6E5, 1.5E7])
            slopeAy = powerlaw(slopeAx)
            blue_line = mpl.lines.Line2D([], [], color="blue", label="Nemoto et al., 2009")
            red_line = mpl.lines.Line2D([], [], color="red", alpha=0.3, label="Numerical")
            fit_line = mpl.lines.Line2D([], [], color="black", linestyle='--', label="120000 x$^{-0.6}$")
            plt.legend(handles=[blue_line, red_line, fit_line])
            ax.plot(slopeAx, slopeAy, linestyle='--', lw=3, color="black")

        figtools.contactdistribution.plot(
            allareas,
            figsize=(8, 4),
            colors=colors,
            alphas=alphas,
            hook=hook,
            save_as=os.path.join(".", "figures", "contact_num_exp.png"),
        )
