"""
Orbit plotter utility for OrbitSim.
"""

from math import pi, cos, sin
from Orbit import polar_to_xy

try:
    import numpy as np
    from matplotlib import pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False
    plt = None


def plot_orbit(orbit, res=1.0, color='white', ls='--', marker=None):
    if not _HAS_MPL:
        print("matplotlib not installed; skipping plot_orbit.")
        return
    n_points = round(res * (80 + 800 ** min(orbit.ecc, 2.0)))
    amin, amax = orbit.angle_range()
    if amax <= amin:
        amin -= 2.0 * pi
    thetas = np.linspace(amin, amax, n_points)
    xlist, ylist = [], []
    for a in thetas:
        r = orbit.hdist_theta(a)
        x, y = polar_to_xy(r, a)
        xlist.append(x)
        ylist.append(y)
    plt.plot(xlist, ylist, ls=ls, color=color)

    x_peri, y_peri = polar_to_xy(orbit.peri, orbit.lop)
    if marker is not None:
        plt.plot(x_peri, y_peri, color=color, marker=marker, ls='None')
        if not orbit.is_hyp:
            x_apo, y_apo = polar_to_xy(orbit.hdist_ta(pi), orbit.lop + pi)
            plt.plot(x_apo, y_apo, color=color, marker=marker, ls='None')


def plot_body(body, xpos=0, ypos=0, parent_arrow=True):
    if not _HAS_MPL:
        return
    to_scale = plt.Circle((xpos, ypos), body.size, color=body.color)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_artist(to_scale)
    plt.plot([xpos], [ypos], color=body.color, marker='o', markersize=4.5)
    plt.text(xpos, ypos, body.name, va='bottom', ha='left', color='white')


def plot_craft(craft, xpos=0, ypos=0):
    if not _HAS_MPL:
        return
    plt.plot([xpos], [ypos], color=craft.color, marker='s', markersize=3.5)
    plt.text(xpos, ypos, craft.name, va='bottom', ha='left', color='grey')


def plot_system(system, time=0.0, labels=True, show_crafts=False):
    if not _HAS_MPL:
        print("matplotlib not installed; skipping plot_system.")
        return
    sats = system.bodies
    for planet in sats.keys():
        p = sats[planet]
        plot_orbit(p.orbit, color=p.color)
        x, y = polar_to_xy(*p.orbit.position(time))
        plot_body(p, x, y, labels)
    plot_body(system)
    plt.title('Time: {}'.format(round(time, 1)))


def plot_show():
    if not _HAS_MPL:
        return
    plt.gca().patch.set_facecolor('#0e0a16')
    plt.axis('equal')
    plt.show()
