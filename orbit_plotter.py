from numpy import linspace, cos, sin
from math import pi, atan2
from matplotlib import pyplot as plt

"""
	Functions for making plots that do not depend on time

units for distance is Mm
units for time is day
1 day = 2.160e4 s in KSP
1 day = 8.64e4 s in RL
mu conversion = 4.6656e-10 in KSP
mu conversion = 7.465e-9 in RL

RL Times:   01/01/18: 6575
            05/20/18: 6714
"""

def plot_orbit(orbit, res = 1.0, color='white', ls='--', marker=None):
    """
    Plot an orbit path on "plt" where the parent is at (0,0). [doesn't show]
    :param orbit: Orbit object
    :param res: resolution [proportional]
    """
    n_points = round(res * (80 + 800**min(orbit.ecc, 2.0)))
    amin, amax = orbit.angle_range()
    if amax <= amin:
        amin -= 2*pi
    thetas = linspace(amin, amax, n_points)
    xlist, ylist = [], []
    for a in thetas:
        r = orbit.hdist_theta(a)
        x, y = polar_to_xy(r, a)
        xlist.append(x)
        ylist.append(y)
    plt.plot(xlist, ylist, ls=ls, color=color)
    #plot markers for periapsis and apoapsis
    x_peri, y_peri = polar_to_xy(orbit.peri, orbit.lop)
    if marker is not None:
        plt.plot(x_peri, y_peri, color=color, marker=marker, ls='None')
        if not orbit.is_hyp:
            x_apo, y_apo = polar_to_xy(orbit.hdist_ta(pi), orbit.lop + pi)
            plt.plot(x_apo, y_apo, color=color, marker=marker, ls='None')
    # correct plot axes if needed
    window(min(xlist), max(xlist),
           min(ylist), max(ylist))

def plot_body(body, xpos=0, ypos=0, parent_arrow=True):
    """
    Plot a body as a to-scale circle on "plt". [Doesn't show]
    :param body: Body object (for size)
    :param xpos: center x-position
    :param ypos: center y-position
    :param labels: Bool, adds the body name text
    """
    to_scale = plt.Circle((xpos, ypos), body.size, color=body.color)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_artist(to_scale)
    plt.plot([xpos, ], [ypos, ], color=body.color, marker='o', markersize=4.5)
    plt.text(xpos, ypos, body.name,
             va='bottom', ha='left', color='white')


def plot_craft(craft, xpos=0, ypos=0):
    plt.plot([xpos, ], [ypos, ], color=craft.color, marker='s', markersize=3.5)
    plt.text(xpos, ypos, craft.name,
             va='bottom', ha='left', color='grey')

def plot_system(system, time=0.0, labels=True, show_crafts=False):
    """
    Plot a parent body and its satellites to "plt" at a particular time. [Doesn't show]
    :param system: Body object for the system
    :param time: time in days
    :param labels: Bool, adds the body name text
    """
    Sats = system.bodies
    print('Plotting satellites...')
    for planet in Sats.keys():
        p = Sats[planet]
        plot_orbit(p.orbit, color=p.color)
        x, y = polar_to_xy(*p.orbit.position(time))
        print(p.name)
        plot_body(p, x, y, labels)
    print('Plotting parent...')
    plot_body(system)
    plt.title('Time: {}'.format(round(time, 1)))
    if show_crafts:
        for craft_name in system.crafts.keys():
            craft = system.crafts[craft_name]
            plot_orbit(craft.orbit, color='#2d2046')
            x, y = polar_to_xy(*craft.orbit.position(time))
            plot_body(craft, x, y)
    if system.parent is not None:
        # Arrow in the direction of the parent
        r, a = system.soi, pi + system.orbit.true_anom(time)
        plt.arrow(0.95 * r * cos(a), 0.95 * r * sin(a),  # at SOI boundary
                  0.05 * r * cos(a), 0.05 * r * sin(a),
                  color='yellow', length_includes_head=True,
                  width=0.005 * r)
        r = system.size
        plt.arrow(0, 0,  # at center
                  r * cos(a), r * sin(a),
                  color='yellow', length_includes_head=True,
                  width=0.1 * r)

def plot_show():
    """
    Standard orbit plot formatting and showing.
    """
    plt.gca().patch.set_facecolor('#0e0a16')
    plt.axis('equal')
    plt.show()


def window(new_xmin, new_xmax, new_ymin, new_ymax):
    axes = plt.gca()
    xmin, xmax = axes.get_xlim()
    ymin, ymax = axes.get_ylim()
    axes.set_xlim(min(xmin, new_xmin),
                        max(xmax, new_xmax))
    axes.set_ylim(min(ymin, new_ymin),
                        max(ymax, new_ymax))

xy_to_polar = lambda x, y: ((x**2 + y**2)**0.5, atan2(y, x))
polar_to_xy = lambda r, theta: (r*cos(theta), r*sin(theta))