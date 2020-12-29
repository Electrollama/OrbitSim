from numpy import linspace, cos, sin
from math import pi, atan2
from matplotlib import pyplot as plt

"""
units for distance is Mm
units for time is day = 2.16e4 s on Kerbin
1 day = 2.16e4 s in KSP
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
    :param res: resolution [arbitrary]
    """
    n_points = round(res * (80 + (800**min(orbit.ecc, 2.0))))
    amin, amax = orbit.angle_range()
    if amax <= amin:
        amin -= 2*pi
    thetas = linspace(amin, amax, n_points)
    xlist, ylist = [], []
    for a in thetas:
        r = orbit.hdist_theta(a)
        xlist.append(r * cos(a))
        ylist.append(r * sin(a))
    plt.plot(xlist, ylist, ls=ls, color=color)
    #plot markers for periapsis and apoapsis
    x_peri, y_peri = polar_to_xy(orbit.peri, orbit.lop)
    if marker is not None:
        if orbit.is_hyp:
            x_marker, y_marker = [x_peri, y_peri]
        else:
            x_apo, y_apo = polar_to_xy(orbit.hdist_ta(pi), orbit.lop + pi)
            x_marker, y_marker = [x_peri, x_apo], [y_peri, y_apo]
        plt.plot(x_marker, y_marker, color=color, marker=marker, ls='None')

def plot_body(body, xpos=0, ypos=0, labels=True):
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
    plt.plot([xpos, ], [ypos, ], color=body.color, marker='.')
    plt.text(xpos, ypos, body.name,
             va='bottom', ha='left', color='white')

def plot_system(system, time=0.0, labels=True):
    """
    Plot a parent body and its satellites to "plt" at a particular time. [Doesn't show]
    :param system: Body object for the system
    :param time: time in days
    :param labels: Bool, adds the body name text
    """
    Sats = system.satellites
    print('Plotting satellites...')
    for planet in Sats.keys():
        p = Sats[planet]
        plot_orbit(p.orbit, color=p.color)
        x, y = polar_to_xy(*p.orbit.position(time))
        print(p.name)
        plot_body(p, x, y, labels)
    print('Plotting parent...')
    plot_body(system)
    print('Displaying plot...')
    plt.title('Time: {}'.format(round(time, 1)))

def plot_show():
    """
    Standard orbit plot formatting and showing.
    """
    ax = plt.gca()
    ax.patch.set_facecolor('#0e0a16')
    plt.axis('equal')
    plt.show()

xy_to_polar = lambda x, y: ((x**2 + y**2)**0.5, atan2(y, x))
polar_to_xy = lambda r, theta: (r*cos(theta), r*sin(theta))