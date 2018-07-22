from math import (sin, cos, tan, pi, exp, acos,
                  sinh, cosh, atan2, isnan, acosh,
                  atan)
from newton_method import zero

epsilon = 1e-6 #Precision level
is_zero = lambda x: abs(x) < epsilon
"""
lop: angle of periapsis
t0: a time when the object was at periapsis
If retrograde, negative angles when plotting and converting to state,
   and in calculating angular momentum (negative).
angles are in radians and in the prograde direction, except in position and state
"""


class Orbit:
    """
    A simple conic orbit based on modified 2D orbital elements:
    eccentricity - Elliptic: 0 - 1.0, Parabolic: 1.0, Hyperbolic: > 1.0
    periapsis [Mm] - closest approach distance from the parent's center
    t0 [days] - time since some particular periapsis
    lop [rad] - angle of periapsis from some reference angle
    mu [Mm^2/day^2] - Standard gravitational parameter for parent.
        If negative when initialized, the orbit is retrograde
    """
    is_hyp = False
    is_parab = False
    padding = 0.0  # for min and max plotted angle

    def __init__(self, load_dict):
        """
        Do not call directly. Instead call from mk_Orbit
        :param load_dict: keys ecc, peri, t0, lop, mu
        """
        self.ecc = load_dict['ecc']
        self.peri = load_dict['peri']
        self.t0 = load_dict['t0']
        self.lop = load_dict['lop'] * pi / 180.0
        if load_dict['mu'] < 0:  # placeholder to show retrograde
            self.mu = -load_dict['mu']
            self.retro = True
        else:
            self.mu = load_dict['mu']
            self.retro = False

    def __repr__(self):
        result = ''
        size = (round(param) for param in [self.peri, self.apo()])
        result += '<ORBIT ({},{})>'.format(*size)
        return result

    def Elements(self):
        """
        Write a dictionary with orbital elements
        """
        if self.retro:
            mu = self.mu * -1.0
        else:
            mu = self.mu
        return {'ecc': self.ecc,
                'peri': self.peri,
                'lop': self.lop*180/pi,
                't0': self.t0,
                'mu': self.mu}

    def get_element(self, param):
        """
        Obtain an element from the orbit.
        """
        return self.Elements()[param]

    def set_element(self, param, new_val):
        """
        Change an element in the orbit
        :param param: name of the parameter to change
        :param new_val: new value
        :return: a new orbit with the changed element
        """
        Elems = self.Elements()
        Elems[param] = new_val
        return mk_Orbit(Elems)

    def apo(self):
        """apoapsis"""
        return self.semi() * (1 + self.ecc)

    def semi(self):
        """Semimajor axis"""
        if self.is_parab:
            return -1 * self.peri / epsilon  # infinity
        result = self.peri / (1.0 - self.ecc)
        return result  # hyperbolic : negative answer

    def period(self):
        """Orbital period"""
        return 2 * pi * (self.semi() ** 3 / self.mu) ** 0.5

    def mean_anom(self, t):
        """mean anomaly at some time"""
        t = t - self.t0
        anom = t / self.period()
        if not self.is_hyp:
            anom = anom % 1
        return 2 * pi * anom

    def ecc_anom(self, t):
        """Eccentric anomaly"""
        m = self.mean_anom(t)
        f = lambda e: e - self.ecc * sin(e) - m
        return zero(f, guess=m, tol=pi * epsilon)

    def true_anom(self, t):
        """Angle between periapsis and current position"""
        ea = self.ecc_anom(t)
        cx = (1 - self.ecc) ** 0.5 * cos(ea / 2)
        cy = (1 + self.ecc) ** 0.5 * sin(ea / 2)
        ta = atan2(cy, cx) * 2
        return ta % (2 * pi)

    def ang_momentum(self):
        """Orbit angular momentum (for calculations)"""
        ang_momentum = (self.mu * self.semi() * (1 - self.ecc ** 2.0)) ** 0.5
        if self.retro:
            return -ang_momentum
        else:
            return ang_momentum

    def slr(self):
        """Semi-latus Rectum (for calculations)"""
        return self.semi() * (1 - self.ecc ** 2)

    def hdist_ta(self, anom):
        """distance for a particular true anomaly"""
        return self.slr() / (1 + self.ecc * cos(anom))

    def hdist_theta(self, theta):
        """distance for a particular polar angle"""
        return self.hdist_ta(theta - self.lop)

    def hdist(self, t):
        """distance from the parent body's center at a particular time"""
        ta = self.true_anom(t)
        return self.hdist_ta(ta)

    def angle_range(self):
        """Range of angles for polar plotting"""
        return 0.0, 2 * pi

    def position(self, t):
        """Polar coordinates at some time"""
        ta = self.true_anom(t)
        theta = ta + self.lop
        if self.retro:
            theta *= -1
        r, anom = self.hdist_ta(ta), theta % (2 * pi)
        if isinstance(r, complex) or isinstance(anom, complex):
            raise ValueError('Complex value encountered in position.')
        return r, anom

    def get_state(self, t):
        """For a particular time, get the object's state vector."""
        mu = self.mu
        r, long = self.position(t)
        v = ((2.0 / r - 1.0 / self.semi()) * mu) ** 0.5
        temp = self.ang_momentum() / (r * v)
        direction = acos(min(max(temp, -1.0), 1.0))
        if self.true_anom(t) > pi:
            direction *= -1  # returning to periapsis
        state_dict = {'r': r,
                      'speed': v,
                      'dir': direction* 180 / pi,
                      'long': long * 180 / pi,  # polar angle
                      'mu': mu,
                      'time': t}
        if isinstance(r, complex) or isinstance(v, complex):
            raise ValueError
        return State(state_dict)


class Hyperbolic(Orbit):
    """
    A special case for hyperbolic orbits, requiring some different calculations
    and considerations.
    """
    is_hyp = True
    padding = pi / 25

    def __init__(self, load_dict):
        Orbit.__init__(self, load_dict)

    def __repr__(self):
        result = ''
        peri = round(self.peri)
        angle = round(self.esc_angle(), 1)
        result += '<FLYBY:({},{})>'.format(peri, angle)
        return result

    def esc_angle(self):
        """Angle between incoming and outgoing velocity, relative to the parent"""
        return 2 * acos(-1 / self.ecc)

    def period(self):
        """Orbital period"""
        return 2 * pi * (-self.semi() ** 3 / self.mu) ** 0.5

    def ecc_anom(self, t):
        """Eccentric anomaly"""
        f = lambda e: self.ecc * sinh(e) - e - self.mean_anom(t)
        return zero(f, guess=self.mean_anom(t),
                    tol=pi * epsilon)

    def true_anom(self, t):
        """Angle between the current position and periapsis"""
        ea = self.ecc_anom(t)
        cx = (self.ecc - 1) ** 0.5 * cosh(ea / 2)
        cy = (self.ecc + 1) ** 0.5 * sinh(ea / 2)
        ta = atan2(cy, cx) * 2
        return ta % (2 * pi)

    def ang_momentum(self):
        """Angular momentum (for calculations)"""
        if self.semi() >= 0:
            raise TypeError('Non-hyperbolic treated as hyperbolic')
        result = (self.mu * self.slr()) ** 0.5
        if self.retro:
            result *= -1
        return result

    def angle_range(self):
        """Polar angle range for plotting. Note a hyperbola does not
        occupy all angles."""
        theta = self.esc_angle()
        amin = (self.lop - theta / 2 + self.padding) % (2 * pi)
        amax = (self.lop + theta / 2 - self.padding) % (2 * pi)
        if self.retro:
            return -amax, -amin
        else:
            return amin, amax


class Parabolic(Hyperbolic):
    """
    A special case for nearly parabolic orbits to include some
    needed, some handy equation modifications.
    """
    is_parab = True
    padding = pi / 10

    def __init__(self, load_dict):
        Orbit.__init__(self, load_dict)

    def semi(self):
        """Semimajor axis"""
        return None  # -1*self.peri / epsilon  # infinity

    def ang_momentum(self):
        """orbit angular momentum (for calculations)"""
        result = (self.mu * self.slr()) ** 0.5
        if self.retro:
            result *= -1
        return result

    def slr(self):
        """Semi-latus Rectum (for calculations)"""
        return 2 * self.peri

    def true_anom(self, t):
        """Angle between the current position and periapsis,
        using Barker's Equation for an analytic solution."""
        A = 3 / 2 * (self.mu / (2 * self.peri ** 3)) ** 0.5 * (t - self.t0)
        B = (A + (A ** 2 + 1) ** 0.5) ** (1 / 3)
        return 2 * atan(B - 1 / B) % (2 * pi)

    def angle_range(self):
        """Range of angles for polar plotting. Note that a parabolic
        orbit occupies all angles, but the altitude is very large
        away from periapsis."""
        amin = (self.lop - pi + self.padding) % (2 * pi)
        amax = (self.lop + pi - self.padding) % (2 * pi)
        return amin, amax


class State:
    """
    The orbital state vectors representing a location along an orbit.
    r - distance [Mm] from the parent's center
    speed - magnitude of velocity [Mm/day]
    dir - angle of velocity [rad]
    long - longitude (angle) [rad] of position relative to the
        reference direction
    mu - parent's standard gravitational parameter [Mm^2/day^2]
    time - time at which the position was taken
    """

    def __init__(self, load_dict):
        self.r = load_dict['r']
        self.speed = load_dict['speed']
        self.dir = load_dict['dir']* pi / 180
        self.long = load_dict['long'] * pi / 180
        self.mu = load_dict['mu']
        self.time = load_dict['time']

    def __repr__(self):
        r = round(self.r)
        long = round(self.long * 180.0 / pi)
        vel = round(self.speed)
        return '<STATE ({}, {}) v: {}>'.format(r, long, vel)

    def Elements(self):
        """
        Write a dictionary with state elements
        """
        return {'r': self.r,
                'speed': self.speed,
                'dir': self.dir * 180.0 / pi,
                'long': self.long * 180 / pi,
                'mu': self.mu,
                'time': self.time}

    def get_element(self, param):
        """
        Obtain an element from the state.
        """
        return self.Elements()[param]

    def set_element(self, param, new_val):
        """
        Change an element in the state
        :param param: name of the parameter to change
        :param new_val: new value
        :return: a new state with the changed element
        """
        Elems = self.Elements()
        Elems[param] = new_val
        return State(Elems)

    def ecc_vector(self):
        """Eccentricity vector in polar (for calculations)"""
        rx, ry, vx, vy = (self.r * cos(self.long), self.r * sin(self.long),
                          - self.speed * sin(self.long - self.dir),
                          self.speed * cos(self.long - self.dir))
        k0 = rx * vx + ry * vy
        ecc_x = ((self.speed ** 2 - self.mu / self.r) * rx - k0 * vx) / self.mu
        ecc_y = ((self.speed ** 2 - self.mu / self.r) * ry - k0 * vy) / self.mu
        return xy_to_polar(ecc_x, ecc_y)

    def get_orbit(self):
        """
        Finds the orbit corresponding to the state.
        Based on Fundamentals of Astrodynamics and Applications, by Vallado, 2007.
        :return: Orbit object
        """
        orbit = dict()
        self.dir = self.dir % (2.0 * pi)
        if abs(pi - self.dir) < pi / 2:  # test if retrograde; make corrections
            self.dir = pi - self.dir
            self.long = -self.long
            orbit['mu'] = self.mu * -1
        else:
            orbit['mu'] = self.mu
        orbit['ecc'], orbit['lop'] = self.ecc_vector()
        orbit['lop'] = orbit['lop'] * 180 / pi
        energy = (self.speed ** 2 / 2.0) - (self.mu / self.r)
        k1 = -2.0 * energy / self.mu  # 1 / semi
        if is_zero(orbit['ecc']):  # Circular
            orbit['lop'] = self.long*180/pi
            orbit['peri'] = 1 / k1  # semimajor axis
            orbit['ecc'] = 0.0
            orbit['t0'] = self.time
        elif is_zero(k1):  # if parabolic
            ang_momentum = self.r * self.speed * cos(self.dir)
            # print('Parabolic')
            orbit['peri'] = ang_momentum ** 2 / (2 * self.mu)  # parabolic only
            ta = self.long - orbit['lop']
            # Apply Barker's equation
            D = tan(ta / 2)
            dt = 0.5 * (ang_momentum ** 3 / self.mu ** 2) * (D + D ** 3 / 3)
            orbit['t0'] = self.time - dt
        elif orbit['ecc'] > 1: # hyperbolic
            semi = 1 / k1
            orbit['peri'] = semi * (1 - orbit['ecc'])
            # convert to mean anomaly
            k2 = self.r * k1
            ecc_anom = acosh((1 - k2) / orbit['ecc'])
            if self.dir < 0 or self.dir >= pi * 3 / 4:
                ecc_anom *= -1
            m_anom = orbit['ecc'] * sinh(ecc_anom) - ecc_anom  # in revolutions
            n = (-self.mu / semi ** 3) ** 0.5
            orbit['t0'] = self.time - m_anom / n
        else:  # Elliptical
            semi = 1.0 / k1
            orbit['peri'] = semi * (1 - orbit['ecc'])
            # convert to mean anomaly
            k2 = self.r * k1
            temp = (1.0 - k2) / orbit['ecc']
            temp = min(max(temp, -1.0), 1.0)
            ecc_anom = acos(temp)  # range 0 to pi
            if self.dir < 0 or self.dir >= pi * 3 / 4:  # test if returning to periapsis
                ecc_anom *= -1
            m_anom = ecc_anom - orbit['ecc'] * sin(ecc_anom)
            period = 2 * pi * (semi ** 3 / self.mu) ** 0.5
            progress = m_anom / (2 * pi)
            orbit['t0'] = self.time - progress * period
        return mk_Orbit(orbit)

    def vel_vect(self):
        """Velocity vector in polar"""
        return self.speed, self.long + pi / 2 - self.dir

    def get_relative(self, parent):
        """
        The state from another state's reference frame
        :param parent: The state of the reference frame
        :return: A new State object
        """
        result = {}
        rel_pos = add_polar(self.r, parent.r,
                            self.long, parent.long,
                            subtr=True)
        result['r'], result['long'] = rel_pos
        v1, av1 = self.vel_vect()
        v2, av2 = parent.vel_vect()
        rel_vel = add_polar(v1, v2, av1, av2, subtr=True)
        result['speed'], result['dir'] = rel_vel
        result['time'] = self.mu
        result['mu'] = self.mu
        return State(result)

    def get_dist(self, parent):
        """
        Distance from a reference state
        :param parent: the reference state object
        :return: Distance [Mm]
        """
        return (self.r ** 2 + parent.r ** 2 -
                2 * self.r * parent.r *
                cos(self.long - parent.long)) ** 0.5

    def boost(self, dv_vect):
        """
        Add a velocity to the current velocity vector
        :param dv_vect: change in velocity (polar)
        Change the state's velocity vector
        """
        v, dir = self.speed, self.dir
        vx, vy = v * cos(dir), v * sin(dir)
        dv, dv_dir = dv_vect
        dvx, dvy = dv * cos(dv_dir), dv * sin(dv_dir)
        vx += dvx
        vy += dvy
        new_speed = (vx ** 2 + vy ** 2) ** 0.5
        new_dir = atan2(vy, vx)
        self.speed = new_speed
        self.dir = new_dir


def mk_Orbit(load_dict):
    """
    Correctly initialize an Orbit object.
    :param load_dict: [See Orbit.__init__]
    :return: An Orbit, Hyperbolic or Parabolic object
    """
    if load_dict is None:
        return
    if is_zero(load_dict['ecc'] - 1.0):  # if parabolic
        result = Parabolic(load_dict)
        return result
    elif load_dict['ecc'] > 1:
        return Hyperbolic(load_dict)
    else:
        return Orbit(load_dict)

def add_polar(r1, r2, a1, a2, subtr=False):
    """Add two polar vectors"""
    if subtr:
        a2 += pi
    a = a2 - a1
    dr = (r1 ** 2 + r2 ** 2 +
          2 * r1 * r2 * cos(a)) ** 0.5
    dy = -r2 * sin(a)
    dx = r1 + r2 * cos(a)
    da = a1 + atan2(dy, dx)
    return dr, da

xy_to_polar = lambda x, y: ((x ** 2 + y ** 2) ** 0.5, atan2(y, x))
polar_to_xy = lambda r, theta: (r * cos(theta), r * sin(theta))
