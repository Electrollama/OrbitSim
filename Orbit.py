from math import (
    sin, cos, tan, pi, exp, acos,
    sinh, cosh, atan2, isnan, acosh,
    atan, sqrt
)

epsilon = 1e-6  # Precision level
is_zero = lambda x: abs(x) < epsilon

try:
    from scipy.optimize import fsolve as _scipy_fsolve, minimize_scalar as _scipy_minimize_scalar
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


def solve_kepler_elliptic(m, ecc, tol=1e-12, max_iter=50):
    """
    Solve Kepler's equation for eccentric anomaly E:
    E - ecc * sin(E) = m
    Uses Halley's method for rapid, third-order convergence.
    """
    m = m % (2 * pi)
    # Good initial guess
    if ecc < 0.8:
        e_val = m
    else:
        e_val = pi if m > pi else m + ecc * sin(m)

    for _ in range(max_iter):
        f = e_val - ecc * sin(e_val) - m
        if abs(f) < tol:
            return e_val
        f_prime = 1.0 - ecc * cos(e_val)
        f_double_prime = ecc * sin(e_val)
        # Halley step
        delta = f / (f_prime - 0.5 * f * f_double_prime / f_prime)
        e_val -= delta
        if abs(delta) < tol:
            return e_val
    return e_val


def solve_kepler_hyperbolic(m, ecc, tol=1e-12, max_iter=60):
    """
    Solve Kepler's equation for hyperbolic anomaly F:
    ecc * sinh(F) - F = m
    """
    # Initial guess
    if abs(m) < 1.0:
        f_val = m / (ecc - 1.0)
    else:
        f_val = (1.0 if m >= 0 else -1.0) * (acosh(abs(m) / ecc + 1.0))

    for _ in range(max_iter):
        f = ecc * sinh(f_val) - f_val - m
        if abs(f) < tol:
            return f_val
        f_prime = ecc * cosh(f_val) - 1.0
        f_double_prime = ecc * sinh(f_val)
        # Halley step
        denom = f_prime - 0.5 * f * f_double_prime / f_prime
        delta = f / denom if abs(denom) > 1e-14 else f / f_prime
        f_val -= delta
        if abs(delta) < tol:
            return f_val
    return f_val


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
    padding = 0.0

    def __init__(self, load_dict):
        self.ecc = float(load_dict['ecc'])
        self.peri = float(load_dict['peri'])
        self.t0 = float(load_dict['t0'])
        self.lop = float(load_dict['lop']) * pi / 180.0
        if load_dict['mu'] < 0:
            self.mu = float(-load_dict['mu'])
            self.retro = True
        else:
            self.mu = float(load_dict['mu'])
            self.retro = False

    def __repr__(self):
        size = (round(param) for param in [self.peri, self.apo()])
        return '<ORBIT ({},{})>'.format(*size)

    def Elements(self):
        mu = -self.mu if self.retro else self.mu
        return {
            'ecc': self.ecc,
            'peri': self.peri,
            'lop': self.lop * 180.0 / pi,
            't0': self.t0,
            'mu': mu
        }

    def get_element(self, param):
        return self.Elements()[param]

    def set_element(self, param, new_val):
        elems = self.Elements()
        elems[param] = new_val
        return mk_Orbit(elems)

    def apo(self):
        return self.semi() * (1.0 + self.ecc)

    def semi(self):
        if self.is_parab:
            return -1.0 * self.peri / epsilon
        return self.peri / (1.0 - self.ecc)

    def period(self):
        return 2.0 * pi * sqrt((self.semi() ** 3) / self.mu)

    def mean_anom(self, t):
        p = self.period()
        t_eff = (t - self.t0) % p
        return 2.0 * pi * (t_eff / p)

    def ecc_anom(self, t):
        ma = self.mean_anom(t)
        return solve_kepler_elliptic(ma, self.ecc)

    def true_anom(self, t):
        ea = self.ecc_anom(t)
        cx = sqrt(max(0.0, 1.0 - self.ecc)) * cos(ea / 2.0)
        cy = sqrt(1.0 + self.ecc) * sin(ea / 2.0)
        return atan2(cy, cx) * 2.0

    def ang_momentum(self):
        h = sqrt(self.mu * self.semi() * (1.0 - self.ecc ** 2))
        return -h if self.retro else h

    def energy(self):
        return -self.mu / (2.0 * self.semi())

    def slr(self):
        return self.semi() * (1.0 - self.ecc ** 2)

    def hdist_ta(self, anom):
        denom = 1.0 + self.ecc * cos(anom)
        return self.slr() / denom

    def hdist_theta(self, theta):
        return self.hdist_ta(theta - self.lop)

    def hdist(self, t):
        ta = self.true_anom(t)
        return self.hdist_ta(ta)

    def dist_time(self, r):
        """Times at a particular altitude, return (None, None) if out of range."""
        if r < self.peri or (not self.is_hyp and r > self.apo()):
            return (None, None)
        # Exact solution via true anomaly
        # r = p / (1 + e*cos(nu)) -> cos(nu) = (p/r - 1) / e
        val = (self.slr() / r - 1.0) / self.ecc
        val = max(-1.0, min(1.0, val))
        nu = acos(val)
        # Convert to eccentric anomaly
        cx = (1.0 - self.ecc) ** 0.5 * cos(nu / 2.0)
        cy = (1.0 + self.ecc) ** 0.5 * sin(nu / 2.0)
        ea = atan2(cy, cx) * 2.0
        ma = ea - self.ecc * sin(ea)
        dt = (ma / (2.0 * pi)) * self.period()
        return (self.t0 + dt, self.t0 - dt)

    def angle_range(self):
        return 0.0, 2.0 * pi

    def position(self, t):
        ta = self.true_anom(t)
        theta = ta + self.lop
        if self.retro:
            theta = -theta
        r = self.hdist_ta(ta)
        return r, theta

    def get_state(self, t):
        r, long = self.position(t)
        v = sqrt(max(0.0, (2.0 / r - 1.0 / self.semi()) * self.mu))
        temp = self.ang_momentum() / (r * v) if (r * v) > 0 else 1.0
        direction = acos(min(max(temp, -1.0), 1.0))
        if self.true_anom(t) % (2.0 * pi) > pi:
            direction *= -1.0
        state_dict = {
            'r': r,
            'speed': v,
            'dir': direction * 180.0 / pi,
            'long': long * 180.0 / pi,
            'mu': self.mu,
            'time': t
        }
        return State(state_dict)


class Hyperbolic(Orbit):
    """
    A special case for hyperbolic orbits.
    """
    is_hyp = True
    padding = pi / 25.0

    def __init__(self, load_dict):
        super().__init__(load_dict)

    def __repr__(self):
        peri = round(self.peri)
        angle = round(self.esc_angle(), 1)
        return '<FLYBY:({},{})>'.format(peri, angle)

    def esc_angle(self):
        return 2.0 * acos(-1.0 / self.ecc)

    def period(self):
        raise ValueError("Hyperbolic orbits have undefined period.")

    def mean_anom(self, t):
        n = sqrt(-self.mu / (self.semi() ** 3))
        return n * (t - self.t0)

    def ecc_anom(self, t):
        ma = self.mean_anom(t)
        return solve_kepler_hyperbolic(ma, self.ecc)

    def true_anom(self, t):
        ea = self.ecc_anom(t)
        cx = sqrt(self.ecc - 1.0) * cosh(ea / 2.0)
        cy = sqrt(self.ecc + 1.0) * sinh(ea / 2.0)
        ta = atan2(cy, cx) * 2.0
        return ta % (2.0 * pi)

    def ang_momentum(self):
        if self.semi() >= 0:
            raise TypeError('Non-hyperbolic treated as hyperbolic')
        h = sqrt(self.mu * self.slr())
        return -h if self.retro else h

    def angle_range(self):
        theta = self.esc_angle()
        amin = (self.lop - theta / 2.0 + self.padding) % (2.0 * pi)
        amax = (self.lop + theta / 2.0 - self.padding) % (2.0 * pi)
        if self.retro:
            return -amax, -amin
        return amin, amax


class Parabolic(Hyperbolic):
    """
    Parabolic orbit using Barker's Equation.
    """
    is_parab = True
    padding = pi / 10.0

    def __init__(self, load_dict):
        super().__init__(load_dict)

    def semi(self):
        return None

    def ang_momentum(self):
        h = sqrt(self.mu * self.slr())
        return -h if self.retro else h

    def energy(self):
        return 0.0

    def slr(self):
        return 2.0 * self.peri

    def true_anom(self, t):
        a_val = 1.5 * sqrt(self.mu / (2.0 * self.peri ** 3)) * (t - self.t0)
        b_val = (a_val + sqrt(a_val ** 2 + 1.0)) ** (1.0 / 3.0)
        return 2.0 * atan(b_val - 1.0 / b_val) % (2.0 * pi)

    def angle_range(self):
        amin = (self.lop - pi + self.padding) % (2.0 * pi)
        amax = (self.lop + pi - self.padding) % (2.0 * pi)
        return amin, amax


class State:
    """
    The orbital state vectors representing a location along an orbit.
    r - distance [Mm] from parent center
    speed - magnitude of velocity [Mm/day]
    dir - flight path angle [rad]
    long - position longitude [rad]
    mu - gravitational parameter [Mm^2/day^2]
    time - epoch time [days]
    """

    def __init__(self, load_dict):
        self.r = float(load_dict['r'])
        self.speed = float(load_dict['speed'])
        self.dir = float(load_dict['dir']) * pi / 180.0
        self.long = float(load_dict['long']) * pi / 180.0
        self.mu = float(load_dict['mu'])
        self.time = float(load_dict['time'])

    def __repr__(self):
        r = round(self.r)
        long = round(self.long * 180.0 / pi)
        vel = round(self.speed)
        return '<STATE ({}, {}) v: {}>'.format(r, long, vel)

    def __add__(self, other, subtr=False):
        result = {}
        r_rel, long_rel = add_polar(self.r, other.r, self.long, other.long, subtr=subtr)
        result['r'] = r_rel
        result['long'] = long_rel * 180.0 / pi
        v1, av1 = self.vel_vect()
        v2, av2 = other.vel_vect()
        v_rel, av_rel = add_polar(v1, v2, av1, av2, subtr=subtr)
        result['speed'] = v_rel
        # Convert absolute velocity direction to flight path angle
        dir_rel = (long_rel + pi / 2.0 - av_rel) % (2.0 * pi)
        if dir_rel > pi:
            dir_rel -= 2.0 * pi
        result['dir'] = dir_rel * 180.0 / pi
        result['time'] = self.time
        result['mu'] = self.mu
        return State(result)

    def __sub__(self, other):
        return self.__add__(other, subtr=True)

    def Elements(self):
        return {
            'r': self.r,
            'speed': self.speed,
            'dir': self.dir * 180.0 / pi,
            'long': self.long * 180.0 / pi,
            'mu': self.mu,
            'time': self.time
        }

    def get_element(self, param):
        return self.Elements()[param]

    def set_element(self, param, new_val):
        elems = self.Elements()
        elems[param] = new_val
        return State(elems)

    def ecc_vector(self):
        rx = self.r * cos(self.long)
        ry = self.r * sin(self.long)
        vx = -self.speed * sin(self.long - self.dir)
        vy = self.speed * cos(self.long - self.dir)
        k0 = rx * vx + ry * vy
        ecc_x = ((self.speed ** 2 - self.mu / self.r) * rx - k0 * vx) / self.mu
        ecc_y = ((self.speed ** 2 - self.mu / self.r) * ry - k0 * vy) / self.mu
        return xy_to_polar(ecc_x, ecc_y)

    def get_orbit(self):
        orbit = dict()
        self.dir = self.dir % (2.0 * pi)
        if abs(pi - self.dir) < pi / 2.0:
            self.dir = pi - self.dir
            self.long = -self.long
            orbit['mu'] = -self.mu
        else:
            orbit['mu'] = self.mu

        orbit['ecc'], orbit['lop'] = self.ecc_vector()
        orbit['lop'] = orbit['lop'] * 180.0 / pi
        energy = (self.speed ** 2 / 2.0) - (self.mu / self.r)
        k1 = -2.0 * energy / self.mu

        rx = self.r * cos(self.long)
        ry = self.r * sin(self.long)
        vx = -self.speed * sin(self.long - self.dir)
        vy = self.speed * cos(self.long - self.dir)
        radial_v = rx * vx + ry * vy  # r . v

        if is_zero(orbit['ecc']):
            orbit['lop'] = self.long * 180.0 / pi
            orbit['peri'] = 1.0 / k1
            orbit['ecc'] = 0.0
            orbit['t0'] = self.time
        elif is_zero(k1):
            ang_momentum = self.r * self.speed * cos(self.dir)
            orbit['peri'] = ang_momentum ** 2 / (2.0 * self.mu)
            ta = self.long - orbit['lop'] * pi / 180.0
            d_val = tan(ta / 2.0)
            dt = 0.5 * (ang_momentum ** 3 / self.mu ** 2) * (d_val + (d_val ** 3) / 3.0)
            orbit['t0'] = self.time - dt
        elif orbit['ecc'] > 1.0:
            semi = 1.0 / k1
            orbit['peri'] = semi * (1.0 - orbit['ecc'])
            k2 = self.r * k1
            ecc_anom = acosh(max(1.0, (1.0 - k2) / orbit['ecc']))
            if radial_v < 0:
                ecc_anom = -abs(ecc_anom)
            else:
                ecc_anom = abs(ecc_anom)
            m_anom = orbit['ecc'] * sinh(ecc_anom) - ecc_anom
            n = sqrt(-self.mu / (semi ** 3))
            orbit['t0'] = self.time - m_anom / n
        else:
            semi = 1.0 / k1
            orbit['peri'] = semi * (1.0 - orbit['ecc'])
            k2 = self.r * k1
            temp = max(-1.0, min(1.0, (1.0 - k2) / orbit['ecc']))
            ecc_anom = acos(temp)
            if radial_v < 0:
                ecc_anom = 2.0 * pi - ecc_anom
            m_anom = ecc_anom - orbit['ecc'] * sin(ecc_anom)
            period = 2.0 * pi * sqrt(semi ** 3 / self.mu)
            progress = m_anom / (2.0 * pi)
            orbit['t0'] = self.time - progress * period
        return mk_Orbit(orbit)

    def vel_vect(self):
        return self.speed, self.long + pi / 2.0 - self.dir

    def get_dist(self, parent):
        return sqrt(
            self.r ** 2 + parent.r ** 2 -
            2.0 * self.r * parent.r * cos(self.long - parent.long)
        )

    def boost(self, dv_vect):
        v, direction = self.speed, self.dir
        vx = v * cos(direction)
        vy = v * sin(direction)
        dv, dv_dir = dv_vect
        vx += dv * cos(dv_dir)
        vy += dv * sin(dv_dir)
        self.speed = sqrt(vx ** 2 + vy ** 2)
        self.dir = atan2(vy, vx)


def mk_Orbit(load_dict):
    if load_dict is None:
        return None
    ecc = load_dict.get('ecc', 0.0)
    if is_zero(ecc - 1.0):
        return Parabolic(load_dict)
    elif ecc > 1.0:
        return Hyperbolic(load_dict)
    else:
        return Orbit(load_dict)


def add_polar(r1, r2, a1, a2, subtr=False):
    x1, y1 = r1 * cos(a1), r1 * sin(a1)
    x2, y2 = r2 * cos(a2), r2 * sin(a2)
    if subtr:
        x = x1 - x2
        y = y1 - y2
    else:
        x = x1 + x2
        y = y1 + y2
    return sqrt(x * x + y * y), atan2(y, x)


xy_to_polar = lambda x, y: (sqrt(x ** 2 + y ** 2), atan2(y, x))
polar_to_xy = lambda r, theta: (r * cos(theta), r * sin(theta))
