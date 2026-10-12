"""
Encounters and Patched Conics Geometry for OrbitSim.

Implements shared-frame separation, search window bounding, 2D geometric MOID,
and root-finding for SOI enter and exit events.
"""

from math import cos, sin, sqrt, pi, acos, atan2, acosh, sinh
from Orbit import polar_to_xy


class EncounterGeometryError(ValueError):
    """Invalid inputs for encounter geometry."""


def _require_orbit(body, label):
    if body is None or getattr(body, 'orbit', None) is None:
        name = getattr(body, 'name', 'unnamed')
        raise EncounterGeometryError(
            '{} {!r} has no orbit.'.format(label, name)
        )


def common_parent(body_a, body_b):
    """True if both satellites share the same parent Body."""
    return (
        body_a is not None and
        body_b is not None and
        body_a.parent is not None and
        body_a.parent is body_b.parent
    )


def assert_common_parent(body_a, body_b):
    if not common_parent(body_a, body_b):
        name_a = getattr(body_a, 'name', 'body_a')
        name_b = getattr(body_b, 'name', 'body_b')
        raise EncounterGeometryError(
            'Bodies must share a parent for separation: {} and {}.'.format(name_a, name_b)
        )


def position_xy(body, t):
    """
    Inertial position in the parent body's frame (parent at origin).
    :param body: Satellite, Craft, or Body with an orbit
    :param t: time [days]
    :return: (x, y) [Mm]
    """
    _require_orbit(body, 'body')
    r, theta = body.orbit.position(t)
    return polar_to_xy(r, theta)


def separation(body_a, body_b, t):
    """
    Separation distance between two bodies sharing the same parent frame.
    :param body_a: Satellite or Craft
    :param body_b: Satellite or Craft
    :param t: time [days]
    :return: separation distance [Mm]
    """
    assert_common_parent(body_a, body_b)
    xa, ya = position_xy(body_a, t)
    xb, yb = position_xy(body_b, t)
    return sqrt((xa - xb) ** 2 + (ya - yb) ** 2)


# Alias for backwards compatibility
distance = separation


def separation_fn(body_a, body_b):
    """
    Return callable f(t) -> separation distance [Mm].
    """
    assert_common_parent(body_a, body_b)
    return lambda t: separation(body_a, body_b, t)


def synodic_period(orbit_a, orbit_b):
    """
    Synodic period for two bound Keplerian orbits about the same focus.
    :return: period [days], or None if either orbit is hyperbolic/parabolic
    """
    if orbit_a is None or orbit_b is None:
        return None
    if getattr(orbit_a, 'is_hyp', False) or getattr(orbit_b, 'is_hyp', False):
        return None
    ta = orbit_a.period()
    tb = orbit_b.period()
    diff = abs(1.0 / ta - 1.0 / tb)
    if diff < 1e-15:
        return min(ta, tb)
    return 1.0 / diff


def relative_period(orbit_a, orbit_b):
    """
    Synodic / relative period for two bound orbits.
    Alias / wrapper for synodic_period.
    """
    return synodic_period(orbit_a, orbit_b)


def ca_search_end(t0, craft, target, craft_orbits=1.0, max_days=None):
    """
    Upper time bound for a closest-approach search starting at t0.
    Uses the shorter of (craft_orbits * craft period) and synodic period.
    """
    assert_common_parent(craft, target)
    _require_orbit(craft, 'craft')
    _require_orbit(target, 'target')
    if craft_orbits <= 0:
        raise EncounterGeometryError('craft_orbits must be positive.')

    orbit_c = craft.orbit
    if orbit_c.is_hyp:
        # Hyperbolic search bounded by target period or max_days
        window = target.orbit.period() if not target.orbit.is_hyp else 30.0
    else:
        spans = [orbit_c.period() * craft_orbits]
        syn = synodic_period(orbit_c, target.orbit)
        if syn is not None:
            spans.append(syn)
        window = min(spans)

    if max_days is not None and max_days > 0:
        window = min(window, max_days)
    return t0 + window


def ca_search_interval(t0, craft, target, craft_orbits=1.0, max_days=None):
    """Return (t0, t_end) for CA search."""
    return t0, ca_search_end(t0, craft, target, craft_orbits, max_days)


# ==============================================================================
# 2D Geometric MOID (Minimum Orbit Intersection Distance)
# ==============================================================================

def geometric_moid_2d(orbit_a, orbit_b):
    """
    Calculate the 2D coplanar Minimum Orbit Intersection Distance (MOID)
    between two orbits sharing the same focus.

    For two conics r1(theta) and r2(theta), intersection equates radii:
        p1 / (1 + e1*cos(theta - lop1)) = p2 / (1 + e2*cos(theta - lop2))
    which reduces to the linear trigonometric form:
        A * cos(theta) + B * sin(theta) = C
    where:
        A = p1 * e2 * cos(lop2) - p2 * e1 * cos(lop1)
        B = p1 * e2 * sin(lop2) - p2 * e1 * sin(lop1)
        C = p2 - p1

    :param orbit_a: Orbit object
    :param orbit_b: Orbit object
    :return: dict with 'moid', 'intersects' (bool), 'thetas' (list of angles)
    """
    p1 = orbit_a.slr()
    e1 = orbit_a.ecc
    lop1 = orbit_a.lop

    p2 = orbit_b.slr()
    e2 = orbit_b.ecc
    lop2 = orbit_b.lop

    A = p1 * e2 * cos(lop2) - p2 * e1 * cos(lop1)
    B = p1 * e2 * sin(lop2) - p2 * e1 * sin(lop1)
    C = p2 - p1

    R = sqrt(A * A + B * B)

    if R < 1e-12:
        if abs(C) < 1e-12:
            return {'moid': 0.0, 'intersects': True, 'thetas': [0.0]}
        else:
            # Concentric circular orbits
            return {'moid': abs(C), 'intersects': False, 'thetas': [0.0]}

    if abs(C) <= R:
        # Analytic geometric intersections exist!
        phi = atan2(B, A)
        delta = acos(max(-1.0, min(1.0, C / R)))
        theta1 = (phi + delta) % (2.0 * pi)
        theta2 = (phi - delta) % (2.0 * pi)
        return {
            'moid': 0.0,
            'intersects': True,
            'thetas': [theta1, theta2]
        }
    else:
        # Orbits do not intersect: find the angle of closest approach
        # Sample theta around the circle to find min |r1(theta) - r2(theta)|
        best_dist = float('inf')
        best_theta = 0.0
        n_samples = 72
        for i in range(n_samples):
            th = (2.0 * pi * i) / n_samples
            r1 = orbit_a.hdist_theta(th)
            r2 = orbit_b.hdist_theta(th)
            d = abs(r1 - r2)
            if d < best_dist:
                best_dist = d
                best_theta = th

        # Fine-tune with golden section search around best_theta
        step = (2.0 * pi) / n_samples
        th_min = best_theta - step
        th_max = best_theta + step

        def obj(th_val):
            return abs(orbit_a.hdist_theta(th_val) - orbit_b.hdist_theta(th_val))

        min_d, opt_th = _golden_section_min(obj, th_min, th_max, tol=1e-7)
        return {
            'moid': min_d,
            'intersects': False,
            'thetas': [opt_th % (2.0 * pi)]
        }


def _golden_section_min(f, a, b, tol=1e-6, max_iter=60):
    """Pure-Python golden section 1D minimization."""
    invphi = (sqrt(5.0) - 1.0) / 2.0
    invphi2 = (3.0 - sqrt(5.0)) / 2.0

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return f((a + b) / 2.0), (a + b) / 2.0

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for _ in range(max_iter):
        if (b - a) < tol:
            break
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

    x_min = (a + b) / 2.0
    return f(x_min), x_min


def _bisection_root(f, a, b, tol=1e-7, max_iter=70):
    """Pure-Python root finding using bracketed bisection/secant hybrid."""
    fa = f(a)
    fb = f(b)
    if abs(fa) <= tol:
        return a
    if abs(fb) <= tol:
        return b
    if fa * fb > 0:
        return None

    for _ in range(max_iter):
        mid = (a + b) / 2.0
        fmid = f(mid)
        if abs(fmid) <= tol or (b - a) / 2.0 < tol:
            return mid
        if fa * fmid < 0:
            b = mid
            fb = fmid
        else:
            a = mid
            fa = fmid
    return (a + b) / 2.0


# ==============================================================================
# Encounter & Escape Search
# ==============================================================================

def find_encounter(craft, time, target_name, n_orbits=1.0, max_days=None):
    """
    Find the time, if it exists, that a craft enters a target body's sphere of influence.
    Uses two-stage detection:
      1. Fast geometric filtering (bounding radii + 2D geometric MOID).
      2. Root finding for d(t) = target.soi.

    :param craft: Craft or Satellite object
    :param time: starting search time [days]
    :param target_name: name of the target Body (must be child of craft.parent)
    :param n_orbits: number of orbits or synodic periods to search
    :param max_days: optional maximum search duration cap [days]
    :return: timestamp [days] of SOI entry, or None
    """
    if craft.parent is None or target_name not in craft.parent.bodies:
        return None
    target = craft.parent.bodies[target_name]

    if getattr(target.orbit, 'is_hyp', False):
        raise ValueError('Cannot complete encounter search for target with hyperbolic orbit.')

    # Gating checks on periapsis / apoapsis
    if craft.orbit.peri > target.orbit.apo() + target.soi:
        return None
    if not craft.orbit.is_hyp and craft.orbit.apo() < target.orbit.peri - target.soi:
        return None

    # Stage 1: Geometric MOID check
    moid_data = geometric_moid_2d(craft.orbit, target.orbit)
    if moid_data['moid'] > target.soi:
        return None

    # Stage 2: Time-domain root finding
    t0 = time
    t_end = ca_search_end(t0, craft, target, craft_orbits=n_orbits, max_days=max_days)
    if t_end <= t0:
        return None

    soi_radius = target.soi

    def g(t):
        return separation(craft, target, t) - soi_radius

    # Check if already inside SOI
    if g(t0) <= 0:
        return t0

    # Subdivide search interval to detect boundary crossings
    n_steps = max(40, int(80 * n_orbits))
    dt = (t_end - t0) / n_steps

    prev_t = t0
    prev_g = g(prev_t)

    for i in range(1, n_steps + 1):
        curr_t = t0 + i * dt
        curr_g = g(curr_t)

        # Entering SOI: previous separation > SOI, current separation <= SOI
        if prev_g > 0 and curr_g <= 0:
            root = _bisection_root(g, prev_t, curr_t, tol=1e-6)
            if root is not None:
                return root

        # Check for local minimum inside interval [prev_t, curr_t]
        mid_t = (prev_t + curr_t) / 2.0
        mid_g = g(mid_t)
        if mid_g < prev_g and mid_g < curr_g:
            # Local dip: check if dip penetrates SOI
            min_val, min_t = _golden_section_min(g, prev_t, curr_t, tol=1e-5)
            if min_val <= 0:
                root = _bisection_root(g, prev_t, min_t, tol=1e-6)
                if root is not None:
                    return root

        prev_t = curr_t
        prev_g = curr_g

    return None


def find_escape(sat, time, max_days=None):
    """
    Find the time, if it exists, that a craft leaves its parent's sphere of influence.
    :param sat: Craft or Satellite object
    :param time: starting time [days]
    :param max_days: maximum search window length [days]
    :return: timestamp [days] of departure, or None
    """
    parent = sat.parent
    if parent is None or parent.parent is None:
        return None

    soi_rad = parent.soi

    # Hyperbolic escape: exact analytical solution
    if sat.orbit.is_hyp:
        semi = sat.orbit.semi()  # negative
        ecc = sat.orbit.ecc
        t0 = sat.orbit.t0
        val = (1.0 - soi_rad / semi) / ecc
        if val < 1.0:
            return None
        F = acosh(val)
        M = ecc * sinh(F) - F
        n = sqrt(-sat.orbit.mu / (semi ** 3))
        dt = M / n
        t_exit = t0 + dt
        if t_exit > time:
            return t_exit
        return time

    # Bound elliptical orbit whose apoapsis is within SOI never escapes
    if sat.orbit.apo() < soi_rad:
        return None

    # Escaping elliptical orbit (apo >= soi_rad)
    val = (sat.orbit.slr() / soi_rad - 1.0) / sat.orbit.ecc
    val = max(-1.0, min(1.0, val))
    nu = acos(val)
    cx = sqrt(1.0 - sat.orbit.ecc) * cos(nu / 2.0)
    cy = sqrt(1.0 + sat.orbit.ecc) * sin(nu / 2.0)
    ea = atan2(cy, cx) * 2.0
    ma = ea - sat.orbit.ecc * sin(ea)
    dt = (ma / (2.0 * pi)) * sat.orbit.period()
    t_cand = sat.orbit.t0 + dt
    p = sat.orbit.period()
    while t_cand <= time:
        t_cand += p
    return t_cand


# ==============================================================================
# SOI Transition State Transformers
# ==============================================================================

def soi_enter(sat, time, target):
    """
    Transform a craft's orbit and parent relationship upon entering target's SOI.
    :param sat: Craft or Satellite
    :param time: time [days] of entry
    :param target: Body being entered
    """
    old_state = sat.orbit.get_state(time)
    target_state = target.orbit.get_state(time)
    new_state = old_state - target_state
    new_state.mu = target.mu

    sat.orbit = new_state.get_orbit()

    # Re-link parents
    if sat.parent is not None and hasattr(sat.parent, 'crafts'):
        sat.parent.crafts.pop(sat.name, None)
    sat.parent = target
    target.crafts[sat.name] = sat


def soi_exit(sat, escape_time):
    """
    Transform a craft's orbit and parent relationship upon leaving parent's SOI.
    :param sat: Craft or Satellite
    :param escape_time: time [days] of exit
    """
    exited_body = sat.parent
    outer_parent = exited_body.parent
    if outer_parent is None:
        return

    old_state = sat.orbit.get_state(escape_time)
    parent_state = exited_body.orbit.get_state(escape_time)
    new_state = parent_state + old_state
    new_state.mu = outer_parent.mu

    sat.orbit = new_state.get_orbit()

    # Re-link parents
    if hasattr(exited_body, 'crafts'):
        exited_body.crafts.pop(sat.name, None)
    sat.parent = outer_parent
    outer_parent.crafts[sat.name] = sat
