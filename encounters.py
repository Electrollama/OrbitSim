"""
Geometry for encounter search: positions and separation in a shared parent frame.

Closest-approach search (later phases) should use ca_search_end() so the time
window stays bounded (one craft orbit and/or synodic period, never unbounded).
"""

from math import cos, sin, sqrt

from Orbit import polar_to_xy

def find_encounter(craft, time, target_name, n_orbits=1):
    """
    Find the time, if it exists, that a satellite enters a target body's sphere
    of influence.
    :param craft: Satellite object of the satellite
    :param target: Body object of the target body; Must not be hyperbolic
    :param time: Starting time for search
    :return: time of encounter, or None
    """
    # test impossibility
    target = craft.parent.bodies[target_name]
    if target.orbit.is_hyp:
        raise ValueError('Cannot complete encounter search for ' +
                         'target with hyperbolic orbit.')
    if craft.orbit.peri > target.orbit.apo() + target.soi:
        return None
    elif not craft.orbit.is_hyp and craft.orbit.apo() < target.orbit.peri - target.soi:
        return None
    # find time of minimum distance to sphere of influence
    objective = lambda t: abs(distance(craft, target, t) - target.soi)
    if craft.orbit.is_hyp:
        end_time = time + n_orbits * target.orbit.period()
    else:
        end_time = time + n_orbits * craft.orbit.period()
    opt_res = minimize(objective, time, method='Bounded',
                       bounds=(time, end_time), tol=1 / (24 * 60))
    if opt_res.success and opt_res.fun < target.soi * 0.01:  # encounter found
        return opt_res.x
    else:
        return None


def find_escape(sat, time):
    """
    Find the time, if it exists, that a satellite leaves its parent's sphere
    of influence.
    :param sat: Satellite object of the satellite
    :param target: Body object of the target body; Must not be hyperbolic
    :param time: Starting time for search
    :return: time of departure, or None
    """
    parent = sat.parent
    # test impossibility
    if not sat.orbit.is_hyp and sat.orbit.apo() < parent.soi:
        return None
    # find time of minimum distance to sphere of influence
    objective = lambda t: abs(distance(sat.orbit, parent.orbit, t) - parent.soi)
    if sat.orbit.is_hyp:
        end_time = time + parent.orbit.period()
    else:
        end_time = time + 2 * sat.orbit.period()
    opt_res = minimize(objective, time,
                       bounds=(time, end_time), tol=1 / (24 * 60))
    if opt_res.success and opt_res.fun < parent.soi * 0.01:  # encounter found
        return opt_res.x
    else:
        return None


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


def _require_orbit(body, label):
    if body.orbit is None:
        raise EncounterGeometryError(
            '{} {!r} has no orbit.'.format(label, body.name))


def distance(craft, target, t):
    xa, ya = position_xy(body_a, t)
    xb, yb = position_xy(body_b, t)
    return sqrt((xa - xb) ** 2 + (ya - yb) ** 2)
    return add_polar(r2, r1, a2, a1, subtr=True)[0]


def soi_enter(sat, time, target):
    """
    Find a satellite's new orbit with respect to the encountered target
    :param sat: Satellite object of the satellite
    :param target: Body object of the target body
    :param time: found time that the satellite enters the target's SOI
    :return: [None; mutable]
    """
    old_state = sat.orbit.get_state(time)
    target_state = target.orbit.get_state(time)
    new_state = old_state - target_state
    new_state.mu = target.mu
    sat.orbit = new_state.get_orbit()
    sat.parent.bodies.pop(sat.name, None)  # remove from old parent's satellites
    sat.parent = target  # switch parents
    target.bodies[sat.name] = sat  # add to new parent's satellites


def soi_exit(sat, escape_time):
    """
    Find a satellite's new orbit with respect to the escaped parent's parent
    :param sat: Satellite object of the satellite
    :param encounter_time: found time that the satellite enters the target's SOI
    :return: [None; mutable]
    """
    if sat.parent.parent is None:  # orbiting sun
        if sat.orbit.is_hyp:
            sat.parent.soi *= 2
        else:
            sat.parent.soi = sat.orbit.apo() * 1.01
    old_state = sat.orbit.get_state(escape_time)
    parent_state = sat.parent.orbit.get_state(escape_time)
    new_state = parent_state + old_state  # also changes mu to that of parent.parent
    sat.orbit = new_state.get_orbit()
    sat.parent.bodies.pop(sat.name, None)  # remove from old parent's satellites
    sat.parent = sat.parent.parent  # switch parents
    sat.parent.bodies[sat.name] = sat  # add to new parent's satellites


class EncounterGeometryError(ValueError):
    """Invalid inputs for encounter geometry."""


def common_parent(body_a, body_b):
    """True if both satellites share the same parent Body."""
    return body_a.parent is not None and body_a.parent is body_b.parent


def assert_common_parent(body_a, body_b):
    if not common_parent(body_a, body_b):
        names = body_a.name, body_b.name
        raise EncounterGeometryError(
            'Bodies must share a parent for separation: {} and {}.'.format(*names))


def separation_fn(body_a, body_b):
    """
    :return: callable f(t) -> separation distance [Mm]
    """
    assert_common_parent(body_a, body_b)

    def distance_at(t):
        return distance(body_a, body_b, t)

    return distance_at


def synodic_period(orbit_a, orbit_b):
    """
    Synodic period for two bound Keplerian orbits about the same focus.
    :return: period [days], or None if either orbit is hyperbolic/parabolic
    """
    if orbit_a.is_hyp or orbit_b.is_hyp:
        return None
    ta = orbit_a.period()
    tb = orbit_b.period()
    diff = abs(1.0 / ta - 1.0 / tb)
    if diff < 1e-15:
        return min(ta, tb)
    return 1.0 / diff


def ca_search_end(t0, craft, target, craft_orbits=1.0, max_days=None):
    """
    Upper time bound for a closest-approach search starting at t0.

    Uses the shorter of (craft_orbits * craft period) and the synodic period
    with the target, so the window does not grow without limit when many years
    of simulation time have passed—only the interval (t0, t_end) matters.

    Separation d(t) is not a single sinusoid, but for two Keplerian paths about
    the same parent it repeats on scales set by each period and the synodic
    period; searching min(craft orbit, synodic) catches the next meaningful
    alignment without scanning arbitrarily far ahead.

    :param t0: start time [days]
    :param craft: Craft (or any Satellite) whose path is propagated
    :param target: other Body under the same parent
    :param craft_orbits: width in multiples of the craft's orbital period
    :param max_days: optional hard cap on window length [days]
    :return: t_end [days]
    """
    assert_common_parent(craft, target)
    _require_orbit(craft, 'craft')
    _require_orbit(target, 'target')
    if craft_orbits <= 0:
        raise EncounterGeometryError('craft_orbits must be positive.')
    orbit_c = craft.orbit
    if orbit_c.is_hyp:
        raise EncounterGeometryError(
            'Hyperbolic craft leg needs a parent-SOI exit time (not implemented).')
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










