"""Phase 0: shared-frame geometry and search-window bounds."""

from math import cos, sqrt
import pytest

from Satellite import Craft

from encounters import (
    ca_search_end,
    common_parent,
    position_xy,
    relative_period,
    separation,
    separation_fn,
    EncounterGeometryError,
)


def test_common_parent(Kerbin):
    mun = Kerbin.bodies['Mun']
    minmus = Kerbin.bodies['Minmus']
    assert common_parent(mun, minmus)
    assert not common_parent(Kerbin, mun)


def test_separation_same_place(Kerbin):
    mun = Kerbin.bodies['Mun']
    t = 0.0
    assert separation(mun, mun, t) == pytest.approx(0.0, abs=1e-9)


def test_separation_fn_matches(Kerbin, craft):
    mun = Kerbin.bodies['Mun']
    t = 3.5
    f = separation_fn(craft, mun)
    assert f(t) == pytest.approx(separation(craft, mun, t))


def test_position_xy_consistent_with_polar(Kerbin, craft):
    t = 1.0
    r, theta = craft.orbit.position(t)
    x, y = position_xy(craft, t)
    assert sqrt(x * x + y * y) == pytest.approx(r, rel=1e-9)
    assert x == pytest.approx(r * cos(theta), rel=1e-9)


def test_separation_requires_shared_parent(Kerbol, Kerbin):
    duna = Kerbol.bodies['Duna']
    mun = Kerbin.bodies['Mun']
    with pytest.raises(EncounterGeometryError):
        separation(duna, mun, 0.0)


def test_relative_period_mun_kerbin_orbits(Kerbin, craft):
    mun = Kerbin.bodies['Mun']
    syn = relative_period(Kerbin.orbit, mun.orbit)
    assert syn == pytest.approx(6.51, rel=0.02)


def test_ca_search_end_uses_shorter_synodic(Kerbin, craft):
    mun = Kerbin.bodies['Mun']
    t0 = 1000.0
    t_end = ca_search_end(t0, craft, mun, craft_orbits=1.0)
    craft_period = craft.orbit.period()
    syn = relative_period(craft.orbit, mun.orbit)
    assert t_end - t0 == pytest.approx(min(craft_period, syn), rel=1e-6)


def test_ca_search_end_respects_max_days(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Craft in high orbit with period ~ 13.8 days so orbital and synodic windows exceed max_days
    high_craft = Craft({
        'name': 'HighSat',
        'orbit': {'ecc': 0.0, 'peri': 20.0, 't0': 0.0, 'lop': 0.0}
    }, parent=Kerbin)
    t0 = 0.0
    t_end = ca_search_end(t0, high_craft, mun, craft_orbits=1.0, max_days=2.0)
    assert t_end - t0 == pytest.approx(2.0)
