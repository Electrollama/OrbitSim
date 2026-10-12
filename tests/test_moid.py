"""Unit tests for 2D Geometric MOID calculation."""

import pytest
from math import pi
from Orbit import mk_Orbit
from encounters import geometric_moid_2d


def test_concentric_circles_moid():
    # Orbit 1: r = 10, e = 0
    orb1 = mk_Orbit({'ecc': 0.0, 'peri': 10.0, 't0': 0.0, 'lop': 0.0, 'mu': 1000.0})
    # Orbit 2: r = 15, e = 0
    orb2 = mk_Orbit({'ecc': 0.0, 'peri': 15.0, 't0': 0.0, 'lop': 0.0, 'mu': 1000.0})

    res = geometric_moid_2d(orb1, orb2)
    assert not res['intersects']
    assert res['moid'] == pytest.approx(5.0, rel=1e-5)


def test_intersecting_orbits_moid():
    # Circular orbit at r = 12
    orb1 = mk_Orbit({'ecc': 0.0, 'peri': 12.0, 't0': 0.0, 'lop': 0.0, 'mu': 1000.0})
    # Elliptic orbit with peri = 10 and apo = 14 (semi = 12, e ~ 0.1667)
    # peri / (1 - e) = 12 -> 10 / (1 - e) = 12 -> 1 - e = 10/12 = 5/6 -> e = 1/6
    orb2 = mk_Orbit({'ecc': 1.0 / 6.0, 'peri': 10.0, 't0': 0.0, 'lop': 0.0, 'mu': 1000.0})

    res = geometric_moid_2d(orb1, orb2)
    assert res['intersects']
    assert res['moid'] == pytest.approx(0.0, abs=1e-7)
    assert len(res['thetas']) == 2


def test_separated_ellipses_moid():
    # Inner ellipse: peri = 2.0, apo = 4.0 (semi = 3.0, e = 1/3) with apoapsis along theta=0
    orb1 = mk_Orbit({'ecc': 1.0 / 3.0, 'peri': 2.0, 't0': 0.0, 'lop': 180.0, 'mu': 1000.0})
    # Outer ellipse: peri = 8.0, apo = 12.0 (semi = 10.0, e = 0.2) with periapsis along theta=0
    orb2 = mk_Orbit({'ecc': 0.2, 'peri': 8.0, 't0': 0.0, 'lop': 0.0, 'mu': 1000.0})

    res = geometric_moid_2d(orb1, orb2)
    assert not res['intersects']
    # Closest approach is along theta=0 where orb1 apo is 4.0 and orb2 peri is 8.0 -> dist = 4.0
    assert res['moid'] == pytest.approx(4.0, rel=1e-3)


def test_identical_orbits_moid():
    orb1 = mk_Orbit({'ecc': 0.1, 'peri': 5.0, 't0': 0.0, 'lop': 45.0, 'mu': 1000.0})
    orb2 = mk_Orbit({'ecc': 0.1, 'peri': 5.0, 't0': 0.0, 'lop': 45.0, 'mu': 1000.0})

    res = geometric_moid_2d(orb1, orb2)
    assert res['intersects']
    assert res['moid'] == pytest.approx(0.0, abs=1e-9)
