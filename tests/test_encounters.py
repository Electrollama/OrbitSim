"""Unit tests for encounter detection and SOI boundary crossings."""

import pytest
from Satellite import Craft
from encounters import find_encounter, find_escape, separation, soi_enter, soi_exit


def test_encounter_gating_out_of_range(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Craft in tiny orbit well below Mun SOI (Mun peri = 12.0, Mun SOI ~ 2.42 Mm)
    tiny_orbit = {'ecc': 0.0, 'peri': 1.0, 't0': 0.0, 'lop': 0.0}
    craft = Craft({'name': 'Probe', 'orbit': tiny_orbit, 'color': 'white'}, parent=Kerbin)

    t_enc = find_encounter(craft, 0.0, 'Mun')
    assert t_enc is None


def test_encounter_mun_transfer(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Setup Hohmann transfer from LKO to Mun
    # Mun is at r = 12.0. At t = 0, Mun is at some position.
    # Craft periapsis = 0.68 Mm, apoapsis = 12.0 Mm
    # semi = (0.68 + 12.0)/2 = 6.34 Mm. ecc = (12.0 - 0.68)/(12.0 + 0.68) ~ 0.8927
    semi = (0.68 + 12.0) / 2.0
    ecc = (12.0 - 0.68) / (12.0 + 0.68)
    peri = 0.68

    # We want craft to hit Mun at apoapsis:
    # At t_apo = craft.period()/2, craft is at r = 12.0, true anomaly = pi.
    # So lop can be chosen so that craft meets Mun at Mun's position.
    # Let's find Mun's position at t_target:
    t_target = 1.8
    _, mun_theta = mun.orbit.position(t_target)
    # Craft reaches apoapsis at t_target with theta = mun_theta
    # For craft, apoapsis position is lop + pi. So lop = mun_theta - pi.
    lop_deg = ((mun_theta - 3.1415926535) * 180.0 / 3.1415926535) % 360.0
    # craft.t0 is set so apoapsis is at t_target: t0 = t_target - period/2
    from Orbit import mk_Orbit
    temp_orb = mk_Orbit({'ecc': ecc, 'peri': peri, 't0': 0.0, 'lop': lop_deg, 'mu': Kerbin.mu})
    t0_val = t_target - temp_orb.period() / 2.0

    craft_orbit = {'ecc': ecc, 'peri': peri, 't0': t0_val, 'lop': lop_deg}
    craft = Craft({'name': 'MunTransfer', 'orbit': craft_orbit, 'color': 'white'}, parent=Kerbin)

    t_enc = find_encounter(craft, 0.0, 'Mun', n_orbits=2.0)
    assert t_enc is not None
    # Verify that at t_enc, separation is approximately Mun SOI
    sep = separation(craft, mun, t_enc)
    assert sep == pytest.approx(mun.soi, rel=1e-3)


def test_escape_hyperbolic_orbit(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Craft in hyperbolic orbit around Mun moving away from periapsis
    hyp_orbit = {
        'ecc': 1.5,
        'peri': 0.3,
        't0': 0.0,
        'lop': 0.0
    }
    craft = Craft({'name': 'Escaper', 'orbit': hyp_orbit, 'color': 'white'}, parent=mun)

    t_esc = find_escape(craft, 0.0)
    assert t_esc is not None
    assert t_esc > 0.0
    # Verify distance at escape time equals Mun SOI
    dist_at_esc = craft.orbit.hdist(t_esc)
    assert dist_at_esc == pytest.approx(mun.soi, rel=1e-3)


def test_soi_enter_and_exit_frame_switch(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Place craft at Mun SOI boundary
    t_event = 1.0
    # Create craft orbiting Kerbin
    craft_orbit = {'ecc': 0.5, 'peri': 5.0, 't0': 0.0, 'lop': 0.0}
    craft = Craft({'name': 'Traveler', 'orbit': craft_orbit, 'color': 'white'}, parent=Kerbin)
    assert craft.name in Kerbin.crafts

    # Perform SOI Enter
    soi_enter(craft, t_event, mun)
    assert craft.parent is mun
    assert craft.name in mun.crafts
    assert craft.name not in Kerbin.crafts

    # Perform SOI Exit
    soi_exit(craft, t_event + 0.5)
    assert craft.parent is Kerbin
    assert craft.name in Kerbin.crafts
    assert craft.name not in mun.crafts
