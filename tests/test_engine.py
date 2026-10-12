"""Integration tests for Physics Engine and automated EventQueue lifecycle."""

import pytest
from Satellite import Craft
from Orbit import mk_Orbit
from Physics_Engine import Engine, check_orbit, burn
from EventQueue import SoiEnter, SoiExit


def test_engine_initialization_and_queue(Kerbin):
    # Add a craft to Kerbin
    craft_orbit = {'ecc': 0.05, 'peri': 1.0, 't0': 0.0, 'lop': 0.0}
    craft = Craft({'name': 'Sputnik', 'orbit': craft_orbit, 'color': 'white'}, parent=Kerbin)

    engine = Engine(Kerbin, init_time=0.0)
    assert not engine.events.is_empty()
    assert engine.events.next_time() > 0.0


def test_automated_soi_enter_and_exit_lifecycle(Kerbin):
    Kerbin.crafts.clear()
    mun = Kerbin.bodies['Mun']
    mun.crafts.clear()

    # Design transfer orbit that encounters Mun
    semi = (0.68 + 12.0) / 2.0
    ecc = (12.0 - 0.68) / (12.0 + 0.68)
    peri = 0.68
    t_target = 1.8
    _, mun_theta = mun.orbit.position(t_target)
    lop_deg = ((mun_theta - 3.1415926535) * 180.0 / 3.1415926535) % 360.0

    temp_orb = mk_Orbit({'ecc': ecc, 'peri': peri, 't0': 0.0, 'lop': lop_deg, 'mu': Kerbin.mu})
    t0_val = t_target - temp_orb.period() / 2.0

    craft_orbit = {'ecc': ecc, 'peri': peri, 't0': t0_val, 'lop': lop_deg}
    craft = Craft({'name': 'MunFlyby', 'orbit': craft_orbit, 'color': 'yellow'}, parent=Kerbin)

    # Initial check_orbit should schedule SoiEnter
    ev = check_orbit(craft, 0.0)
    assert isinstance(ev, SoiEnter)
    t_enter = ev.time
    assert t_enter > 0.0
    assert ev.target is mun

    # Initialize Engine
    engine = Engine(Kerbin, init_time=0.0)

    # Elapse up to slightly after t_enter
    engine.elapse(t_enter + 0.01)

    # Craft should now be parented by Mun!
    assert craft.parent is mun
    assert craft.name in mun.crafts

    # And the queue should automatically contain the upcoming SoiExit event!
    assert not engine.events.is_empty()
    next_ev = engine.events.pop()
    assert isinstance(next_ev, SoiExit)
    assert next_ev.time > t_enter

    # Put it back and elapse past exit time
    engine.events.insert(next_ev)
    t_exit = next_ev.time
    engine.elapse(t_exit + 0.1)

    # Craft should now have exited Mun and returned to Kerbin frame!
    assert craft.parent is Kerbin
    assert craft.name in Kerbin.crafts


def test_craft_burn_maneuver(Kerbin):
    craft_orbit = {'ecc': 0.0, 'peri': 1.0, 't0': 0.0, 'lop': 0.0}
    craft = Craft({'name': 'BurnTest', 'orbit': craft_orbit, 'color': 'cyan'}, parent=Kerbin)
    craft.adjust_res('delta-v', 500.0, 1000.0)

    old_speed = craft.orbit.get_state(0.1).speed
    # Apply prograde burn
    next_ev = burn(craft, 0.1, 50.0, 0.0)
    new_speed = craft.orbit.get_state(0.1).speed

    assert new_speed > old_speed
    assert next_ev is not None
