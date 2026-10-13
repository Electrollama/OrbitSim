"""Unit tests for Physics_Engine and EventQueue integration."""

try:
    import pytest
except ImportError:
    class _MockPytest:
        class Approx:
            def __init__(self, expected, rel=1e-3, abs=1e-6):
                self.expected = expected
                self.rel = rel
                self.abs = abs
            def __eq__(self, actual):
                diff = abs(self.expected - actual)
                tol = max(self.abs, abs(self.expected) * self.rel)
                return diff <= tol
        approx = Approx
    pytest = _MockPytest()

from Satellite import Craft
from Physics_Engine import Engine, check_orbit, burn
from EventQueue import EventQueue, SoiEnter, SoiExit


def test_engine_initialization(Kerbin):
    engine = Engine(Kerbin)
    assert engine.system is Kerbin
    assert engine.time == 0.0
    assert engine.queue.is_empty()
    assert engine.next_time() is None


def test_check_orbit_escape_detection(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Craft on hyperbolic escape from Mun
    hyp_orbit = {
        'ecc': 1.5,
        'peri': 0.3,
        't0': 0.0,
        'lop': 0.0
    }
    craft = Craft({'name': 'Escaper', 'orbit': hyp_orbit, 'color': 'white'}, parent=mun)
    event = check_orbit(craft, 0.0)
    assert event is not None
    assert isinstance(event, SoiExit)
    assert event.time > 0.0


def test_check_orbit_encounter_detection(Kerbin):
    mun = Kerbin.bodies['Mun']
    # Setup transfer orbit that encounters Mun
    semi = (0.68 + 12.0) / 2.0
    ecc = (12.0 - 0.68) / (12.0 + 0.68)
    peri = 0.68
    t_target = 1.8
    _, mun_theta = mun.orbit.position(t_target)
    lop_deg = ((mun_theta - 3.1415926535) * 180.0 / 3.1415926535) % 360.0

    from Orbit import mk_Orbit
    temp_orb = mk_Orbit({'ecc': ecc, 'peri': peri, 't0': 0.0, 'lop': lop_deg, 'mu': Kerbin.mu})
    t0_val = t_target - temp_orb.period() / 2.0

    craft_orbit = {'ecc': ecc, 'peri': peri, 't0': t0_val, 'lop': lop_deg}
    craft = Craft({'name': 'MunTransfer', 'orbit': craft_orbit, 'color': 'white'}, parent=Kerbin)

    event = check_orbit(craft, 0.0, n_orbits=2.0)
    assert event is not None
    assert isinstance(event, SoiEnter)
    assert event.target.name == 'Mun'


def test_burn_alters_orbit(Kerbin):
    orbit_dict = {'ecc': 0.01, 'peri': 1.0, 't0': 0.0, 'lop': 0.0}
    craft = Craft({'name': 'Probe', 'orbit': orbit_dict, 'color': 'white'}, parent=Kerbin)
    craft.adjust_res('delta-v', 500.0, 500.0)

    initial_peri = craft.orbit.peri
    burn(craft, [50.0, 0.0], 0.0)
    # Orbit should have changed
    assert craft.orbit.peri != initial_peri or craft.orbit.ecc != 0.01


def test_engine_step_and_advance(Kerbin):
    mun = Kerbin.bodies['Mun']
    hyp_orbit = {
        'ecc': 1.5,
        'peri': 0.3,
        't0': 0.0,
        'lop': 0.0
    }
    craft = Craft({'name': 'Escaper', 'orbit': hyp_orbit, 'color': 'white'}, parent=mun)
    engine = Engine(mun)
    engine.add_craft(craft, 0.0)

    assert not engine.queue.is_empty()
    t_next = engine.next_time()
    assert t_next is not None

    event = engine.step()
    assert event is not None
    assert engine.time == pytest.approx(t_next)
    # Craft should now have switched to Kerbin as parent
    assert craft.parent is Kerbin
