"""Test runner that executes pytest fixtures and tests without requiring external pytest binary."""

__test__ = False

import math
import sys
import os

# Ensure project root is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tests import conftest
import tests.test_positions as tp
import tests.test_encounters as te
import tests.test_engine as teng


class Approx:
    def __init__(self, expected, rel=1e-3, abs=1e-6):
        self.expected = expected
        self.rel = rel
        self.abs = abs

    def __eq__(self, actual):
        diff = abs(self.expected - actual)
        tol = max(self.abs, abs(self.expected) * self.rel)
        return diff <= tol

    def __repr__(self):
        return f"approx({self.expected})"


# Mock pytest.approx for testing
class MockPytest:
    approx = Approx

tp.pytest = MockPytest
te.pytest = MockPytest
teng.pytest = MockPytest

passed = 0
failed = 0
errors = []

Kerbol_fixture = conftest.Kerbol()
Kerbin_fixture = conftest.Kerbin()

tests_to_run = [
    # test_positions
    ("test_Kerbin_period", lambda: tp.test_Kerbin_period(Kerbin_fixture)),
    ("test_Mun_period", lambda: tp.test_Mun_period(Kerbin_fixture)),
    ("test_Kerbin_MA", lambda: tp.test_Kerbin_MA(Kerbin_fixture)),
    ("test_planet_pos", lambda: tp.test_planet_pos(Kerbol_fixture)),
    # test_encounters
    ("test_encounter_gating_out_of_range", lambda: te.test_encounter_gating_out_of_range(Kerbin_fixture)),
    ("test_encounter_mun_transfer", lambda: te.test_encounter_mun_transfer(Kerbin_fixture)),
    ("test_escape_hyperbolic_orbit", lambda: te.test_escape_hyperbolic_orbit(Kerbin_fixture)),
    ("test_soi_enter_and_exit_frame_switch", lambda: te.test_soi_enter_and_exit_frame_switch(Kerbin_fixture)),
    # test_engine
    ("test_engine_initialization", lambda: teng.test_engine_initialization(Kerbin_fixture)),
    ("test_check_orbit_escape_detection", lambda: teng.test_check_orbit_escape_detection(Kerbin_fixture)),
    ("test_check_orbit_encounter_detection", lambda: teng.test_check_orbit_encounter_detection(Kerbin_fixture)),
    ("test_burn_alters_orbit", lambda: teng.test_burn_alters_orbit(Kerbin_fixture)),
    ("test_engine_step_and_advance", lambda: teng.test_engine_step_and_advance(Kerbin_fixture)),
]

print("Running test suite...")
for name, test_fn in tests_to_run:
    try:
        test_fn()
        passed += 1
        print(f"  PASS: {name}")
    except Exception as e:
        failed += 1
        errors.append((name, e))
        print(f"  FAIL: {name} -> {e}")

print(f"\nResults: {passed} passed, {failed} failed.")
if failed > 0:
    sys.exit(1)
