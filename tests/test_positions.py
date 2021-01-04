import pytest
from math import pi

# list of initial planet true anomalies, +- 5 degrees, relative to Kerbin
PlanetPos = {"Moho": 84.7,
              "Eve": 15.0,
             "Kerbin": 0.0,
             "Duna": 135.5,
             "Dres": 10.2,
             "Jool": -122.0,
             "Eeloo": -50.5}

precision = 0.02 #radians

# true anomalies of bodies around Kerbin, relative to Kerbol
MoonPos = {"Mun": 101,
           "Minmus": 110}

def test_Kerbin_period(Kerbin):
    assert abs(Kerbin.orbit.period() - 426.09) < 1

def test_Mun_period(Kerbin):
    assert abs(Kerbin.bodies['Mun'].orbit.period() - 6.417) < precision

def test_Kerbin_MA(Kerbin):
    actual =  Kerbin.orbit.mean_anom(0.0)%(2*pi)
    expected = 3.14
    delta = abs(actual - expected)
    assert abs(delta) < precision or abs(delta - 2*pi) < precision

def test_planet_pos(Kerbol):
    Kerbin_pos = Kerbol.bodies['Kerbin'].orbit.position(0.0)[1]
    for planet_name in Kerbol.bodies.keys():
        planet = Kerbol.bodies[planet_name]
        actual = (planet.orbit.position(0.0)[1] - Kerbin_pos)%(2*pi)
        expected = (PlanetPos[planet_name]*pi/180)%(2*pi)
        assert (abs(actual - expected) < precision
                or abs(actual - expected - (2*pi)) < precision
                or abs(actual - expected + (2*pi)) < precision)
        print(planet_name)
