import pytest

from Satellite import Craft
from diagrams import plot_system
from load import load_system

from math import pi
from numpy import arange
from os.path import dirname, abspath

project_dir = dirname(dirname(abspath(__file__)))
System = load_system(project_dir + '/System Tables/Bodies_KSP')

@pytest.fixture()
def Kerbol():
    return System

@pytest.fixture()
def Kerbin():
    return System.bodies['Kerbin']

@pytest.fixture()
def craft(Kerbin):
    # circular orbit
    orbit_dict = {'ecc': 0.015, 'peri': .080 + home.size, 't0': -0.1, 'lop': 0.0}
    return Craft({'name': 'KerbSat',
                  'color': 'red',
                  'orbit': orbit_dict},
                  parent=home)
