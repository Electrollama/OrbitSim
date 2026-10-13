try:
    import pytest
except ImportError:
    class _MockPytest:
        @staticmethod
        def fixture(*args, **kwargs):
            def decorator(f):
                return f
            return decorator
    pytest = _MockPytest()

from os.path import dirname, abspath
from Physics.Satellite import Craft
from Load.load import load_system

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
    # Craft in orbit around Kerbin
    orbit_dict = {'ecc': 0.015, 'peri': 0.080 + Kerbin.size, 't0': -0.1, 'lop': 0.0}
    return Craft({
        'name': 'KerbSat',
        'color': 'red',
        'orbit': orbit_dict
    }, parent=Kerbin)
