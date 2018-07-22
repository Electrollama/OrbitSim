from load import load_system
from Satellite import Craft
from Orbit import Orbit, State
from orbit_plotter import plot_body, plot_orbit, plot_show, polar_to_xy, plt
from numpy import arange
from math import pi

def orbit_family(Start, param, step, n=5, init_val=None):
    """
    Plot orbits while varying a single parameter.
    :param Start: Initial orbit
    :param param: parameter to vary
    :param step: step size between test values
    :param n: number of steps above and below
    """
    if init_val is None:
        init_val = Start.get_element(param)
    minval = init_val - step * n
    maxval = init_val + step * (n+1)
    vals = arange(minval, maxval, step)
    for x in vals:
        # get the orbit instance
        if isinstance(Start, State):
            St = Start.set_element(param, x)
            Orb = St.get_orbit()
        elif isinstance(Start, Orbit):
            Orb = Start.set_element(param, x)
        else:
            raise TypeError
        # color according to comparison with the initial value
        if abs(x - init_val) < step*0.9: #x == init_val
            c = 'yellow'
        elif x < init_val:
            c = 'red'
        else: #x > init_val:
            c = 'green'
        plot_orbit(Orb, color=c, ls='-', marker='s')
    plot_show()

def test_burn(t_burn, v_burn):
    """
    Test a burn to a new orbit
    :param t_burn: time to perform the burn
    :param v_burn: vector (polar) representing the burn
    """
    craft = Craft({'name': 'test',
                   'size': 10,
                   'orbit': {'mu': 547e6,
                             'peri': 20e3,
                             'ecc': 0.234,
                             't0': 0.0,
                             'lop': 0.0},
                   'color': 'white'})
    plot_orbit(craft.orbit, color='#f1ffb7')
    craft.adjust_res('delta-v', 1000.0, 1500.0)
    x_burn, y_burn = polar_to_xy(craft.orbit.position(t_burn))
    craft.burn([v_burn, 0.0], t_burn)
    plot_orbit(craft.orbit, color='#ffd1d1', ls='-')
    plt.plot([x_burn,], [y_burn,], 'ro')
    plot_show()

def planet_chart(System, start=0):
    """
    Make a toy-scale chart of the system.
    :param System: The system's parent Body object
    :param start: Start scaling according to the Nth planet.
            (useful for outer planets)
    """
    Sats = System.satellites
    planets = [Sats[p_name] for p_name in Sats.keys()]
    semis = [(p.orbit.semi(), p) for p in planets]
    semis_sorted = list(semis)
    semis_sorted.sort()
    semis_sorted = semis_sorted[start:]
    scaling = 0.0
    for i in range(len(semis_sorted) - 1):
        curr, next = semis_sorted[i], semis_sorted[i+1]
        d1, p1 = curr
        d2, p2 = next
        curr_scale = (p1.size + p2.size) * 1.5 / (d2-d1)
        scaling = max(scaling, curr_scale)
    for p in planets:
        plot_body(p, p.orbit.semi()*scaling, 0.0)
    plot_body(System, -System.size, 0.0)
    plot_show()

def test_param(param):
    """
    Test the orbit-family function
    :param param: name of the parameter to change
    :return:
    """
    System = load_system('Celestial_Bodies_KSP')
    Test_Orbit = System.satellites['Kerbin'].orbit
    Test_State = Test_Orbit.get_state(0.0)
    print(Test_State.Elements())
    print(Test_Orbit.Elements())
    plot_body(System)
    orbit_family(Test_State, 'param', 10)

System = load_system('Celestial_Bodies_RL')
planet_chart(System)