from scipy.optimize import minimize_scalar as minimize
from Orbit import add_polar
from encounters import *
from EventQueue import *

"""
Controls events in a game, including crafts changing sphere of influence

To do:
* add method docstrings
? allow user to create a check_orbit for multiple orbits
* implement array heap
* distinguish between bodies and crafts (check type or rename satellite dict?)
    - bodies can be a parent
    - crafts can have hyperbolic orbits
    - rename function inputs to "craft" instead of "sat"
* add under-scores to private methods and variables
* test a print_system event
* test an orbiting craft
* test an encounter
* implement burns (change orbits, delete future events?)
* test a burn into an encounter
* test an insertion burn
* test an ejection burn
* test multiple encounters and burns
* video showing example flight
* figure out folder organization
* move on to the user interface
Alpha complete
* user tools for orbit planning
* user tools for winging it
* move on to minor features
    - limited fuel
    - rendezvous?
    - easy commands?
        + ex. enter into this orbit
        + burn until apoapsis is some value
        + relative angle for transfer
    
"""

#### Event Handling ####
class Engine: #
    """
    A handler of craft events in a game state
    """
    def __init__(self, system, init_time=0.0):
        """
        Initialize an Engine of a new or loaded game state
        :param system: the overall parent Body, ex. the Sun
        """
        self.system = system
        self.time = init_time
        # add events
        self.events = EventQueue()
        for name in system.crafts.keys():
            sat = system.crafts[name]
            new_event = check_orbit(sat, init_time)
            self.events.insert(new_event)

    def elapse(self, t_end):
        """
        Run Events chronologically from the queue
        :param t_end: time to stop running events
        """
        while self.events.next_time() < t_end:
            curr = self.events.pop()
            if curr.func_type == 'check_orbit':
                new = check_orbit(curr.craft)
            else:
                curr.func(*curr.args) # complete event
                new = check_orbit(curr.craft)
            self.events.insert(new) # queue a new event for the same craft
        self.time = t_end

"""
Typical events that the engine handles
"""

def check_encounter(craft1, craft2, time_start):
    """
    Check for future changes of orbit, then run the craft program
    :param time: current time
    :param craft: a Craft object
    :return: the next event to queue
    """
    #find close approaches
    #find SOI changes
    return #events

def check_orbit(craft, time):
    """
    Check for future changes of orbit, then run the craft program
    :param time: current time
    :param craft: a Craft object
    :return: the next event to queue
    """
    # default is the next check_orbit
    print('Checking' ,craft.name)
    new_event = None
    # check for encounters
    Targets = craft.parent.bodies
    target_list = [Targets[name] for name in Targets.keys()]
    next_enc = None
    for target in target_list:
        t_enc = find_encounter(craft, time, target.name)
        if (t_enc is not None and # encounter found
                (next_enc is None or t_enc < next_enc[0])): # first one or the soonest one
            next_enc = (t_enc, target)
    if next_enc is None:
        # check for escape (an elliptical orbit can escape)
        t_esc = find_escape(craft, t_enc)
        if t_esc is not None:
            new_event = Event(craft, t_esc, soi_exit, ())  # escape
    else: # encounter found (soonest one)
        t_enc, target = next_enc
        new_event = Event(craft, t_enc, soi_enter, (target,))
    # run program
    if new_event is None:
        new_event = Event(craft, time, None) #next check_orbit
    prog_event = craft.program(time, new_event.time, new_event.event_type)
    if prog_event is None:
        pass
    elif prog_event.time < new_event.time:
        new_event = prog_event
    else:  # programmed event is too late
        msg = 'Programmed event time of {} is too late.'
        raise ValueError(msg.format(craft))
    return new_event

def burn(craft, time, dv, dir):
    """
    Execute a craft burn and plan the consequences of the new orbit
    :param craft:
    :param time:
    :param dv:
    :param dir:
    :return:
    """
    # if there's an encounter before then, add the encounter and return an error
    # else, add the burn event
    craft.burn(time, dv, dir)
    #new orbit
    check_orbit(craft, time)
    return

#### Changing Sphere of Influence ####


