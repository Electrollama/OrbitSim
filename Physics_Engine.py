"""
Physics Engine for OrbitSim.

Manages the simulation clock, chronological event queue across all crafts,
burn maneuvers, and automated SOI boundary transition scheduling.
"""

from EventQueue import EventQueue, Event, SoiEnter, SoiExit
from encounters import find_encounter, find_escape, soi_enter, soi_exit


class Engine:
    """
    Simulation handler of craft events in a multi-body game system.
    """

    def __init__(self, system, init_time=0.0):
        """
        Initialize an Engine with a system tree (e.g. Kerbol or Kerbin).
        :param system: the overall root parent Body
        :param init_time: starting simulation timestamp [days]
        """
        self.system = system
        self.time = init_time
        self.events = EventQueue()

        # Register all crafts across the entire body hierarchy
        self._register_crafts(self.system, init_time)

    def _register_crafts(self, body, time):
        if hasattr(body, 'crafts'):
            for craft_name, craft in list(body.crafts.items()):
                new_event = check_orbit(craft, time)
                if new_event is not None:
                    self.events.insert(new_event)

        if hasattr(body, 'bodies'):
            for child_name, child_body in list(body.bodies.items()):
                self._register_crafts(child_body, time)

    def elapse(self, t_end):
        """
        Run Events chronologically from the queue until t_end.
        :param t_end: target simulation timestamp [days]
        """
        while not self.events.is_empty():
            next_t = self.events.next_time()
            if next_t > t_end:
                break

            curr = self.events.pop()
            continuation = curr._do()
            if continuation is not None:
                self.events.insert(continuation)

        self.time = t_end


def check_orbit(craft, time):
    """
    Search for future changes of orbit (SOIEnter, SOIExit), then invoke craft program.
    :param craft: a Craft object
    :param time: current simulation time [days]
    :return: the next Event to queue
    """
    if craft.parent is None:
        return None

    new_event = None
    targets = craft.parent.bodies
    next_enc = None

    # Search for earliest SOI encounter across all sibling bodies
    for target_name in targets.keys():
        target = targets[target_name]
        t_enc = find_encounter(craft, time, target.name)
        if t_enc is not None:
            if next_enc is None or t_enc < next_enc[0]:
                next_enc = (t_enc, target)

    if next_enc is not None:
        # Soonest encounter found
        t_enc, target = next_enc
        new_event = SoiEnter(craft, t_enc, target=target)
    else:
        # Check for escape from current parent's SOI
        t_esc = find_escape(craft, time)
        if t_esc is not None:
            new_event = SoiExit(craft, t_esc)
        else:
            # No upcoming SOI changes: schedule regular orbit check for next period
            if not craft.orbit.is_hyp:
                t_next = time + craft.orbit.period()
            else:
                t_next = time + 10.0
            new_event = Event(craft, t_next, check_orbit, (craft, t_next))

    # Evaluate craft's programmed mission plan if available
    if hasattr(craft, 'program') and callable(craft.program):
        prog_event = craft.program(time, new_event.time, getattr(new_event, 'event_type', None))
        if prog_event is not None:
            if prog_event.time <= new_event.time:
                new_event = prog_event
            else:
                msg = 'Programmed event time of {} is too late.'
                raise ValueError(msg.format(craft))

    return new_event


def burn(craft, time, dv, dir_angle=0.0):
    """
    Execute a craft burn maneuver and plan consequences for the new orbit.
    :param craft: Craft object
    :param time: burn execution time [days]
    :param dv: delta-v magnitude [Mm/day]
    :param dir_angle: burn direction angle [rad]
    :return: next Event for the craft
    """
    craft.burn([dv, dir_angle], time)
    return check_orbit(craft, time)
