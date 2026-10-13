"""
Physics Engine for OrbitSim.

Responsible for simulation time orchestration, stepping events from EventQueue,
and executing burns.
"""

__test__ = False

from EventQueue import EventQueue, Event, SoiEnter, SoiExit, check_orbit


class Engine:
    """
    Discrete-event simulation engine for OrbitSim.
    Advances simulation time, manages crafts and celestial systems,
    and executes events in chronological order.
    """

    def __init__(self, system, queue=None):
        self.system = system
        self.time = 0.0
        self.queue = queue if queue is not None else EventQueue()
        self.crafts = {}

    def add_craft(self, craft, time=None):
        """Add craft to the engine and predict its initial event."""
        if time is None:
            time = self.time
        self.crafts[craft.name] = craft
        event = check_orbit(craft, time)
        if event is not None:
            self.queue.insert(event)

    def next_time(self):
        """Next scheduled event time, or None if queue is empty."""
        if self.queue.is_empty():
            return None
        return self.queue.next_time()

    def step(self):
        """
        Execute the next earliest event in the queue.
        Returns the executed event, or None if queue is empty.
        """
        if self.queue.is_empty():
            return None
        event = self.queue.pop()
        self.time = event.time
        next_event = event._do()
        if next_event is not None:
            self.queue.insert(next_event)
        return event

    def advance_to(self, target_time):
        """
        Advance simulation time up to target_time, processing all intervening events.
        """
        while not self.queue.is_empty() and self.queue.next_time() <= target_time:
            self.step()
        self.time = max(self.time, target_time)

    def run_until(self, target_time):
        """Alias for advance_to."""
        self.advance_to(target_time)


def burn(craft, dv_vect, time, queue=None):
    """
    Perform an instantaneous maneuver burn on a craft at `time`.
    Updates the craft's orbit and returns the next predicted event (if any).
    If a queue is provided, inserts the new event into the queue.
    :param craft: Craft object
    :param dv_vect: [dv_magnitude, dv_angle] in Mm/day
    :param time: time [days] of burn
    :param queue: optional EventQueue
    :return: next predicted Event, or None
    """
    craft.burn(dv_vect, time)
    event = check_orbit(craft, time)
    if queue is not None and event is not None:
        queue.insert(event)
    return event


__all__ = ['Engine', 'check_orbit', 'burn']
