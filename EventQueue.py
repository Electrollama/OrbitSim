import bisect
import heapq
from itertools import count

class EventQueue:
    """Global chronological queue across all crafts. AI assisted."""

    def __init__(self):
        self._heap = []            # entries: (time, seq, event)
        self._queues = {}          # craft -> CraftEventQueue
        self._seq = count()        # unique tiebreak; events themselves never get compared

    def is_empty(self):
        return not any(self._queues.values())

    def insert(self, event):
        """Queue an event for event.craft."""
        cq = self._queues.get(event.craft)
        if cq is None:
            cq = self._queues[event.craft] = CraftEventQueue()
        cq.insert(event)
        if cq.see_next() is event:            # new earliest for this craft
            heapq.heappush(self._heap, (event.time, next(self._seq), event))
            # the craft's previous head entry is now stale; next() skips it

    def next(self):
        """Remove and return the earliest event across all crafts."""
        heap = self._heap
        while heap:
            _, _, event = heap[0]
            cq = self._queues[event.craft]
            if not cq or cq.see_next() is not event:   # stale entry
                heapq.heappop(heap)
                continue
            cq.pop()
            if cq:                                      # craft's next event takes its place
                head = cq.see_next()
                heapq.heapreplace(heap, (head.time, next(self._seq), head))
            else:
                heapq.heappop(heap)
            return event
        raise IndexError('Event queue is empty.')

    def next_time(self):
        """Time of the earliest event without removing it."""
        heap = self._heap
        while heap:
            time, _, event = heap[0]
            cq = self._queues[event.craft]
            if cq and cq.see_next() is event:
                return time
            heapq.heappop(heap)
        raise IndexError('Event queue is empty.')


def _neg_time(e):
    return -e.time

class CraftEventQueue:
    """Events for one craft, sorted so the earliest is at the back (cheap pop)."""

    def __init__(self):
        self.queue = []

    def insert(self, event):
        bisect.insort(self.queue, event, key=_neg_time)   # needs Python 3.10+

    def pop(self):
        return self.queue.pop()

    def see_next(self):
        return self.queue[-1]

    def __bool__(self):
        return bool(self.queue)

class Event:
    """
    An operation related to a craft at a future time; ex. change in orbit, alarm
    variables: event_type [str], time [float], func[function],
        args[func arguments], craft [Craft]
    """
    def __init__(self, craft, time, args=()):
        self.craft = craft
        self.args = args
        self.time = time

    def __repr__(self):
        return "<{} of {} at {}d>".format(self.event_type,
                                       self.craft.name,
                                       round(self.time))

    def __lt__(self, other):
        return self.time < other.time

    def __eq__(self, other):
        return self.time == other.time

    def __gt__(self, other):
        return self.time > other.time

    def _do(self):
        """
        Action to complete at self.time
        :return: new event if relevant
        """
        return

class CloseApproach(Event):
    """
    An event representing a found close approach.
    Draw as a line connecting the craft position and target position at the closest approach time.
    If it's an approach to a body, draw as a 2-color line showing the SOI crossing.
    """
    def __init__(self, craft, time, func, t_search_end, dist, target):
        Event.__init__(self, craft, time, func)
        self.t_search_end = t_search_end
        self.dist = dist
        self.target = target

    def _do(self):
        """
        Find the next close approach
        :return:
        """
        return

class MOID(Event):
    """
    Minimum Orbital Intersection Distance, the geometric close approaches. Does not change per orbital period.
    Can be a single event or a pair of events (ascending or descending).
    Draw as dotted lines connecting the crossing points (or geometric minimum) and the target at that position.
    """
    def __init__(self, craft, time, target, moid_type):
        """
        A geometric close approach. ToDo: determine if hyperbolic orbits need their own implementation
        :param target: "min" (non-intersecting), "asc" (crossing ascending), "desc" (crossing descending)
        :param type:
        """
        Event.__init__(self, craft, time)
        self.target = target
        self.moid_type = moid_type

    def _do(self):
        if self.craft.orbit.is_hyp:
            return None
        else:
            self.time += self.craft.period() #recycle


class SoiEnter(Event):
    """
    An event representing a found encounter start
    """
    def __init__(self, craft, time, func, target, enter=True):
        Event.__init__(self, craft, time, func)
        self.target = target

    def new_orbit(self): #run when the event occurs
        # get state relative to target at the time
        rel_state = self.craft.orbit.get_state(self.time) - self.target.orbit.get_state(self.time)
        # convert to orbit around the target
        rel_state.mu = self.target.mu
        orbit2 = rel_state.get_orbit()
        self.craft.orbit = orbit2
        return

    def _do(self):
        return None #delete when done


class SoiExit(Event):
    """
    An event representing a found encounter end
    """
    def __init__(self, craft, time, func):
        Event.__init__(self, craft, time, func)

    def new_orbit(self):
        # get state relative to target at the time
        rel_state = self.craft.orbit.get_state(self.time) - self.craft.parent.orbit.get_state(self.time)
        # convert to orbit around the target
        rel_state.mu = self.target.mu
        orbit2 = rel_state.get_orbit()
        self.craft.orbit = orbit2
        return

    def _do(self):
        return None #delete when done

