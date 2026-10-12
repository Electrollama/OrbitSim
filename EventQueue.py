import bisect
import heapq
from itertools import count
from Physics_Engine import check_orbit
from encounters import find_escape


class EventQueue:
    """Global chronological queue across all crafts."""

    def __init__(self):
        self._heap = []            # entries: (time, seq, event)
        self._queues = {}          # craft -> CraftEventQueue
        self._seq = count()        # unique tiebreak

    def is_empty(self):
        return not any(self._queues.values())

    def insert(self, event):
        """Queue an event for event.craft."""
        if event is None:
            return
        cq = self._queues.get(event.craft)
        if cq is None:
            cq = self._queues[event.craft] = CraftEventQueue()
        cq.insert(event)
        if cq.see_next() is event:
            heapq.heappush(self._heap, (event.time, next(self._seq), event))

    def next(self):
        """Remove and return the earliest event across all crafts."""
        heap = self._heap
        while heap:
            _, _, event = heap[0]
            cq = self._queues.get(event.craft)
            if not cq or cq.see_next() is not event:
                heapq.heappop(heap)
                continue
            cq.pop()
            if cq:
                head = cq.see_next()
                heapq.heapreplace(heap, (head.time, next(self._seq), head))
            else:
                heapq.heappop(heap)
            return event
        raise IndexError('Event queue is empty.')

    def pop(self):
        """Alias for next()."""
        return self.next()

    def next_time(self):
        """Time of the earliest event without removing it."""
        heap = self._heap
        while heap:
            time, _, event = heap[0]
            cq = self._queues.get(event.craft)
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
        bisect.insort(self.queue, event, key=_neg_time)

    def pop(self):
        return self.queue.pop()

    def see_next(self):
        return self.queue[-1]

    def __bool__(self):
        return bool(self.queue)

    def __len__(self):
        return len(self.queue)


class Event:
    """
    An operation related to a craft at a future time.
    """
    event_type = 'generic'

    def __init__(self, craft, time, func=None, args=()):
        self.craft = craft
        self.time = time
        self.func = func
        self.args = args

    @property
    def func_type(self):
        return self.event_type

    def __repr__(self):
        return "<{} of {} at {}d>".format(
            self.event_type,
            self.craft.name if self.craft else 'unknown',
            round(self.time, 2)
        )

    def __lt__(self, other):
        return self.time < other.time

    def __eq__(self, other):
        return self.time == other.time

    def __gt__(self, other):
        return self.time > other.time

    def _do(self):
        """Action to complete at self.time. Returns new event if relevant."""
        if self.func is not None:
            return self.func(*self.args)
        return None


class CloseApproach(Event):
    """
    An event representing a found close approach.
    """
    event_type = 'close_approach'

    def __init__(self, craft, time, func=None, t_search_end=None, dist=0.0, target=None):
        super().__init__(craft, time, func)
        self.t_search_end = t_search_end
        self.dist = dist
        self.target = target

    def _do(self):
        return super()._do()


class MOID(Event):
    """
    Minimum Orbital Intersection Distance geometric event.
    """
    event_type = 'moid'

    def __init__(self, craft, time, target, moid_type='min'):
        super().__init__(craft, time)
        self.target = target
        self.moid_type = moid_type

    def _do(self):
        if self.craft.orbit.is_hyp:
            return None
        self.time += self.craft.period()
        return self


class SoiEnter(Event):
    """
    An event representing an SOI entry.
    """
    event_type = 'soi_enter'

    def __init__(self, craft, time, target, func=None):
        super().__init__(craft, time, func)
        self.target = target

    def new_orbit(self):
        """Run when the event occurs: switch coordinates to target body frame."""
        old_state = self.craft.orbit.get_state(self.time)
        target_state = self.target.orbit.get_state(self.time)
        rel_state = old_state - target_state
        rel_state.mu = self.target.mu
        self.craft.orbit = rel_state.get_orbit()

        if self.craft.parent is not None and hasattr(self.craft.parent, 'crafts'):
            self.craft.parent.crafts.pop(self.craft.name, None)
        self.craft.parent = self.target
        self.target.crafts[self.craft.name] = self.craft

    def _do(self):
        self.new_orbit()
        # Automatically transition frame and schedule SOIExit upon entry

        t_exit = find_escape(self.craft, self.time)
        if t_exit is not None:
            return SoiExit(self.craft, t_exit, target=self.target)
        else:
            return check_orbit(self.craft, self.time)


class SoiExit(Event):
    """
    An event representing an SOI departure.
    """
    event_type = 'soi_exit'

    def __init__(self, craft, time, target=None, func=None):
        super().__init__(craft, time, func)
        self.target = target

    def new_orbit(self):
        """Run when event occurs: switch coordinates to parent's parent frame."""
        exited_body = self.craft.parent
        outer_parent = exited_body.parent
        if outer_parent is None:
            # Sun or top-level body
            return

        old_state = self.craft.orbit.get_state(self.time)
        parent_state = exited_body.orbit.get_state(self.time)
        new_state = parent_state + old_state
        new_state.mu = outer_parent.mu
        self.craft.orbit = new_state.get_orbit()

        exited_body.crafts.pop(self.craft.name, None)
        self.craft.parent = outer_parent
        outer_parent.crafts[self.craft.name] = self.craft

    def _do(self):
        self.new_orbit()
        return check_orbit(self.craft, self.time)
