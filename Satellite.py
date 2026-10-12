from Orbit import mk_Orbit


class Satellite:
    def __init__(self, body_dict, parent=None):
        """
        A generic object in space with an orbit.
        :param body_dict: parameters
            keys: name, color, orbit (dict, see Orbit.py)
        :param parent: the parent Body object
        """
        self.name = body_dict['name']
        self.parent = parent
        if parent is None:
            self.orbit = None
        else:
            orbit_data = dict(body_dict['orbit'])
            orbit_data['mu'] = parent.mu
            self.orbit = mk_Orbit(orbit_data)
        self.color = body_dict.get('color', 'white')

        # Double-link with parent
        if self.parent is not None and isinstance(self.parent, Body):
            if isinstance(self, Body):
                self.parent.bodies[self.name] = self
            elif isinstance(self, Craft):
                self.parent.crafts[self.name] = self

    def __repr__(self):
        return '<{}>'.format(self.name)


class Body(Satellite):
    """
    A Satellite with its own gravity and satellites.
    """

    def __init__(self, body_dict, parent=None):
        # Initialize dictionary caches before super().__init__ links them
        self.bodies = {}
        self.crafts = {}
        self.satellites = {}
        self.size = float(body_dict['size'])
        self.mu = float(body_dict['mu'])

        super().__init__(body_dict, parent)

        if parent is None:
            self.soi = self.size * 2.0
        else:
            semi_maj = abs(self.orbit.semi())
            self.soi = semi_maj * (self.mu / self.parent.mu) ** (2.0 / 5.0)
            if self.parent.parent is None:
                if self.orbit.is_hyp:
                    self.parent.soi = max(self.parent.soi, 1.5 * self.orbit.peri)
                else:
                    self.parent.soi = max(self.parent.soi, 1.5 * self.orbit.apo())


class Craft(Satellite):
    """
    A spacecraft that can change orbit via burns and maneuver planning.
    """
    res_types = ['delta-v']

    def __init__(self, body_dict, parent=None):
        super().__init__(body_dict, parent)
        self.resources = {}
        for r in self.res_types:
            self.resources[r] = (0.0, 0.0)

    def period(self):
        if self.orbit is None or self.orbit.is_hyp:
            return None
        return self.orbit.period()

    def program(self, time, next_time, next_type):
        """
        Optional mission plan callback. Returns an Event or None.
        """
        return None

    def adjust_res(self, name, quant_adj, cap_adj=0.0):
        """
        Change a resource quantity.
        """
        if name in self.resources:
            curr, cap = self.resources[name]
            cap += cap_adj
        else:
            if name not in self.res_types:
                self.res_types.append(name)
            curr, cap = 0.0, cap_adj

        self.resources[name] = (curr, cap)

        prev = curr
        if quant_adj == 'full':
            self.resources[name] = (cap, cap)
        elif quant_adj == 'empty':
            self.resources[name] = (0.0, cap)
        else:
            new = min(max(curr + quant_adj, 0.0), cap)
            self.resources[name] = (new, cap)

        curr = self.resources[name][0]
        return curr - prev

    def burn(self, dv_vect, t):
        """
        Apply a burn to reach a new orbit.
        :param dv_vect: vector representing the burn [magnitude, angle]
        :param t: time of burn [days]
        """
        dv_vect = list(dv_vect)
        dv_available = self.adjust_res('delta-v', -dv_vect[0])
        dv_vect[0] = -dv_available
        state = self.orbit.get_state(t)
        state.boost(dv_vect)
        self.orbit = state.get_orbit()
