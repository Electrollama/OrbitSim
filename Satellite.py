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
        if parent is None:
            self.orbit = None
        else:
            body_dict['orbit']['mu'] = parent.mu
            self.orbit = mk_Orbit(body_dict['orbit'])
        self.color = body_dict['color']
        # Double-link with parent
        self.parent = parent
        if (self.parent is not None and
                isinstance(self.parent, Body)):
            if isinstance(self, Body):
                self.parent.bodies[self.name] = self
            elif isinstance(self, Craft):
                self.parent.crafts[self.name] = self

    def __repr__(self):
        return '<{}>'.format(self.name)


class Body(Satellite):
    """
    A Satellite with its own gravity and satellites.
    surface gravity > .01g i.e.
    density[g/cm^2] * radius[m] > 3.5e5
    If density is not known, assume 1.0
    plus Gilly as an exception?
    """

    def __init__(self, body_dict, parent=None):
        """
        load a Body object
        :param body_dict:
            keys: name, color, size, mu, orbit (dict, see Orbit.py)
        :param parent: a loaded
        """
        Satellite.__init__(self, body_dict, parent)
        self.bodies = {}  # populated when a satellite is initialized
        self.satellites = {}
        self.size = body_dict['size']
        self.mu = body_dict['mu']
        if parent is None:
            self.soi = self.size * 2.
        else:
            self.soi = self.orbit.semi() * (self.mu / self.parent.mu) ** (2 / 5)
            if self.parent.parent is None:  # adjust relevant range of the system
                # update parent soi
                if self.orbit.is_hyp:
                    parent.soi = max(parent.soi, 1.5 * self.orbit.peri)
                else:
                    parent.soi = max(parent.soi, 1.5 * self.orbit.apo())


class Craft(Satellite):
    """
    A man-made satellite (space craft) that can change orbit
    """
    res_types = ['delta-v', ]

    def __init__(self, body_dict, parent=None):
        Satellite.__init__(self, body_dict, parent)
        self.resources = {}
        for r in self.res_types:
            self.resources[r] = (0, 0)

    def adjust_res(self, name, quant_adj, cap_adj=0):
        """
        Change a resource quantity
        :param name: name of the resource
        :param quant_adj: change in quantity
                (value, 'full' or 'empty')
        :param cap_adj: change in max quantity
        :return: The resulting change in quantity
        """
        # adjust capacity
        if name in self.res_types:
            curr, cap = self.resources[name]
            cap += cap_adj
        else:
            self.res_types.append(name)
            curr, cap = 0.0, cap_adj
        self.resources[name] = (curr, cap)
        # adjust current quantity
        prev = curr
        if quant_adj == 'full':
            self.resources[name] = (cap, cap)
        elif quant_adj == 'empty':
            self.resources[name] = (0.0, cap)
        else:
            new = min(max(curr + quant_adj, 0.0), cap)
            self.resources[name] = (new, cap)
        curr = self.resources[name][0]
        return curr - prev  # true change in the resource

    def burn(self, dv_vect, t):
        """
        Apply a burn to reach a new orbit.
        :param dv_vect: vector representing the burn (2-value list)
        :param t: time of burn [days]
        """
        # manage resources
        dv_vect = list(dv_vect)  # make non-mutable
        dv_available = self.adjust_res('delta-v', -dv_vect[0])
        dv_vect[0] = -dv_available
        # change the orbit
        state = self.orbit.get_state(t)
        state.boost(dv_vect)
        self.orbit = state.get_orbit()
