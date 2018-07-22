from Orbit import mk_Orbit

class Satellite:
    def __init__(self, body_dict, parent=None):
        """
        A generic object in space with an orbit.
        :param body_dict: See load.py
        :param parent: the parent Body object
        """
        self.name = body_dict['name']
        self.orbit = mk_Orbit(body_dict['orbit'])
        self.color = body_dict['color']
        #Double-link with parent
        self.parent = parent
        if (self.parent is not None and
                isinstance(self.parent, Body)):
            self.parent.satellites[self.name] = self

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
        Satellite.__init__(self, body_dict, parent)
        self.satellites = {} #populated when a satellite is initialized
        self.size = body_dict['size']
        self.mu = body_dict['mu']


class Craft(Satellite):
    res_types = ['delta-v',]

    def __init__(self, body_dict, parent=None):
        Satellite.__init__(self, body_dict, parent)
        self.resources ={}
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
        #adjust capacity
        if name in self.res_types:
            curr, cap = self.resources[name]
            cap += cap_adj
        else:
            self.res_types.append(name)
            curr, cap = 0.0, cap_adj
        self.resources[name] = (curr, cap)
        #adjust current quantity
        prev = curr
        if quant_adj == 'full':
            self.resources[name] = (cap, cap)
        elif quant_adj == 'empty':
            self.resources[name] = (0.0, cap)
        else:
            new = min(max(curr + quant_adj, 0.0), cap)
            self.resources[name] = (new, cap)
        curr = self.resources[name][0]
        return curr - prev #true change in the resource

    def burn(self, dv_vect, t):
        """
        Apply a burn to reach a new orbit.
        :param dv_vect: vector representing the burn (2-value list)
        :param t: time of burn [days]
        """
        #manage resources
        dv_vect = list(dv_vect) #make non-mutable
        dv_available = self.adjust_res('delta-v', -dv_vect[0])
        dv_vect[0] = -dv_available
        #change the orbit
        state = self.orbit.get_state(t)
        state.boost(dv_vect)
        self.orbit = state.get_orbit()