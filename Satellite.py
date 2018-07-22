from Orbit_class import mk_Orbit

class Satellite:
    def __init__(self, body_dict, parent=None):
        """
        A generic celestial body
        :param params: A dictionary?
        """
        self.name = body_dict['name']
        self.orbit = mk_Orbit(body_dict['orbit'])
        self.size = body_dict['size']
        self.color = body_dict['color']
        self.parent = parent
        self.mu = body_dict['mu']
        if self.parent is not None:
            self.parent.satellites[self.name] = self

    def __repr__(self):
        return '<{}>'.format(self.name)


class Body(Satellite):
    def __init__(self, body_dict, parent=None):
        Satellite.__init__(self, body_dict, parent)
        self.satellites = {}


class Craft(Satellite):
    res_types = ['delta-v']

    def __init__(self, body_dict, parent=None):
        Satellite.__init__(self, body_dict, parent)
        self.resources ={}
        for r in self.res_types:
            self.resources[r] = (0, 0)
        self.parts = []
        self.mass = 0

    def build(self, parts):
        for part in parts:
            self.parts.append(part)
            for res in part.res_types:
                part_curr, part_cap = part.resources[res]
                self.adjust_res(res, part_curr, part_cap)

    def adjust_res(self, name, curr_adj, cap_adj=0):
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
        if curr_adj == 'full':
            self.resources[name] = (cap, cap)
        elif curr_adj == 'empty':
            self.resources[name] = (0.0, cap)
        else:
            new = min(max(curr + curr_adj, 0.0), cap)
            self.resources[name] = (new, cap)
        curr = self.resources[name][0]
        return curr - prev #true change in the resource

    def burn(self, dv_vect, t):
        """
        Apply a burn to reach a new orbit.
        :param dv_vect: polar 2-item list
        :param t: time in days
        """
        #manage resources
        dv_vect = list(dv_vect) #make non-mutable
        dv_available = self.adjust_res('delta-v', -dv_vect[0])
        dv_vect[0] = -dv_available
        #change the orbit
        state = self.orbit.get_state(t)
        state.boost(dv_vect)
        self.orbit = state.get_orbit()


class Part(Craft):

    def __init__(self, name, resources):
        self.res_types = resources.keys()
        self.resources = {}
        for res in self.res_types:
            self.resources[res] = (0, 0)

