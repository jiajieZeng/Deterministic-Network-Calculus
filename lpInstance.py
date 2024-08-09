
class LPInstance():
    def __init__(self, objective, flow_constraints, time_constraints, topo_list, C, C_prime, id, time_map, flow_map, obj, ne, neb, eq, eqb, bound,map):
        self.objective = objective
        self.flow_constraints = list(flow_constraints)
        self.time_constraints = list(time_constraints)
        self.topo_list = topo_list
        self.C = C
        self.C_prime = C_prime
        self.id = id
        self.time_map = time_map
        self.flow_map = flow_map
        self.obj = obj
        self.ne = ne
        self.neb = neb
        self.eq = eq
        self.eqb = eqb
        self.bound = bound
        self.map = map

