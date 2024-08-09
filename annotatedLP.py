class AnnotatedLP():
    def __init__(self, val, topo_list, C, C_prime, C_id, time_map, flow_map, new_edges, map, digits, time_constraints, flow_constraints):
        self.val = val
        self.topo_list = topo_list
        self.C = C
        self.C_prime = C_prime
        self.C_id = C_id
        self.time_map = time_map
        self.new_edges = new_edges
        self.map = map
        self.digits = digits
        self.flow_map = flow_map
        self.time_constraints = time_constraints
        self.flow_constraints = flow_constraints

