import networkx as nx

class CandidateDesign():
    def __init__(self, graph, flow_specific_AC, paths, SC, s_net,delay_constraints):
        self.graph = graph  # Physical network abstraction
        # self.all_flows = all_flows  # Set of all the flows
        self.flow_specific_AC = flow_specific_AC  # Flow_specific arrival curves
        self.paths = paths  # Set of paths for all flows
        self.SC = SC  # Set of servers with service curves
        self.s_net = s_net
        # self.flow_specific = flow_specific  # Flow_specific
        # self.cumulative = cumulative  # cumulative functions at servers
        # self.flow_through_node = flow_through_node  # Set of flows passing through a node
        self.delay_constraints = delay_constraints   # Set of flow-specific delay constraints
        self.time_variable_map = dict()  # t1 : t^si_j
        self.time_variable_rmap = dict()  # t^si_j : t1
        self.time_variable_index = 0
        self.flow_variable_map = dict()  # F
        self.flow_variable_rmap = dict()
        self.flow_variable_index = 0
        self.departure_variable_map = dict()
        self.departure_variable_rmap = dict()
        self.entry_variable_map = dict()
        self.entry_variable_rmap = dict()