import networkx as nx
from candidateDesign import CandidateDesign

class AnnotatedDesign():
    def __init__(self, C, T_C, F_C, O, E_C, D_C):
        self.C = C
        self.T_C = T_C
        self.F_C = F_C
        self.O = O
        self.E_C = E_C
        self.D_C = D_C