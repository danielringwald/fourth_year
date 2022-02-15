from numpy import *
from scipy import *
from scipy.linalg import *

class OptimizationProblem:
    
    def __init__(self, f, grad=None):
        self.f = f
        self.grad = grad
                