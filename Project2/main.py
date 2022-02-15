from method import*
from problem import OptimizationProblem
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
from chebyquad_problem import *
import pretty_errors

class Main:
    
    def __init__(self, func='rosen'):
        self.func = func
    
    def __call__(self, line_search=True, initial_guess=[0,0], tol=1e-5):
        if self.func == 'rosen':
            #Rosenbrock
            f = lambda x: 100 * (x[1] - x[0]**2)**2 + ( 1 - x[0])**2
            grad = lambda x: [2*100 * (-2*x[0])*(x[1]-x[0]**2) -2*(1-x[0]), 2*100 * (x[1] - x[0]**2)]
        
            problem = OptimizationProblem(f, grad)
        
            #Classic Newton method
            method = NewtonMethod(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'Classic Newton')
            print('Min classic Newton: ', minimum)
             
            #Good Broyden
            method = GoodBroyden(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'Good Broyden')
            print("Min Good Broyden: ", minimum)
            
            #Bad Broyden
            method = BadBroyden(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'Bad Broyden')
            print("Min Bad Broyden: ", minimum)
            
            #Symmetric Broyden
            method = SymmetricBroyden(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'Symmetric Broyden')
            print("Min Symmetric Broyden: ", minimum)
            
            #DFP
            method = DFP(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'DFP')
            print("Min DFP: ", minimum)
            
            #BFGS
            method = BFGS(line_search)
            minimum, steps = method(problem, initial_guess)
            steps = array(steps)
            self._plot(steps, 'BFGS')
            print("Min BFGS: ", minimum)
            
        
        if func == 'cheby':
            f = chebyquad
            grad = gradchebyquad
            problem = OptimizationProblem(f, grad)
            
            #BFGS
            method = BFGS(line_search)
            minimum, steps = method(problem, initial_guess=linspace(0,1,4))
            steps = array(steps)
            
            print("Min BFGS: ", minimum)
            
        
    def _plot(self, steps, method):
        v_func = plt.vectorize(plot_fun)
        x, y = plt.meshgrid(linspace(-0.5, 2, 1000), linspace(-0.5, 4, 1000))
        fig, ax = plt.subplots(1)
        ax.plot(steps[1:-1,0] , steps[1:-1,1], '.',c='black')
        ax.plot(steps[0,0] , steps[0,1], '.', c='r')
        ax.plot(steps[-1,0] , steps[-1,1], '.',c='g')
        CS = ax.contour(x, y, v_func(x, y), [1,3.831,14.678,56.234,215.443,825.404])
        ax.clabel(CS, inline=1, fontsize=9)
        ax.set_title(f'Method: {method}')
        plt.show()
     
if __name__ == '__main__':
     rosen_plot = Main()
     rosen_plot()
     cheby_plot = Main('cheby')
     cheby_plot()
     