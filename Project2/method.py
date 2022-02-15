from scipy import *
from scipy.optimize import minimize_scalar
from scipy.linalg import *
from numpy import *
import matplotlib.pyplot as plt
from numdifftools import Hessian

"""
  
"""

class OptimizationMethod:
    """
    Find minimum of function.
    
    Parameters:     opt_problem: OptimizationProblem object
                    
                    x0: Initial guess
    
    Returns:        out: 
    """
    
    def __init__(self, line_search, tol=1e-5, show_hessian=False):
        
        self.line_search_bool = line_search      #Bool: True for exact, False for inexact
        self.tol = tol
        self.show_hessian = show_hessian
        
    def __call__(self, opt_problem, initial_guess):
        """
        

        Parameters
        ----------
        opt_problem : OptimizationProblem
            Problem to be optimized.
        initial_guess : Array of float
            Point of the initial guess.

        Returns
        -------
        sol : Array of float
            The point of the minimum.

        """
        self.opt_problem = opt_problem
        self.point = initial_guess
        sol = self.solver()
        
        return sol
    
    def get_gradient(self, point):
        """
        

        Parameters
        ----------
        point : Array of float
            Location of the point.

        Returns
        -------
        gradient : Array of float
            Gradient at the specified point.

        """
        if self.opt_problem.grad != None:
            return self.opt_problem.grad(point)
        
        n = len(point)
        h = self.told
        f = self.opt_problem.f
        gradient = zeros(n)
        
        for i in range(n):
            
            x_copy = point
            x_plus_h = x_copy
            x_plus_h[i] += h
            x_minus_h = x_copy
            x_minus_h[i] -= h
            
            gradient[i] = ( f(x_plus_h) - f(x_minus_h) ) / (2 * h)
        
        return gradient
    
    def direction(self, point):
        """
        

        Parameters
        ----------
        point : Array of float
            Location of the point.

        Returns
        -------
        Array of float
            Direction.

        """
        
        hi = self.tol
        hj= self.tol
        if hasattr(point, "_len_"):
            n = len(point)
        else:
            n = 1
        hessian = zeros((n,n))
        
        if self.opt_problem.grad == None:
            gradient = self.get_gradient(f, initial_guess)
        else:
            gradient = self.opt_problem.grad     #if grad is not approximated it is 
        hessian = Hessian(self.opt_problem.f)
        hessian = hessian(point)
                    
        return - inv( 1 / 2 * (hessian + transpose(hessian))) @ gradient(point)
            
    def line_search(self, point, direction):
        """
        

        Parameters
        ----------
        point : Array of float
            DESCRIPTION.
        direction : Array of float
            DESCRIPTION.

        Raises
        ------
        Exception
            If scpipy.optimize.minimize_scalar does not return succes == True.

        Returns
        -------
        Float
            alpha.

        """
        f = self.opt_problem.f
        max_iterations = 10
        point = array(point)
        direction = array(direction)
        
        if self.line_search_bool:    #Exact line search
            def alpha_f(alpha):    
                return f(point + alpha * direction)
            
            minimized = minimize_scalar(alpha_f)
            if minimized.success:
                return minimized.x
            else:
                raise Exception
        
        #Inexact line search using Powell-wolfe algorithm
        sigma = 0.01         #[0,1/2]
        rho = 0.9           #[sigma,1]
        alphaLB = 2         #alpha lower boundary
        iterations = 0    
        h = self.tol
        
        def phi(alpha):
            return f(point + alpha*direction)
        
        def get_derivative(alpha):
            return (phi(alpha + h) - phi(alpha - h) ) / (2 * h)
        
        #while alphaLB does not fulfill Armijo-rule
        while not phi(alphaLB) <= phi(0) + sigma*alphaLB*get_derivative(0):
            alphaLB = alphaLB/2
        
        alphaUB = alphaLB       #alpha upper boundary
        
        #while alphaUB fulfills Armijo-rule
        while phi(alphaUB) <= phi(0) + sigma*alphaUB*get_derivative(0):
            alphaUB = 2*alphaUB
           
        #while alphaLB does not fulfill the second Powell-wolfe rule
        while not get_derivative(alphaLB) >= rho*get_derivative(0):
           
            if iterations >= max_iterations:
                break
                #raise RuntimeError
            iterations += 1
            
            alpha0 = (alphaUB + alphaLB)/2
            #If alpha0 satisfies Armijo-rule
            if phi(alpha0) <= phi(0) + sigma*alpha0*get_derivative(0):
                alphaLB = alpha0
            else:
                alphaUB = alpha0

        return alphaLB
    
#Classic Newton Method
class NewtonMethod(OptimizationMethod):
    """
    Implements step() and solver() in OptimizationMethod with classical Newton
    """
    
    def solver(self):
       """
        

        Returns
        -------
        minimum : Array of float
            Minimum point.
        steps : Array of float
            Point of each step.

        """
        
        gradient = self.opt_problem.grad
        tol = self.tol
        residual_crit = 1
        cauchy_crit = 1
        x_new = self.point
        steps = []
        steps.append(x_new)
        criterion_bool = True
        while criterion_bool:
            x_old = x_new
            x_new = self.step(x_old)
            residual_crit = norm(gradient(x_new))
            cauchy_crit = norm(x_new - x_old)
            criterion_bool = residual_crit > tol and cauchy_crit > tol
            steps.append(x_new)
            
        minimum = x_new
        return minimum, steps
    
    def step(self, x_old):
        """
        

        Parameters
        ----------
        x_old : Array of float
            Old point.

        Returns
        -------
        Array of float
            New point.

        """
        step = self.direction(x_old)
        
        return x_old + self.line_search(x_old, self.direction(x_old)) * step
    
#Quasi Newton Methods
class QN(OptimizationMethod):
            
    def step(self, point):
        """
        

        Parameters
        ----------
        point : Array of float
            Initial point.

        Returns
        -------
        x_new : Array of float
            New point.

        """
        x = array(point)
        gradient = self.get_gradient(x)
        try:
            dirr = -self.Hinv @ gradient
        except:
            dirr = -self.Hinv * gradient
        alpha = self.line_search(point,dirr)
        x_new = x + alpha*dirr
        gradient_new = self.get_gradient(x_new)
        delta = x_new - x
        gamma = array(gradient_new) - array(gradient)
        
        self.Hinv = self.getHinv(self.Hinv, delta, gamma)
        
        return x_new
        
    def solver(self):
        """
        

        Returns
        -------
        minimum : Array of float
            Minimum point.
        steps : Array of float
            Point of each step.

        """
        point = self.point
        
        #Can handle 1d
        try:
            self.Hinv = eye(len(point))
        except TypeError:
            self.Hinv = 1
        cHess = Hessian(self.opt_problem.f)
        self.Hinv = inv(cHess(self.point))
        
        gradient = self.get_gradient
        tol = self.tol
        residual_crit = 1
        cauchy_crit = 1
        x_new = point
        steps = []
        steps.append(x_new)
        criterion_bool = True
        
        hessians = []
        counter = 0
        
        while criterion_bool:
            if self.show_hessian:
                correct_hess = inv(cHess(x_new))
                hessians.append(norm(self.Hinv-correct_hess))
                counter += 1
                
            x_old = x_new
            x_new = self.step(x_old)
            
            residual_crit = norm(gradient(x_new))
            cauchy_crit = norm(x_new - x_old)
            
            criterion_bool = residual_crit > tol and cauchy_crit > tol
            steps.append(x_new)
        
        if self.show_hessian:
            plt.plot(arange(counter)[1:], log(hessians[1:]))
            plt.title('Norm of difference between BFGS and finite difference calculation of Hessian.')
            plt.xlabel('k')
            plt.ylabel('log of norm')
            plt.show()
        
        minimum = x_new
            
        return minimum, steps    
    
class GoodBroyden(QN):
    
    def getHinv(self, Hinv, delta, gamma):
        """
        

        Parameters
        ----------
        Hinv : Array of float
            Old inverse hessian.
        delta : float
            point_new - point_old.
        gamma : float
            gradient_new - gradient_old.

        Returns
        -------
        Array of float
            Updated inverse hessian.

        """
        
        return Hinv + (outer(delta-Hinv@gamma,delta)/(delta@(Hinv@gamma)))@Hinv

class BadBroyden(QN):
    
    def getHinv(self,Hinv, delta, gamma):
        """
        

        Parameters
        ----------
        Hinv : Array of float
            Old inverse hessian.
        delta : float
            point_new - point_old.
        gamma : float
            gradient_new - gradient_old.

        Returns
        -------
        Array of float
            Updated inverse hessian.

        """
        t = outer(delta - Hinv @ gamma, gamma)
        n = dot(gamma,gamma)
        
        return Hinv + t/n
    
class SymmetricBroyden(QN):
    
    def getHinv(self,Hinv, delta, gamma):
        """
        

        Parameters
        ----------
        Hinv : Array of float
            Old inverse hessian.
        delta : float
            point_new - point_old.
        gamma : float
            gradient_new - gradient_old.

        Returns
        -------
        Array of float
            Updated inverse hessian.

        """
        u = delta - Hinv @ gamma
        a = 1 / (u @ gamma)
        return Hinv + a*outer(u, u)
    
class DFP(QN):
    
    def getHinv(self,Hinv, delta, gamma):
        """
        

        Parameters
        ----------
       Hinv : Array of float
            Old inverse hessian.
        delta : float
            point_new - point_old.
        gamma : float
            gradient_new - gradient_old.

        Returns
        -------
        Array of float
            Updated inverse hessian.

        """
        temp1 = outer(delta,delta)/dot(delta,gamma)
        temp2 = outer(gamma,gamma)
        temp3 = Hinv@temp2@Hinv
        temp4 = dot(gamma,(Hinv @ gamma))
        
        return Hinv + temp1 - temp3/temp4
    
class BFGS(QN):
   
    def getHinv(self,Hinv, delta, gamma):
        """
        

        Parameters
        ----------
        Hinv : Array of float
            Old inverse hessian.
        delta : float
            point_new - point_old.
        gamma : float
            gradient_new - gradient_old.

        Returns
        -------
        Array of float
            Updated inverse hessian.

        """
        temp3 = (outer(delta,gamma)@Hinv+outer(Hinv@gamma,delta))/(delta@gamma)
        temp2 = outer(delta,delta)/(delta@gamma)
        temp1 = 1 + (dot(gamma,dot(Hinv,gamma)))/(dot(delta,gamma))
        Hinv_new = Hinv + dot(temp1,temp2) - temp3    
        return Hinv_new