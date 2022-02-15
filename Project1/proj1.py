from numpy import *
from scipy.linalg import *
import matplotlib.pyplot as plt

class cubic_spline:
    def __init__(self, grid, control_points):
        self.control_points = asarray(control_points)
        self.grid = grid
        
    def __call__(self):
        self.plot()
        pass
    
    def blossom(self,u,grid):
        index = self.hot_interval(grid,u)
        d = self.get_control_points(index)
        dA = self.interpolation(d[0], d[1], u, index-2,index+1)
        dB = self.interpolation(d[1], d[2], u, index-1, index+2)
        dC = self.interpolation(d[2], d[3], u, index, index+3)
        dAB = self.interpolation(dA, dB, u, index-1,index+1)
        dBC = self.interpolation(dB, dC, u,index,index+2)
        su = self.interpolation(dAB, dBC, u,index,index+1)
             
        return su
    
    def interpolation(self,d1,d2,u,left_most,right_most):
        if self.grid[right_most] == self.grid[left_most]:
            alpha = 0        
        else:
            alpha = (self.grid[right_most] - u)/(self.grid[right_most] - self.grid[left_most])
        
        return alpha*d1 + (1-alpha)*d2
    
    def hot_interval(self, grid, u):      
        for i in range(len(grid)):
             if u == self.grid[2]:
                    return 2
             if u <= self.grid[i]:
                    return i-1
    
    def get_control_points(self, index):
        return self.control_points[index-2:index+2]

    #Recursive function for the basis functions
    def bFunc(self,k,j,u):  #Returns function for the j:th basis-function, k is the degree of N and u is the knot sequence
        grid = self.grid.copy()
        
        if(k==0):  
            if ((self.grid[j-1]==self.grid[j])):
                return 0
            elif (self.grid[j-1]<=u)*(u<=self.grid[j]):
                return 1
            else:
                return 0
        else:
            n1=0
            n2=0
            if j > 0:
                if(self.grid[j+k-1] != self.grid[j-1]):
                    n1 = ((u-self.grid[j-1])/(self.grid[j+k-1]-self.grid[j-1]))*self.bFunc(k-1,j,u)
            if(j+k < len(self.grid)):
                if(self.grid[j+k]!=self.grid[j]):
                    n2 = ((self.grid[j+k]-u)/(self.grid[j+k]-self.grid[j]))*self.bFunc(k-1,j+1,u)                    
            N = n1 + n2 
        return N
        
    def Nj(self):
        function = lambda index,u: self.bFunc(3,index,u)
        return function

    def plot(self):
        
        plt.plot(self.control_points[:,0], self.control_points[:,1], '-.*')#Plottar controllpunkerna + polygon
        
        u = linspace(self.grid[0], self.grid[-1], num = 100*len(self.grid))
        u[1] = u[2] = u[0]
        u[-3] = u[-2] = u[-1]
        su = len(u)*[0]
        for k in range(len(su)):
                
                su[k] = self.blossom(u[k],self.grid)
                
        plt.plot(asarray(su)[:,0],asarray(su)[:,1])
        plt.show()
        
if __name__ == '__main__':
    CONTROL = [(-12.73564, 9.03455),
        (-26.77725, 15.89208),
        (-42.12487, 20.57261),
        (-15.34799, 4.57169),
        (-31.72987, 6.85753),
        (-49.14568, 6.85754),
        (-38.09753, -1e-05),
        (-67.92234, -11.10268),
        (-89.47453, -33.30804),
        (-21.44344, -22.31416),
        (-32.16513, -53.33632),
        (-32.16511, -93.06657),
        (-2e-05, -39.83887),
        (10.72167, -70.86103),
        (32.16511, -93.06658),
        (21.55219, -22.31397),
        (51.377, -33.47106),
        (89.47453, -33.47131),
        (15.89191, 0.00025),
        (30.9676, 1.95954),
        (45.22709, 5.87789),
        (14.36797, 3.91883),
        (27.59321, 9.68786),
        (39.67575, 17.30712)]
    KNOTS = linspace(0, 1, len(CONTROL)+2)
    KNOTS[ 1] = KNOTS[ 2] = KNOTS[ 0]
    KNOTS[-3] = KNOTS[-2] = KNOTS[-1]
    
    spline = cubic_spline(KNOTS, CONTROL)
    spline.plot()