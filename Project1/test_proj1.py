from numpy import *
import numpy as np
from numpy.linalg import *
import unittest
#from test_matrices import test_matrix1
from proj1 import cubic_spline

class Test(unittest.TestCase):
    
    def setUp(self):
        self.testCP = [(-12.73564, 9.03455),
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
        self.testGrid = linspace(0, 1, 26)
        self.testGrid[1] = self.testGrid[2] = self.testGrid[0]
        self.testGrid[-3] = self.testGrid[-2] = self.testGrid[-1]
        self.spline = cubic_spline(self.testGrid, self.testCP)

    def test_evaluation(self):
        expected = array([-31.90219167, 6.47655833])
        u = 0.2
        evaluated_spline = self.spline.blossom(u, self.testGrid)
        for i in range(len(expected)):
            self.assertAlmostEqual(evaluated_spline[i], expected[i])

    def test_basis_func(self):
        basis_func = self.spline.Nj()
        u = linspace(self.testGrid[0], self.testGrid[-1], num=100*len(self.testGrid))
        result = zeros(len(u))
        
        for i in range(len(self.testCP)):
        #for i in range(5):
            for j in range(len(u)):
                result[j] = result[j] + basis_func(i,u[j])
        expected = ones(len(result))
        
        for i in range(len(expected)):
            self.assertAlmostEqual(expected[i], result[i])

    def test_basis_equal_blossom(self):
        u = linspace(self.testGrid[0], self.testGrid[-1], num=100*len(self.testGrid))
        s = array([])
        basis_func = self.spline.Nj()
        for i in range(len(self.testCP)):
            su = self.spline.blossom(u[i], self.testGrid)
            N = zeros(len(u))
            for j in range(len(u)):
                N[j] = basis_func(i, u[j])
            
            x_val_basis = sum(self.testCP[i][0]*N)
            y_val_basis = sum(self.testCP[i][1]*N)
            
            self.assertAlmostEqual(su[0], x_val_basis)
            self.assertAlmostEqual(su[1], y_val_basis)

if __name__ == '__main__':
    unittest.main()

# if _name_ == '_main_':     
#     testGrid = array([0.  , 0.  , 0.  , 0.12, 0.16, 0.2 , 0.24, 0.28, 0.32, 0.36, 0.4 ,
#                       0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.72, 0.76, 0.8 , 0.84,
#                       0.88, 1.  , 1.  , 1.  ])
#     testCP = [(-12.73564, 9.03455),
#       (-26.77725, 15.89208),
#       (-42.12487, 20.57261),
#       (-15.34799, 4.57169),
#       (-31.72987, 6.85753),
#       (-49.14568, 6.85754),
#       (-38.09753, -1e-05),
#       (-67.92234, -11.10268),
#       (-89.47453, -33.30804),
#       (-21.44344, -22.31416),
#       (-32.16513, -53.33632),
#       (-32.16511, -93.06657),
#       (-2e-05, -39.83887),
#       (10.72167, -70.86103),
#       (32.16511, -93.06658),
#       (21.55219, -22.31397),
#       (51.377, -33.47106),
#       (89.47453, -33.47131),
#       (15.89191, 0.00025),
#       (30.9676, 1.95954),
#       (45.22709, 5.87789),
#       (14.36797, 3.91883),
#       (27.59321, 9.68786),
#       (39.67575, 17.30712)]
             
     
     
 
#     testSpline = cubic_spline(testGrid,testCP)
#     test_N_grid = testGrid[0:5]
#     f = testSpline.Nj()
 
#     u = linspace(test_N_grid[0],test_N_grid[-1],num=100*len(test_N_grid))
 
#     N = (len(u))*[0]
 
#     # for j in range(5):
#     #     for i in range(len(u)):
#     #         N[i] = f(j+1,u[i])
 
#     #     plt.plot(u,N)
#     # plt.show()
#     N = (len(u))*[0]
#     for i in range(len(u)):
#         N[i] = (f(0,u[i]))
#     plt.plot(u,N)
#     plt.show()
#     testSpline.plot()