"""
Created on Wed Oct  6 15:05:44 2021

@author: Rungeng
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve 
from mpi4py import MPI 

class Dirichlet():
    def __init__(self, u_left, u_right, vk, **kargs):
        """
        --------------------------------------------------
        Parameter:
            u_left: Boundary from left room
            u_right: Boudary from right room 
            vk: Previous v(k) to solve v(k+1) using relaxation
            num: Number of mesh grid points
            l: Lenth of the square, default = 2
            w: Width of the suqare, default = 1
            T_top: Tempreture of the top wall
            T_bottom: Tempreture of the bottom wall
            T_left: Tempreture of the left wall
            T_right: Tempreture of the right wall
            omega: Relaxation parameter
            size: Csr_martrix size m*m
        --------------------------------------------------
        """
        self.l = 2
        self.w = 1
        self.u_left = u_left
        self.u_right = u_right
        self.vk = vk
        self.num = kargs.setdefault('num', 20)
        self.T_top = kargs.setdefault('T_top', 40)
        self.T_bottom = kargs.setdefault('T_bottom', 5)
        self.T_left = kargs.setdefault('T_left', 15)
        self.T_right = kargs.setdefault('T_right', 15)
        self.omega = kargs.setdefault('omega', 0.8)
    def _csr_matrix(self):
        rowind = []
        colind = []
        values = []
        num = self.num 
        size = (self.l * num - 1)  * (self.w * num - 1)
        # Depth-First-Search
        def dfs(i,j):
            nonlocal rowind, colind, values
            if i >= size and j >= size:
                return 
            dfs(i+1, j+1)
            rowind.append(i)
            colind.append(j)
            values.append(-4)
            
            if 0 <= 1 + j < size and (j + 1) % (num - 1) != 0: 
                rowind.append(i)
                colind.append(j + 1)
                values.append(1)
            if 0 <= -1 + j < size and j  % (num - 1) != 0:
                rowind.append(i)
                colind.append(j - 1)
                values.append(1)        
            for tmp_j in [self.w * num - 1, 1 - self.w * num]:
                if 0 <= tmp_j + j < size:
                    rowind.append(i)
                    colind.append(j + tmp_j)
                    values.append(1)
            return 
        dfs(0,0)
        return csr_matrix((values, (rowind, colind)))
    def derivative(self, u_left_next, u_right_next):
        u_left = self.u_left 
        u_right = self.u_right
        num = self.num 
        derivative_l = abs((u_left_next - u_left)) * num
        derivative_r = abs((u_right_next - u_right)) * num 
        return derivative_l, derivative_r
    def _get_vkk(self):
        omega = self.omega
        vk = self.vk
        num = self.num 
        # find A matrix
        A =  num**2 * self._csr_matrix()
        # find delta_u matrix
        end_row = self.l * num - 1
        end_column = self.w * num - 1
        delta_u = np.zeros((end_row, end_column))
        delta_u[0, :] -= self.T_top
        delta_u[end_row - 1, :] -= self.T_bottom
        delta_u[0 : end_row//2, 0] -= self.T_left
        delta_u[end_row//2 : end_row, end_column - 1] -= self.T_right
        delta_u[end_row//2 : end_row, 0] -= self.u_left
        delta_u[0 : end_row//2, end_column - 1] -= self.u_right
        delta_u = delta_u.flatten()
        # find v(k+1) matrix
        vkk = spsolve(A, (num**2) * delta_u)
        # relaxation
        vkk = vkk * omega + (1 - omega) * vk
        A = A.toarray()
        u_new = (np.dot(A, vkk) / num**2 ).reshape(end_row, end_column)
        u_right_next = -u_new[0 : end_row//2, end_column - 1]
        u_right_next[0] -= self.T_top
        u_left_next = -u_new[end_row//2: end_row, 0]
        u_left_next[-1] -= self.T_bottom
        derivative_l, derivative_r = self.derivative(u_left_next, u_right_next)
        return vkk, derivative_l, derivative_r
    
class Neumann():
    def __init__(self, vk, u, **kargs):
        self.vk = vk
        self.n = kargs.setdefault('num', 20)
        self.m = len(u)
        self.l = 1
        self.w = 1
        self.T_normal = kargs.setdefault('T_normal', 15)
        self.T_heat = kargs.setdefault('T_heat', 40)
        self.u = u  
        self.omega = kargs.setdefault('omega', 0.8)
        self.size = (self.w * (self.n - 1)) * (self.l * (self.m))
    def _csr_matrix(self):
        n = self.n
        m= self.m
        size = (self.w * (n)) * (self.l * m)
        rowind = []
        colind = []
        values = []
         
        # Depth-First-Search
        def dfs(i,j):
            nonlocal rowind, colind, values
            if i >= size and j >= size:
                return 
            dfs(i+1, j+1)
            if (i + 1) % n  == 0:
                rowind.append(i)
                colind.append(j)
                values.append(-3)
            else:
                rowind.append(i)
                colind.append(j)
                values.append(-4)
            if 0 <= 1 + j < size and (j + 1) % n != 0: 
                rowind.append(i)
                colind.append(j + 1)
                values.append(1)
            if 0 <= -1 + j < size and j  % n != 0:
                rowind.append(i)
                colind.append(j - 1)
                values.append(1)        
            for tmp_j in [self.w * n, - self.w * n]:
                if 0 <= tmp_j + j < size:
                    rowind.append(i)
                    colind.append(j + tmp_j)
                    values.append(1)
            return 
        dfs(0,0)
        return csr_matrix((values, (rowind, colind)))
  
        
    def _get_vkk(self):
        omega = self.omega
        vk = self.vk
        n = self.n
        m = self.m
        # find A matrix
        A = self._csr_matrix() * (n ** 2)
        # find delta_u matrix
        end_row = self.l * m 
        end_column = self.w * n 
    
        delta_u = np.zeros((end_row, end_column))
        delta_u[0, :] -= self.T_normal * n ** 2
        delta_u[end_row - 1, :] -= self.T_normal * n ** 2
        delta_u[:, 0] -= self.T_heat * n ** 2
        delta_u[:, end_column - 1] -= (self.u) * n
        delta_u = delta_u.flatten()
        vkk = spsolve (A , delta_u)
        vkk = vkk * omega + (1 - omega) * vk
        
        u_next = vkk.reshape((end_row, end_column))[:, end_column - 1]
        return vkk, u_next

def _plot(Room, fig, size, origin):
    ax = fig.add_axes(size)
    plt.imshow(Room, origin = origin, interpolation = 'gaussian')
    ax.set_aspect('auto')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpi = True
    deltax = 1 / 20
    num = int(1 / deltax)
    iteration = 10
    parameter = {'num' : num,
                 'l' : 2,
                 'w' : 1,
                 'T_top' : 40,
                 'T_bottom' : 5,
                 'T_left' : 15,
                 'T_right' : 15,
                 'omega' : 0.8
                 }
    u_right = np.array([15.] * (num - 1), dtype = float)
    u_left = np.array([15.] * num , dtype = float)
    v0 = np.full([2 * num - 1, num - 1], 10.).flatten()
    v1 = np.full([num, num ], 10.).flatten() 
    v2 = np.full([num - 1, num], 10.).flatten() 
    room0 = v0
    rooml = v1
    roomr = v2 
    
    if not mpi:
        for i in range(iteration):
            D = Dirichlet(u_left,u_right,v0, num = num)
            b = D._get_vkk()
            N1 = Neumann(v1, b[1],num = num)
            N2 = Neumann(v2, b[2],num = num)   
            d1 = N1._get_vkk()
            d2 = N2._get_vkk()
            t1 = N1._csr_matrix().toarray()
            t2 = N2._csr_matrix().toarray()
            u_left = d1[1]
            u_right = d2[1]
            v0 = b[0]
            v1 = d1[0]
            v2 = d2[0]
            room0 = v0.reshape((2 * num - 1, num - 1))
            rooml =v1.reshape((num, num))
            roomr =np.fliplr(v2.reshape((num - 1, num)))
            fig = plt.figure(figsize=(6, 4))                
            ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
            ax.set_xticks([0,1,2,3])
            ax.set_yticks([0,1,2])
            ax.set_title('Temperature Distribution' + ' %s' %(i))
            size = [0.101,0.1,0.267,0.4]  #left,bottom,width,height 
            _plot(rooml, fig, size, origin = 'lower')
            size = [0.368,0.1,0.267,0.8]  #left,bottom,width,height 
            _plot(room0, fig, size, origin = 'upper')
            size = [0.635,0.5,0.333,0.4]  #left,bottom,width,height 
            _plot(roomr, fig, size, origin = 'lower')
            plt.colorbar(orientation='vertical')
            plt.show()
    if mpi:
        for i in range(iteration):
            if rank == 0:
                D = Dirichlet(u_left, u_right, v0, num = num)
                b = D._get_vkk()
                v0 = np.array(b[0])
                comm.Send([b[1], MPI.DOUBLE], dest = 1, tag = 31)
                comm.Send([b[2], MPI.DOUBLE], dest = 2, tag = 41)
                
            if rank == 1:
                comm.Recv([u_left, MPI.DOUBLE], source = 0, tag = 31)
                N1 = Neumann(v1, u_left, num = num)
                d1 = N1._get_vkk()
                v1 = np.array(d1[0])
                u_left =  np.array(d1[1])
                comm.Send([u_left, MPI.DOUBLE], dest = 0, tag = 33)
                comm.Send([v1, MPI.DOUBLE], dest = 0, tag = 666)
            if rank == 2:
                comm.Recv([u_right, MPI.DOUBLE], source = 0, tag = 41)
                N2 = Neumann(v2, u_right, num = num)
                d2 = N2._get_vkk()
                v2 = np.array(d2[0])
                u_right = np.array(d2[1])
                comm.Send([u_right, MPI.DOUBLE], dest = 0, tag = 43)
                comm.Send([v2, MPI.DOUBLE], dest = 0, tag = 888)
            if rank == 0:
                comm.Recv([u_left, MPI.DOUBLE], source = 1, tag = 33)
                comm.Recv([u_right, MPI.DOUBLE], source = 2, tag = 43)
                comm.Recv([rooml, MPI.DOUBLE], source = 1, tag = 666)
                comm.Recv([roomr, MPI.DOUBLE], source = 2, tag = 888)
                if i == iteration - 1:
                    room0 = v0
                    room0 = room0.reshape((2 * num - 1, num - 1))
                    rooml = rooml.reshape((num, num))
                    roomr = np.fliplr(roomr.reshape((num - 1, num))) 
                    fig = plt.figure(figsize=(6, 4))                
                    ax = fig.add_axes([0.1,0.1,0.8,0.8]) 
                    ax.set_xticks([0,1,2,3])
                    ax.set_yticks([0,1,2])
                    ax.set_title('Temperature Distribution' + ' %s' %(i))
                    size = [0.101,0.1,0.267,0.4]  #left,bottom,width,height 
                    _plot(rooml, fig, size, origin = 'lower')
                    size = [0.368,0.1,0.267,0.8]  #left,bottom,width,height 
                    _plot(room0, fig, size, origin = 'upper')
                    size = [0.635,0.5,0.333,0.4]  #left,bottom,width,height 
                    _plot(roomr, fig, size, origin = 'lower')
                    plt.colorbar(orientation='vertical')
                    plt.show()