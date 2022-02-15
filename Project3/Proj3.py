from scipy import*
from scipy.linalg import*
from numpy import*
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpi4py import MPI


class Room:

    def __init__(self, deltax, gamma_normal=15, gamma_H=40, gamma_WF=5, omega=0.8):
        self.deltax = deltax
        self.gamma_normal = gamma_normal
        self.gamma_H = gamma_H
        self.gamma_WF = gamma_WF
        grid_x = linspace(0, 1, int(1/deltax + 1))
        grid_y = linspace(0, 2, int(2/deltax + 1))
        self.grid_room2 = [grid_x, grid_y]
        self.grid_room13 = [grid_x, grid_x]
        self.Nx = int(1/self.deltax - 1)
        self.Ny = int(2/self.deltax - 1)

    def __call__(self, room_nbr, test_neumann=None):
        """
        Plots and returns solution for specified room

        Parameters
        ----------
        room_nbr : INT, 
        test_neumann : array-like, optional
            To plot room 1 and 3 a Neumann condition has to be added. The default is None.

        Raises
        ------
        Exception
         

        Returns
        -------
        room_solution : ndarray
            Solution for temperatures in specified room.

        """
        if not type(room_nbr) == type(None):
            if room_nbr == 2:
                room_solution, condition = self.solve_room(room_nbr)
                
            if not type(test_neumann) == type(None):
                if room_nbr == 1 or room_nbr == 3:
                    room_solution, condition = self.solve_room(room_nbr, test_neumann)
            self.plot(room_solution, room_nbr, condition)
            plt.show()
            return room_solution       
        else:
            raise Exception
           

    def neumann_dirichlet_method(self, begin_dirichlet=None):
        '''         
        Returns the solution u for the linear system Au = b using direchlet conditions.

            Parameters:
                A (matrix): Matrix A
                b (vector): Bouandary values for the room

            Returns:
                u (vector): Temperature distribution
        '''

        sol2, cond2 = self.solve_room(2)
        sol1, cond1 = self.solve_room(1, cond2[0])
        sol3, cond3 = self.solve_room(3, cond2[1])
        sol2, cond2 = self.solve_room(2, [cond1, cond3])
        
        return [sol1, cond1], [sol2,cond2], [sol3, cond3]

    def solve_room(self, room_nbr, condition=None):
        Ny = self.Ny
        Nx = self.Nx

        if type(condition) == type(None):
            condition = [[0]*int(floor((Ny/2)))]*2

        if room_nbr == 2:
            A = self.create_matrix_room(room_nbr)
            b = self.rhs(room_nbr, condition)
            solution = solve(A, b)
            
            neumann_condition1 = (condition[0] - 
                            solution[:int(floor((Ny/2)))]) / self.deltax 
            neumann_condition3 = (condition[1] - 
                            solution[-int(floor((Ny/2))):]) / self.deltax
            
            new_condition = [neumann_condition1, neumann_condition3]

        elif room_nbr == 1 or room_nbr == 3:
            A = self.create_matrix_room(room_nbr)
            b = self.rhs(room_nbr, condition)
            solution = solve(A, b)
            if room_nbr == 1:
                new_condition = solution[-Nx:]
            else:
                new_condition = solution[:Nx]
        else:
            raise ValueError

        return solution, new_condition

    def create_matrix_room(self, room_nbr):
        '''
        Returns the matrix for the linear system Au = b, where b is the bouandary temperatures for the room

            Parameters:
                    room_nbr (int): Index for the chosen room, must be 1,2 or 3.

            Returns:
                    room_matrix (numpy.array): Matrix for the chosen room
        '''
        if room_nbr > 3 and room_nbr < 1:
            raise ValueError

        Nx = self.Nx
        Ny = self.Ny

        if room_nbr == 1 or room_nbr == 3:
            n = Nx**2
            
            neu_length = n + Nx
            neu_main_diag = diag([-3]*(neu_length))
            sub_diag = diag([1]*(neu_length-1), -1)
            sup_diag = diag([1]*(neu_length-1), 1)
            subsub_diag = diag([1]*(n), -Nx)
            supsup_diag = diag([1]*(n), Nx)

            neumann_cond_matrix = neu_main_diag + sub_diag + \
                sup_diag + subsub_diag + supsup_diag

            main_diag = diag([-4]*n)
            sub_diag = diag([1]*(n-1), -1)
            sup_diag = diag([1]*(n-1), 1)
            subsub_diag = diag([1]*(n - Nx), -Nx)
            supsup_diag = diag([1]*(n - Nx), Nx)

            inner_points_matrix = main_diag + sub_diag + \
                sup_diag + supsup_diag + subsub_diag
            neumann_cond_matrix[:n, :n] = inner_points_matrix

            for i in range(1, Nx+1):
                neumann_cond_matrix[i*Nx, i*Nx - 1] = 0
                neumann_cond_matrix[i*Nx - 1, i*Nx] = 0
            if room_nbr == 3:
                return neumann_cond_matrix[::-1, ::-1] / self.deltax**2
            else:
                return neumann_cond_matrix / self.deltax**2

        elif room_nbr == 2:
            n = int(Nx * Ny)

            main_diag = diag([-4]*n)
            sub_diag = diag([1]*(n-1), -1)
            sup_diag = diag([1]*(n-1), 1)
            subsub_diag = diag([1]*(n - Ny), -Ny)
            supsup_diag = diag([1]*(n - Ny), Ny)

            room_matrix = main_diag + sub_diag + sup_diag + supsup_diag + subsub_diag

            for i in range(1, Nx):
                room_matrix[i*Ny, i*Ny - 1] = 0
                room_matrix[i*Ny - 1, i*Ny] = 0

        return room_matrix/self.deltax**2

    def rhs(self, room_nbr, condition=None):
       
        Nx = self.Nx
        Ny = self.Ny
        if type(condition) == type(None):
            condition = [0] * int(1 / self.deltax - 1)
        if room_nbr == 3 or room_nbr == 1:
            b = zeros((Nx*(Nx + 1),))

            # adds left Dirichlet condition
            b[0: Nx] = self.gamma_H
            #add Neumann condition
            b[-Nx:] = condition 
            b[-Nx:] = b[-Nx:] * self.deltax**2
            
            # adds top and bottom Dirichlet condition
            for i in range(Nx+1):
                b[i*Nx] = b[i*Nx] + self.gamma_normal
                b[Nx - 1 + i*(Nx)] = b[Nx - 1 + i*(Nx)] + self.gamma_normal
                
            if room_nbr == 3:
                return -b[::-1] / self.deltax**2
            else:
                return -b / self.deltax**2

        elif room_nbr == 2:
            b = zeros((Nx*Ny,))
            b[int(floor((Ny/2))): Ny] = self.gamma_normal
            b[:int(floor((Ny/2)))] = condition[0]
            b[-int(floor((Ny/2))):] = condition[1]
            b[-Ny:-int(floor((Ny/2)))] = self.gamma_normal

            for i in range(Nx):
                b[i*Ny] = b[i*Ny] + self.gamma_WF
                b[Ny - 1 + i*(Ny)] = b[Ny - 1 + i*(Ny)] + self.gamma_H

            return -b / self.deltax**2
        else:
            raise ValueError

    def plot(self, room_sol, room_nbr, cond=None, plot_figure=True):
        Nx = self.Nx
        Ny = self.Ny

        if room_nbr == 2:
            
            temps = zeros((Ny + 2, Nx + 2))
            temps[1:-1, 1:-1] = reshape(room_sol, (Nx, Ny)).T
            temps[0, :] = self.gamma_WF
            temps[-1, :] = self.gamma_H
            
            normal_temp = zeros(Ny)
            normal_temp[:int(floor(len(normal_temp)/2)) +
                        1] = self.gamma_normal
            if type(cond) != type(None):
                normal_temp[int(floor(len(normal_temp)/2)) +1:] = cond[0]
            temps[1:-1, 0] = normal_temp[::-1]
           
            if type(cond) != type(None):
                normal_temp[int(floor(len(normal_temp)/2)) +1:] = cond[1]
            temps[1:-1, -1] = normal_temp
            [T, X] = meshgrid(self.grid_room2[0], self.grid_room2[1])
            #print(temps)
        elif room_nbr == 1 or room_nbr == 3:
            temps = zeros((Nx + 2, Nx + 2))
            Nx = Nx + 1
            if room_nbr == 3:
                room_sol = room_sol[::-1]
            temps[1:-1, 1:] = reshape(room_sol, (Nx, Nx-1)).T
            temps[0, :] = self.gamma_normal
            temps[-1, :] = self.gamma_normal
            temps[:, 0] = (Nx + 1) * [self.gamma_H]
            [T, X] = meshgrid(self.grid_room13[0], self.grid_room13[1])
            
        ax = plt.axes(projection='3d')
        ax.plot_surface(T, X, temps, cmap=cm.plasma)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('Temperature')
        name = 'Room: ' + str(room_nbr)
        ax.set_title(name)

if __name__ == '__main__':
    deltax = 1 / 3
    Nx = int(1 / deltax - 1)
    Ny = int(2 / deltax - 1)
    omega = 0.8
    room = Room(deltax)
    comm = MPI.Comm.Clone(MPI.COMM_WORLD)
    rank = comm.Get_rank()
    iteration = 10
   
    #First solve of room 2
    dirichlet_r1 = array([15.]*Nx)
    dirichlet_r3 = array([15.]*Nx)
    temp_init_r2, cond_init = room.solve_room(2, [dirichlet_r1, dirichlet_r3])
    gamma_1_cond = cond_init[0]
    gamma_3_cond = cond_init[1]
    
    temp_init_r1, cond_init_r1 = room.solve_room(1, gamma_1_cond)
    temp_init_r3, cond_init_r3 = room.solve_room(3, gamma_3_cond)
    
    room1 = temp_init_r1
    room2 = temp_init_r2
    room3 = temp_init_r3
    
    
    for i in range(iteration):
      
        if rank == 0:
            comm.Send([gamma_1_cond, MPI.DOUBLE], dest=1)
            comm.Send([gamma_3_cond, MPI.DOUBLE], dest=2)
        
        #Room 1
        if rank == 1:
            comm.Recv([gamma_1_cond, MPI.DOUBLE], source=0)
            temps_r1, dirichlet_r1 = room.solve_room(1, gamma_1_cond)
            room1 = room1*(1-omega) + temps_r1*omega
            comm.Send([dirichlet_r1, MPI.DOUBLE], dest=0)
        
        #Room 3
        if rank == 2:
            comm.Recv([gamma_3_cond, MPI.DOUBLE], source=0)
            temps_r3, dirichlet_r3 = room.solve_room(3, gamma_3_cond)
            room3 = room3*(1-omega) + temps_r3*omega
            comm.Send([dirichlet_r3, MPI.DOUBLE], dest=0)
        
        #Room 2
        if rank == 0: 
            comm.Recv([dirichlet_r1, MPI.DOUBLE], source=1)
            comm.Recv([dirichlet_r3, MPI.DOUBLE], source=2)
            temps_r2, neumann_conds = room.solve_room(2, [dirichlet_r1, dirichlet_r3])
            #Relaxation
            room2 = room2 * (1 - omega) + temps_r2 * omega
            gamma_1_cond = neumann_conds[0]
            gamma_3_cond = neumann_conds[1]
            
            if i == iteration-1:
                
                x_dim_big = int(3 / deltax - 1) + 2
                y_dim_big = Ny+2
                big_room = zeros((y_dim_big, x_dim_big))
                
                top = array([0]*(Nx+1) + [room.gamma_H]*(Nx +2)+ [room.gamma_normal]*(Nx+1))
                bottom = array([room.gamma_normal]*(Nx+1) + [room.gamma_WF]*(Nx+2) + [0]*(Nx+1))
                left = array([room.gamma_H]*(int(floor(y_dim_big/2))+1) + [0]*(int(floor(y_dim_big/2))))
                right = array([0]*(int(floor(y_dim_big/2))) + [room.gamma_H]*(int(floor(y_dim_big/2))+1))
                middle_vert = array([room.gamma_normal]*(Nx+1))
                middle_hori = array([room.gamma_normal]*Nx)
                
                #Insert the boundaries
                big_room[0,:] = bottom
                big_room[-1,:] = top
                big_room[:,0] = left
                big_room[:,-1] = right
                print(big_room)
                big_room[(Nx+1):-1, Nx+1] = middle_vert
                big_room[1:(Nx+2), 2*Nx+2] = middle_vert
                big_room[Nx+1, 1:Nx+1] = middle_hori
                big_room[Nx+1, -(Nx+1):-1] = middle_hori
                #Insert solutions
                big_room[1:Nx+1, 1:Nx+2] = reshape(room1, (Nx+1,Nx)).T
                big_room[-(Nx+1):-1, -(Nx+2):-1] = reshape(room3, (Nx+1,Nx)).T
                big_room[1:-1, Nx+2:-(Nx+2)] = reshape(room2,(Nx,Ny)).T
                
                #Plot
                print(big_room)
                big_room[big_room==0] = nan
                plt.imshow(big_room, cmap='cool', origin = 'lower', interpolation = 'gaussian', extent=[0,3,0,2])
                plt.xlabel('x')
                plt.ylabel('y')
                plt.title('Kitchen area')
                plt.colorbar()
                               
                plt.show()
        
        
    
        
    #https://www.hemhyra.se/nyheter/sa-varmt-har-du-ratt-till-att-ha-hemma/
    
    
    
    #test_neumann = array([-5.315, -7.1875]) / deltax #[2] * int(1 / deltax - 1)
    #Test room 1
    # neu_matrix = room.create_matrix_room(1)
    # print('Room 1')
    # print(neu_matrix)
    # b = room.rhs(1, test_neumann)
    # print(b)
    # plt.figure()
    # sol1 = solve(neu_matrix,b)
    # room.plot(sol1, 1)
    
    #Test room 3
    # neu_matrix = room.create_matrix_room(3)
    # b = room.rhs(3, test_neumann)
    # print('Room 3')
    # print(neu_matrix)
    # print(b)
    # plt.figure()
    # sol3 = solve(neu_matrix, b)
    # room.plot(sol3, 3)
    
    #Test room 2
    # matrix = room.create_matrix_room(2)
    # b = room.rhs(2)
    # plt.figure()
    # sol2 = solve(matrix,b)
    # room.plot(sol2, 2)
    
    #Test solve_room
    # sol = room.solve_room(2)
    # print(sol[1])
    
    #Test iteration without method
    # sol = room.solve_room(2)
    # print('Room 2 condition')
    # print(sol[1])
    # room.plot(sol[0],2)
    # sol = room.solve_room(1, sol[1][0])
    # print('\nRoom 1 solution')
    # print(sol)
    # print('Dirichlet conditions')
    # print(sol[0][-room.Nx:])
    
    #Test one iteration
    # room1, room2, room3 = room.neumann_dirichlet_method()
    # plt.figure(1)
    # room.plot(room1[0], 1)
    # plt.figure(2)
    # room.plot(room2[0], 2, [room1[1], room3[1]])
    # plt.figure(3)
    # room.plot(room3[0], 3)
    
    #Show for all
    #plt.show()
