from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
iteration = 10


for i in range(iteration):
    if rank == 1:
        if i == 0:
            data0 = np.array([0]*3, dtype='d')
            data2 = np.array([2]*3, dtype='d')
            comm.Send([data0, MPI.DOUBLE], dest=0)
            comm.Send([data2, MPI.DOUBLE], dest=2)
            print("Sent data to process 0")
        else:
            data0 = np.empty(10, dtype='d')
            data2 = np.empty(10, dtype='d')
            comm.Recv(data0, source=0)
            print("Recieved data from process 0")
            comm.Recv(data2, source=2)
            print("Recieved data from process 2")
            data0 = data0 + 1
            comm.Send([data0, MPI.DOUBLE], dest=0)
            print("Sent data: ", data0)
            data2 = data2 - 1
            comm.Send([data2, MPI.DOUBLE], dest=2)
            print("Sent data: ", data2)
    if rank == 0:
         data0 = np.empty(10, dtype='d')
         comm.Recv(data0, source=1)
         print("Data recieved: ", data0, ", from process 1")
         comm.Send([data0*2**2, MPI.DOUBLE], dest=1)
    if rank == 2:
        data2 = np.empty(10, dtype='d')
        comm.Recv(data2, source=1)
        print("Data recieved: ", data2, ", from process 1")
        comm.Send([data2*i**2, MPI.DOUBLE], dest=1)
