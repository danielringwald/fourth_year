from mpi4py import MPI

comm = MPI.COMM_WORLD

print('Hello World: process', comm.Get_rank(), ' out of', comm.Get_size(), ' is report for duty!')
