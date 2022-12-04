
from mpi4py import MPI

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

import numpy as np
from math import *
import random
from numpy.random import rand

import neighb
import pickle


import time

start = time.process_time()


Trank = 4
Lx = Ly = Lz = L = 4
Ns = 16 * Lx * Ly * Lz
xd = 0.3

unit_cell = [None] * 16
unit_cell[0] = np.array([0, 0, 0])
unit_cell[1] = np.array([1, 1, 0])
unit_cell[2] = np.array([0, 1, 1])
unit_cell[3] = np.array([1, 0, 1])
unit_cell[4] = np.array([2, 2, 0])
unit_cell[5] = np.array([3, 3, 0])
unit_cell[6] = np.array([2, 3, 1])
unit_cell[7] = np.array([3, 2, 1])
unit_cell[8] = np.array([0, 2, 2])
unit_cell[9] = np.array([1, 3, 2])
unit_cell[10] = np.array([0, 3, 3])
unit_cell[11] = np.array([1, 2, 3])
unit_cell[12] = np.array([2, 0, 2])
unit_cell[13] = np.array([3, 1, 2])
unit_cell[14] = np.array([2, 1, 3])
unit_cell[15] = np.array([3, 0, 3])

coordinates = []
for z in range(0, int(Lx)):
    for y in range(0, int(Ly)):
        for x in range(0, int(Lz)):
            coordinates.append(unit_cell + np.array([4 * x, 4 * y, 4 * z]))
coordinates = np.array(coordinates)
coordinates_temp = coordinates
coordinates = coordinates.reshape(coordinates.shape[0] * coordinates.shape[1], 3)
coordinates_list = coordinates.tolist()

matrix = np.zeros((4 * int(L), 4 * int(L), 4 * int(L)))
for i in range(int(Ns)):
    matrix[coordinates[i][0], coordinates[i][1], coordinates[i][2]] = i


def random_spin():
    temp = np.array(
        [np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)],
        order="F",
    )
    spin = temp / (2 * np.linalg.norm(temp))
    return spin


def initialstate(Ns):
    """ generates a random spin configuration for initial condition"""
    state = np.array([random_spin() for i in range(int(Ns))], order="F")
    return state


neigh = neighb.neighbour(coordinates, matrix)
neigh = neigh.astype(int)




Spins = pickle.load( open( "Spins_%s_xd%s.p" %(Trank,xd), "rb" ) )


ME_PDF = np.zeros([len(Spins),int(Ns/2),2])


import time
start = time.process_time()

for j in range(len(Spins)):
    if my_rank == 0:
       print(j)
    mepdf = np.zeros([int(Ns/2),2])

    config = Spins[j]
    count =0
    for i in range(0,Ns,4):

#        print(i,count)
        mepdf[count][0] = (-2.0*config[i,0] + config[i,1]+ config[i,2] \
                        -2.0*config[i+2,0] - config[i+2,1]- config[i+2,2] \
                        +2.0*config[i+3,0] + config[i+3,1]- config[i+3,2] \
                        +2.0*config[i+1,0] - config[i+1,1]+ config[i+1,2])/(2.0*sqrt(6.0))
                    
        mepdf[count][1] = ( -config[i,1] + config[i,2]+ config[i+2,1]-config[i+2,2] \
                         -config[i+3,1]- config[i+3,2] +config[i+1,1]+config[i+1,2] ) /(2.0*sqrt(2.0))  

        
        mepdf[count+int(Ns/4)][0] = (-2.0*config[i,0] + config[i,1]+ config[i,2] \
                        -2.0*config[neigh[i,4],0] - config[neigh[i,4],1]- config[neigh[i,4],2] \
                        +2.0*config[neigh[i,5],0] + config[neigh[i,5],1]- config[neigh[i,5],2] \
                        +2.0*config[neigh[i,3],0] - config[neigh[i,3],1]+ config[neigh[i,3],2])/(2.0*sqrt(6.0))
                    
        mepdf[count+int(Ns/4)][1] = ( -config[i,1] + config[i,2]+ config[neigh[i,4],1]-config[neigh[i,4],2] \
                         -config[neigh[i,5],1]- config[neigh[i,5],2] +config[neigh[i,3],1]+config[neigh[i,3],2] ) /(2.0*sqrt(2.0))  
        
        
        
          
        count = count + 1

    
        





    ME_PDF[j] = mepdf
        
pickle.dump(ME_PDF, open("ME_PDF_xd%s_%s.p" % (xd,Trank), "wb"))

print(time.process_time() - start,count)


