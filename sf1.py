from mpi4py import MPI
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()




import numpy as np
from math import *
import random
from numpy.random import rand
import pickle
import sf

Trank = 39
xd = 0.0
Lx = Ly = Lz = L = 10                                                               #LENGTH
Ns = 16*Lx*Ly*Lz          
unit_cell = [None] * 16
unit_cell[0] = np.array([0,0,0])
unit_cell[1] = np.array([1,1,0])
unit_cell[2] = np.array([0,1,1])
unit_cell[3] = np.array([1,0,1])
unit_cell[4] = np.array([2,2,0])
unit_cell[5] = np.array([3,3,0])
unit_cell[6] = np.array([2,3,1])
unit_cell[7] = np.array([3,2,1])
unit_cell[8] = np.array([0,2,2])                                                        #UNIT CELL
unit_cell[9] = np.array([1,3,2])
unit_cell[10] = np.array([0,3,3])
unit_cell[11] = np.array([1,2,3])
unit_cell[12]= np.array([2,0,2])
unit_cell[13] = np.array([3,1,2])
unit_cell[14] = np.array([2,1,3])
unit_cell[15] = np.array([3,0,3])


coordinates = []
for z in range(0,int(Lx)):
    for y in range(0,int(Ly)):
        for x in range(0,int(Lz)):
            coordinates.append(unit_cell + np.array([4*x,4*y,4*z]))
                                                                                        #COORDINATES of SPINS in LATTICE
coordinates = np.array(coordinates)
coordinates_temp = coordinates
coordinates = coordinates.reshape(coordinates.shape[0]*coordinates.shape[1],3 )
coordinates_list = coordinates.tolist()

matrix = np.zeros((4*int(L),4*int(L),4*int(L)))
for i in range(int(Ns)):
    matrix[coordinates[i][0],coordinates[i][1],coordinates[i][2]] = i                   #SPIN INDEXING

def random_spin():
    temp = np.array([np.random.uniform(-1,1), np.random.uniform(-1,1),np.random.uniform(-1,1)],order='F')
    spin = temp/(2*np.linalg.norm(temp))
    return spin                                                                         #GENERATING a SPIN

def initialstate(Ns):
    ''' generates a random spin configuration for initial condition'''
    state = np.array([random_spin() for i in range(int(Ns))],order='F')
    return state     


h  =  np.arange(-4*L,4*L+1, 1)
k = np.arange(-4*L,4*L+1, 1)
Spinsall = pickle.load( open( "SpinsDIS20_%s_xd%s.p" %(Trank,xd), "rb" ) )
Spins = Spinsall[my_rank*10:(my_rank+1)*10]

#Spins = Spins[0:1]


SfT = np.zeros((len(h),len(k)))  
S_q = np.zeros((len(h),len(k))) 


import time
start = time.process_time()


for i in range(len(Spins)): 
    config = Spins[i] 
    for i1 in range(len(h)): 
        for i2 in range(len(k)): 
            q = ((2.0*pi)/(4.0*L))*np.array([h[i1],h[i1],k[i2]]) 
            S_q[i1,i2] =  sf.structurefactor(config,coordinates,q) 
    SfT = SfT + S_q 
     
SfT = SfT/len(Spins) 
 
print(time.process_time() - start,my_rank)

Str_fact =  np.array(p*[S_q])
if my_rank!=0:
        comm.send(SfT, dest=0,tag=1)
else :
        
        Str_fact[0] = SfT
        for procid in range(1,p):
            Str_fact[procid] = comm.recv(source=procid,tag=1)
    
    
        Str_fact = np.mean(Str_fact,axis=0)  
        pickle.dump(Str_fact, open( "STRF_%s_%s.p" %(Trank,xd), "wb" ) )
        
        


#acceptance_ratio = (accept/count)*100
#MPI.Finalize()




from mpi4py import MPI
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()




import numpy as np
from math import *
import random
from numpy.random import rand
import pickle
import sf

Trank = 49
xd = 0.0
Lx = Ly = Lz = L = 10                                                               #LENGTH
Ns = 16*Lx*Ly*Lz          
unit_cell = [None] * 16
unit_cell[0] = np.array([0,0,0])
unit_cell[1] = np.array([1,1,0])
unit_cell[2] = np.array([0,1,1])
unit_cell[3] = np.array([1,0,1])
unit_cell[4] = np.array([2,2,0])
unit_cell[5] = np.array([3,3,0])
unit_cell[6] = np.array([2,3,1])
unit_cell[7] = np.array([3,2,1])
unit_cell[8] = np.array([0,2,2])                                                        #UNIT CELL
unit_cell[9] = np.array([1,3,2])
unit_cell[10] = np.array([0,3,3])
unit_cell[11] = np.array([1,2,3])
unit_cell[12]= np.array([2,0,2])
unit_cell[13] = np.array([3,1,2])
unit_cell[14] = np.array([2,1,3])
unit_cell[15] = np.array([3,0,3])


coordinates = []
for z in range(0,int(Lx)):
    for y in range(0,int(Ly)):
        for x in range(0,int(Lz)):
            coordinates.append(unit_cell + np.array([4*x,4*y,4*z]))
                                                                                        #COORDINATES of SPINS in LATTICE
coordinates = np.array(coordinates)
coordinates_temp = coordinates
coordinates = coordinates.reshape(coordinates.shape[0]*coordinates.shape[1],3 )
coordinates_list = coordinates.tolist()

matrix = np.zeros((4*int(L),4*int(L),4*int(L)))
for i in range(int(Ns)):
    matrix[coordinates[i][0],coordinates[i][1],coordinates[i][2]] = i                   #SPIN INDEXING

def random_spin():
    temp = np.array([np.random.uniform(-1,1), np.random.uniform(-1,1),np.random.uniform(-1,1)],order='F')
    spin = temp/(2*np.linalg.norm(temp))
    return spin                                                                         #GENERATING a SPIN

def initialstate(Ns):
    ''' generates a random spin configuration for initial condition'''
    state = np.array([random_spin() for i in range(int(Ns))],order='F')
    return state     


h  =  np.arange(-4*L,4*L+1, 1)
k = np.arange(-4*L,4*L+1, 1)
Spinsall = pickle.load( open( "SpinsDIS20_%s_xd%s.p" %(Trank,xd), "rb" ) )
Spins = Spinsall[my_rank*10:(my_rank+1)*10]

#Spins = Spins[0:1]


SfT = np.zeros((len(h),len(k)))  
S_q = np.zeros((len(h),len(k))) 


import time
start = time.process_time()


for i in range(len(Spins)): 
    config = Spins[i] 
    for i1 in range(len(h)): 
        for i2 in range(len(k)): 
            q = ((2.0*pi)/(4.0*L))*np.array([h[i1],h[i1],k[i2]]) 
            S_q[i1,i2] =  sf.structurefactor(config,coordinates,q) 
    SfT = SfT + S_q 
     
SfT = SfT/len(Spins) 
 
print(time.process_time() - start,my_rank)

Str_fact =  np.array(p*[S_q])
if my_rank!=0:
        comm.send(SfT, dest=0,tag=1)
else :
        
        Str_fact[0] = SfT
        for procid in range(1,p):
            Str_fact[procid] = comm.recv(source=procid,tag=1)
    
    
        Str_fact = np.mean(Str_fact,axis=0)  
        pickle.dump(Str_fact, open( "STRF_%s_%s.p" %(Trank,xd), "wb" ) )
        
        


#acceptance_ratio = (accept/count)*100
#MPI.Finalize()




from mpi4py import MPI
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()




import numpy as np
from math import *
import random
from numpy.random import rand
import pickle
import sf

Trank = 59
xd = 0.0
Lx = Ly = Lz = L = 10                                                               #LENGTH
Ns = 16*Lx*Ly*Lz          
unit_cell = [None] * 16
unit_cell[0] = np.array([0,0,0])
unit_cell[1] = np.array([1,1,0])
unit_cell[2] = np.array([0,1,1])
unit_cell[3] = np.array([1,0,1])
unit_cell[4] = np.array([2,2,0])
unit_cell[5] = np.array([3,3,0])
unit_cell[6] = np.array([2,3,1])
unit_cell[7] = np.array([3,2,1])
unit_cell[8] = np.array([0,2,2])                                                        #UNIT CELL
unit_cell[9] = np.array([1,3,2])
unit_cell[10] = np.array([0,3,3])
unit_cell[11] = np.array([1,2,3])
unit_cell[12]= np.array([2,0,2])
unit_cell[13] = np.array([3,1,2])
unit_cell[14] = np.array([2,1,3])
unit_cell[15] = np.array([3,0,3])


coordinates = []
for z in range(0,int(Lx)):
    for y in range(0,int(Ly)):
        for x in range(0,int(Lz)):
            coordinates.append(unit_cell + np.array([4*x,4*y,4*z]))
                                                                                        #COORDINATES of SPINS in LATTICE
coordinates = np.array(coordinates)
coordinates_temp = coordinates
coordinates = coordinates.reshape(coordinates.shape[0]*coordinates.shape[1],3 )
coordinates_list = coordinates.tolist()

matrix = np.zeros((4*int(L),4*int(L),4*int(L)))
for i in range(int(Ns)):
    matrix[coordinates[i][0],coordinates[i][1],coordinates[i][2]] = i                   #SPIN INDEXING

def random_spin():
    temp = np.array([np.random.uniform(-1,1), np.random.uniform(-1,1),np.random.uniform(-1,1)],order='F')
    spin = temp/(2*np.linalg.norm(temp))
    return spin                                                                         #GENERATING a SPIN

def initialstate(Ns):
    ''' generates a random spin configuration for initial condition'''
    state = np.array([random_spin() for i in range(int(Ns))],order='F')
    return state     


h  =  np.arange(-4*L,4*L+1, 1)
k = np.arange(-4*L,4*L+1, 1)
Spinsall = pickle.load( open( "SpinsDIS20_%s_xd%s.p" %(Trank,xd), "rb" ) )
Spins = Spinsall[my_rank*10:(my_rank+1)*10]

#Spins = Spins[0:1]


SfT = np.zeros((len(h),len(k)))  
S_q = np.zeros((len(h),len(k))) 


import time
start = time.process_time()


for i in range(len(Spins)): 
    config = Spins[i] 
    for i1 in range(len(h)): 
        for i2 in range(len(k)): 
            q = ((2.0*pi)/(4.0*L))*np.array([h[i1],h[i1],k[i2]]) 
            S_q[i1,i2] =  sf.structurefactor(config,coordinates,q) 
    SfT = SfT + S_q 
     
SfT = SfT/len(Spins) 
 
print(time.process_time() - start,my_rank)

Str_fact =  np.array(p*[S_q])
if my_rank!=0:
        comm.send(SfT, dest=0,tag=1)
else :
        
        Str_fact[0] = SfT
        for procid in range(1,p):
            Str_fact[procid] = comm.recv(source=procid,tag=1)
    
    
        Str_fact = np.mean(Str_fact,axis=0)  
        pickle.dump(Str_fact, open( "STRF_%s_%s.p" %(Trank,xd), "wb" ) )
        
        


#acceptance_ratio = (accept/count)*100
#MPI.Finalize()
