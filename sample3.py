from mpi4py import MPI

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
p = comm.Get_size()

import numpy as np
from math import *
import random
from numpy.random import rand

import monte
import order
import over
import neighb
import order
import pickle
import inter

import time

start = time.process_time()

num = 3
Lx = Ly = Lz = L = 4
Ns = 16 * Lx * Ly * Lz
xd = 1.2
w = 10
nt = 100
T = np.linspace(0.01, 0.5, nt)
gap = 50
mcSteps = 10**6
eqSteps = int(0.2 * mcSteps)
mcSteps_min = int(0.2 * mcSteps)

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
neigh_T = pickle.load(open("neigh_T_%s.p" % L, "rb"))

fac = 11.604
a = fac * np.array([0.07, 0.08, -0.11, 0.04])
b = fac * np.array([0.11, -0.06, -0.10, -0.003])
ab = (a + b) / 2.0

if my_rank == 0:
    J_combined = np.array([np.array([np.zeros(4)] * 6)] * Ns)
    J_total = np.array([J_combined] * w)
    for iw in range(0, w, 2):
        arr = np.empty(Ns)
        for i in range(Ns):
            if random.random() < xd / 2.0:
                arr[i] = 0
            else:
                arr[i] = 1
        for i in range(Ns):
            for j in range(6):
                if np.linalg.norm(J_combined[i][j]) == 0.0:
                    wtr = neigh_T[i][j]
                    arr_sum = arr[int(wtr[0])] + arr[int(wtr[1])]
                    if arr_sum == 0:
                        J_total[iw][i][j] = a
                        J_total[iw][int(neigh[i][j])][
                            np.where(neigh[int(neigh[i][j])] == i)[0][0]
                        ] = a
                    elif arr_sum == 2:
                        J_total[iw][i][j] = b
                        J_total[iw][int(neigh[i][j])][
                            np.where(neigh[int(neigh[i][j])] == i)[0][0]
                        ] = b
                    else:
                        J_total[iw][i][j] = ab
                        J_total[iw][int(neigh[i][j])][
                            np.where(neigh[int(neigh[i][j])] == i)[0][0]
                        ] = ab
        J_total[iw + 1] = J_total[iw]
    for procid in range(1, p):
        comm.send(J_total, dest=procid, tag=420)
else:
    J_total = comm.recv(source=0, tag=420)
Et, Mt, Ct, Xt = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
MA2t, MEt, MT1Bt, MT2t, METHETAt, X_MEt, X_MT2t = (
    np.zeros(nt),
    np.zeros(nt),
    np.zeros(nt),
    np.zeros(nt),
    np.zeros(nt),
    np.zeros(nt),
    np.zeros(nt),
)
X_qt, g4t, g2t, Gt = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
Ew = Mw = Cw = Xw = MA2w = MEw = MT1Bw = MT2w = METHETAw = X_MEw = X_MT2w = 0.0
Q2 = Q4 = Q2_sq = 0.0

iT = 1.0 / T[my_rank]
iT2 = iT * iT
w_count = 0
T_0 = 1.0
slope = (T[my_rank] - T_0) / mcSteps_min


for iw in range(w):
    config = initialstate(Ns)
    J_combined = J_total[iw]
    Jij_matrix = inter.interaction(neigh, J_combined, L, Ns)
    E1 = M1 = E2 = M2 = ma2 = me = mt1b = mt2 = metheta = me_2 = mt2_2 = 0.0

    if iw % 2 == 0:
        S = np.array([np.empty([Ns,3])]*int(mcSteps/gap))
    if iw % 2 == 1:
        q_p2 = 0.0
        q_p4 = 0.0
    for i in range(mcSteps_min):
        T_anneal = slope * i + T_0
        iT_anneal = 1 / T_anneal
        monte.montecarlo(config, coordinates, matrix, neigh, Jij_matrix, iT_anneal)
    for i in range(eqSteps):
        monte.montecarlo(config, coordinates, matrix, neigh, Jij_matrix, iT)
    for i in range(mcSteps):
        if i % 1000 == 0 and my_rank == 0:
            print(i)
        monte.montecarlo(config, coordinates, matrix, neigh, Jij_matrix, iT)
        over.overrelaxation(config, coordinates, matrix, neigh, Jij_matrix, iT)
        EM = order.calcenergy(config, coordinates, matrix, neigh, Jij_matrix)

        if i % 20 == 10 and i != 0:
            if np.mod(my_rank, 2) == 0:
                comm.Send(EM, dest=my_rank + 1, tag=1)
                comm.Send(config, dest=my_rank + 1, tag=2)
                var = comm.recv(source=my_rank + 1, tag=111)
                if var == 1:
                    EM = np.empty(7)
                    comm.Recv(EM, source=my_rank + 1, tag=3)
                    config = np.empty([Ns, 3], order="F")
                    comm.Recv(config, source=my_rank + 1, tag=4)
            else:
                EM_prev = np.empty(7)
                comm.Recv(EM_prev, source=my_rank - 1, tag=1)
                config_prev = np.empty([Ns, 3], order="F")
                comm.Recv(config_prev, source=my_rank - 1, tag=2)
                delta = (1 / T[my_rank - 1] - 1 / T[my_rank]) * (EM_prev[0] - EM[0])
                if random.random() < np.exp(delta):
                    var = 1
                    comm.send(var, dest=my_rank - 1, tag=111)
                    comm.Send(EM, dest=my_rank - 1, tag=3)
                    comm.Send(config, dest=my_rank - 1, tag=4)
                    EM = EM_prev
                    config = config_prev
                else:
                    var = 0
                    comm.send(var, dest=my_rank - 1, tag=111)
        if i % 20 == 0 and i != 0:
            if np.mod(my_rank, 2) == 1 and my_rank < p - 1:
                comm.Send(EM, dest=my_rank + 1, tag=21)
                comm.Send(config, dest=my_rank + 1, tag=22)
                var = comm.recv(source=my_rank + 1, tag=222)
                if var == 1:
                    EM = np.empty(7)
                    comm.Recv(EM, source=my_rank + 1, tag=23)
                    config = np.empty([Ns, 3], order="F")
                    comm.Recv(config, source=my_rank + 1, tag=24)
            elif np.mod(my_rank, 2) == 0 and my_rank < p - 1 and my_rank != 0:
                EM_prev = np.empty(7)
                comm.Recv(EM_prev, source=my_rank - 1, tag=21)
                config_prev = np.empty([Ns, 3], order="F")
                comm.Recv(config_prev, source=my_rank - 1, tag=22)
                delta = (1 / T[my_rank - 1] - 1 / T[my_rank]) * (EM_prev[0] - EM[0])
                if random.random() < np.exp(delta):
                    var = 1
                    comm.send(var, dest=my_rank - 1, tag=222)
                    comm.Send(EM, dest=my_rank - 1, tag=23)
                    comm.Send(config, dest=my_rank - 1, tag=24)
                    EM = EM_prev
                    config = config_prev
                else:
                    var = 0
                    comm.send(var, dest=my_rank - 1, tag=222)
        if iw % 2 == 0 and i % gap == 0:
            S[int(i/gap)] = config
        if iw % 2 == 1 and i % gap == 0:
            q_p = np.sum(np.sum(S[int(i / gap)] * config, axis=1)) / Ns
            q_p2 = q_p ** 2 + q_p2
            q_p4 = q_p ** 4 + q_p4
        Ene = EM[0]
        Mag = EM[1]
        Ma2 = EM[2]
        Me = EM[3]
        Mt1b = EM[4]
        Mt2 = EM[5]
        Metheta = EM[6]

        E1 = E1 + Ene
        E2 = E2 + Ene * Ene
        M1 = M1 + Mag
        M2 = M2 + Mag * Mag
        ma2 = ma2 + Ma2
        me = me + Me
        mt1b = mt1b + Mt1b
        mt2 = mt2 + Mt2
        metheta = metheta + Metheta
        me_2 = me_2 + Me * Me
        mt2_2 = mt2_2 + Mt2 * Mt2
    if iw % 2 == 1:
        w_count = w_count + 1
        q_p2 = q_p2 / (mcSteps / gap)
        q_p4 = q_p4 / (mcSteps / gap)

        Q2 = Q2 + q_p2
        Q4 = Q4 + q_p4
        Q2_sq = Q2_sq + q_p2 ** 2
    E = E1 / (mcSteps * Ns)
    M = M1 / (mcSteps * Ns)
    C = (E2 / (mcSteps * Ns) - (E1 ** 2) / (mcSteps * mcSteps * Ns)) * iT * iT
    X = (M2 / (mcSteps * Ns) - (M1 ** 2) / (mcSteps * mcSteps * Ns)) * iT
    MA2 = ma2 / (mcSteps * Ns)
    ME = me / (mcSteps * Ns)
    MT1B = mt1b / (mcSteps * Ns)
    MT2 = mt2 / (mcSteps * Ns)
    METHETA = metheta / (mcSteps * Ns)
    X_ME = (me_2 / (mcSteps * Ns) - (me ** 2) / (mcSteps * mcSteps * Ns)) * iT
    X_MT2 = (mt2_2 / (mcSteps * Ns) - (mt2 ** 2) / (mcSteps * mcSteps * Ns)) * iT

    Ew = Ew + E
    Mw = Mw + M
    Cw = Cw + C
    Xw = Xw + X
    MA2w = MA2w + MA2
    MEw = MEw + ME
    MT1Bw = MT1Bw + MT1B
    MT2w = MT2w + MT2
    METHETAw = METHETAw + METHETA
    X_MEw = X_MEw + X_ME
    X_MT2w = X_MT2w + X_MT2
Q2 = Q2 / (w / 2)
Q4 = Q4 / (w / 2)
Q2_sq = Q2_sq / (w / 2)

X_q = Ns * Q2
g4 = 1.5 - 0.5 * Q4 / Q2 ** 2
g2 = (Q2_sq - Q2 ** 2) / Q2 ** 2
G = g2 / (2.0 - 2.0 * g4)

Ew = Ew / w
Mw = Mw / w
Cw = Cw / w
Xw = Xw / w
MA2w = MA2w / w
MEw = MEw / w
MT1Bw = MT1Bw / w
MT2w = MT2w / w
METHETAw = METHETAw / w
X_MEw = X_MEw / w
X_MT2w = X_MT2w / w

if my_rank != 0:
    comm.send(Ew, dest=0, tag=5)
    comm.send(Mw, dest=0, tag=6)
    comm.send(Cw, dest=0, tag=7)
    comm.send(Xw, dest=0, tag=8)
    comm.send(MA2w, dest=0, tag=9)
    comm.send(MEw, dest=0, tag=10)
    comm.send(MT1Bw, dest=0, tag=11)
    comm.send(MT2w, dest=0, tag=12)
    comm.send(METHETAw, dest=0, tag=13)
    comm.send(X_MEw, dest=0, tag=14)
    comm.send(X_MT2w, dest=0, tag=15)
    comm.send(X_q, dest=0, tag=16)
    comm.send(g4, dest=0, tag=17)
    comm.send(g2, dest=0, tag=18)
    comm.send(G, dest=0, tag=19)
else:
    Et[0] = Ew
    Mt[0] = Mw
    Ct[0] = Cw
    Xt[0] = Xw
    MA2t[0] = MA2w
    MEt[0] = MEw
    MT1Bt[0] = MT1Bw
    MT2t[0] = MT2w
    METHETAt[0] = METHETAw
    X_MEt[0] = X_MEw
    X_MT2t[0] = X_MT2w
    X_qt[0] = X_q
    g4t[0] = g4
    g2t[0] = g2
    Gt[0] = G
    for procid in range(1, p):
        Et[procid] = comm.recv(source=procid, tag=5)
        Mt[procid] = comm.recv(source=procid, tag=6)
        Ct[procid] = comm.recv(source=procid, tag=7)
        Xt[procid] = comm.recv(source=procid, tag=8)
        MA2t[procid] = comm.recv(source=procid, tag=9)
        MEt[procid] = comm.recv(source=procid, tag=10)
        MT1Bt[procid] = comm.recv(source=procid, tag=11)
        MT2t[procid] = comm.recv(source=procid, tag=12)
        METHETAt[procid] = comm.recv(source=procid, tag=13)
        X_MEt[procid] = comm.recv(source=procid, tag=14)
        X_MT2t[procid] = comm.recv(source=procid, tag=15)
        X_qt[procid] = comm.recv(source=procid, tag=16)
        g4t[procid] = comm.recv(source=procid, tag=17)
        g2t[procid] = comm.recv(source=procid, tag=18)
        Gt[procid] = comm.recv(source=procid, tag=19)
    Order_p = np.array(
        [
            Et,
            Mt,
            Ct,
            Xt,
            MA2t,
            MEt,
            MT1Bt,
            MT2t,
            METHETAt,
            X_MEt,
            X_MT2t,
            X_qt,
            g4t,
            g2t,
            Gt,
        ]
    )
    pickle.dump(Order_p, open("Order_%s_%s.p" % (xd, num), "wb"))
print(time.process_time() - start, my_rank, w_count)
