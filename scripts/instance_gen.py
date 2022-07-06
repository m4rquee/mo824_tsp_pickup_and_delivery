import numpy as np
from math import ceil

np.random.seed(42)
square_dim = 1000
tests = 3
N = [5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 75, 100, 150, 200, 300]
P = [1, 2, np.inf]


def dist(p, q, ord): return ceil(np.linalg.norm(p - q, ord=ord))


def print_to_file(n, points, ord):
    with open(f'../instances/l{ord}-{n:0>3d}.dig', 'w') as f_out:
        print('nnodes narcs type', file=f_out)
        v = len(points)
        print(f'{v} {v * (v - 1)} digraph', file=f_out)

        print('nodename posx posy', file=f_out)
        print(f'{0} {points[0][0]} {points[0][1]}', file=f_out)  # source
        print(f'{v - 1} {points[v - 1][0]} {points[v - 1][1]}', file=f_out)  # target
        for i in range(1, n + 1):
            print(f'{i} {points[i][0]} {points[i][1]}', file=f_out)
            print(f'{i + n} {points[i + n][0]} {points[i + n][1]}', file=f_out)

        print('tail head weight', file=f_out)
        for i, tail in enumerate(points):
            for j, head in enumerate(points):
                if i == j: continue
                print(f'{i} {j} {dist(tail, head, ord)}', file=f_out)


for n in N:
    for p in P:
        points = np.random.randint(square_dim + 1, size=(2 * n + 2, 2))
        print_to_file(n, points, p)
