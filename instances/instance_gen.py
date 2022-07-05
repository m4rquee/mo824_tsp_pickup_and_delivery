import numpy as np
from math import ceil

np.random.seed(42)
square_dim = 1000
tests = 3
N = [5, 6, 7, 8, 9, 10, 15, 20, 25, 50, 75, 100, 150, 200, 300]


def dist(p, q):
    return ceil(np.linalg.norm(p - q))


def print_to_file(n, label, points):
    with open(f'l2{n}{label}.dig', 'w') as f_out:
        print('nnodes narcs type', file=f_out)
        v = len(points)
        print(f'{v} {v * (v - 1)} digraph', file=f_out)

        print('nodename posx posy', file=f_out)
        for i, point in enumerate(points):
            print(f'{i} {point[0]} {point[1]}', file=f_out)

        print('tail head weight', file=f_out)
        for i, tail in enumerate(points):
            for j, head in enumerate(points):
                if i == j: continue
                print(f'{i} {j} {dist(tail, head)}', file=f_out)


a_num = ord('a')
for n in N:
    for label in range(a_num, a_num + tests):
        points = np.random.randint(square_dim + 1, size=(2 * n + 2, 2))
        print_to_file(n, chr(label), points)
