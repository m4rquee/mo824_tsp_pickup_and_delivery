import glob
import subprocess

import numpy as np

EPS = [0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 5, 10]
P = [1, 2, np.inf]
for eps in EPS:
    output = []
    for p in P:
        for instance in glob.glob(f'../instances/l{p}*.dig'):
            ret = subprocess.run(['../bin/pickup_delivery_heuristic', instance, str(eps)],
                                 stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            LB, UB_H, UB_L = map(float, ret.stderr.split())
            GAP = (UB_L - LB) / UB_L
            print(f'eps={eps}; instance={instance}; LB={LB}; '
                  f'heuristic UB={UB_H}; UB with local-search={UB_L}; GAP={GAP}')
            output.append(f'{instance},{LB},{UB_H},{UB_L},{GAP}')

    with open(f'../output/eps-{eps}.csv', 'w') as f_out:
        output.sort()
        print('instance,LB,UB_H,UB_L,GAP', file=f_out)
        print(*output, sep='\n', file=f_out)
