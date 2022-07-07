import glob
import subprocess
import sys

solver = sys.argv[1]
T = eval(sys.argv[2])
test_type = sys.argv[3]
output = []

for instance in glob.glob(f'../instances/{test_type}/*.dig'):
    ret = subprocess.run([f'../bin/{solver}', instance, str(T)],
                         stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    LB, UB = map(float, ret.stderr.split())
    GAP = (UB - LB) / UB
    print(f'instance={instance}; LB={LB}; UB={UB}; GAP={GAP}')
    output.append(f'{instance},{LB},{UB},{GAP}')

with open(f'../output/{solver}-{test_type}.csv', 'w') as f_out:
    output.sort()
    print('instance,LB,UB,GAP', file=f_out)
    print(*output, sep='\n', file=f_out)
