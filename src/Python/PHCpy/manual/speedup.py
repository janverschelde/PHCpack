"""
The script demonstrates the speedup of multithreading.
To measure speedup, we compare wall clock times.
In both runs, the seed is fixed to work 
with the same random constants in the homotopies.
"""
from datetime import datetime
from phcpy.dimension import get_core_count, set_seed
from phcpy.polynomials import set_double_system
from phcpy.families import cyclic
from phcpy.volumes import mixed_volume
from phcpy.solver import solve

nbcores = get_core_count()
print('number of available cores :', nbcores)
print('the cyclic 7-roots problem :')
c7 = cyclic(7)
for pol in c7:
    print(pol)
print('the mixed volume :', mixed_volume(c7))
print('solving on one core ...')
set_seed(2024)
timestart = datetime.now()
s = solve(c7)
timestop = datetime.now()
elapsed_onecore = timestop - timestart
print('computed', len(s), 'solutions in time', elapsed_onecore)
print('solving on', nbcores, 'cores ...')
set_seed(2024)
timestart = datetime.now()
s = solve(c7, tasks=nbcores)
timestop = datetime.now()
elapsed_manycores = timestop - timestart
print('computed', len(s), 'solutions in time', elapsed_manycores)
speedup = elapsed_onecore/elapsed_manycores
print('the speedup : {:2.2f}'.format(speedup))
