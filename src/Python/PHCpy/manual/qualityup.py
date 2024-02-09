"""
Can multithreading compensate for the overhead of double double arithmetic?
If we can afford the time for a sequential run, by how much can we increase
the precision in a multithreaded run in the same time or less?
Wall clock times are compared.  In both runs, the seed is fixed to work 
with the same random constants in the homotopies.
"""
from datetime import datetime
from phcpy.dimension import get_core_count, set_seed, get_seed
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
s = solve(c7, precision='d')
timestop = datetime.now()
elapsed_onecore = timestop - timestart
print('computed', len(s), 'solutions in time', elapsed_onecore)
print('solving on', nbcores, 'cores ...')
set_seed(2024)
timestart = datetime.now()
s = solve(c7, tasks=nbcores, precision='dd')
timestop = datetime.now()
elapsed_manycores = timestop - timestart
print('computed', len(s), 'solutions in time', elapsed_manycores)
if elapsed_manycores < elapsed_onecore:
    print('-> doubled precision with', nbcores, 'cores in less time!')
else:
    print('-> need more than', nbcores, 'to compensate for overhead.')
