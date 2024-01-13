"""
Excerpt from the user manual of phcpy.
"""
import matplotlib.pyplot as plt
from phcpy.dimension import set_seed, get_seed
from phcpy.solutions import strsol2dict
from phcpy.starters import total_degree_start_system
from phcpy.trackers import initialize_double_tracker
from phcpy.trackers import initialize_double_solution
from phcpy.trackers import next_double_solution

set_seed(12871)
print('the seed :', get_seed())
p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
print('constructing a total degree start system ...')
q, qsols = total_degree_start_system(p)
print('number of start solutions :', len(qsols))
initialize_double_tracker(p, q, False)
plt.ion()
fig = plt.figure()
for k in range(len(qsols)):
    if(k == 0):
        axs = fig.add_subplot(221)
    elif(k == 1):
        axs = fig.add_subplot(222)
    elif(k == 2):
        axs = fig.add_subplot(223)
    elif(k == 3):
        axs = fig.add_subplot(224)
    startsol = qsols[k]
    initialize_double_solution(len(p),startsol)
    dictsol = strsol2dict(startsol)
    xpoints =  [dictsol['x']]
    ypoints =  [dictsol['y']]
    for k in range(300):
        ns = next_double_solution()
        dictsol = strsol2dict(ns)
        xpoints.append(dictsol['x'])
        ypoints.append(dictsol['y'])
        tval = dictsol['t'].real
        if(tval == 1.0):
            break
    print(ns)
    xre = [point.real for point in xpoints]
    yre = [point.real for point in ypoints]
    axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
    axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
    dots, = axs.plot(xre,yre,'r-')
    dots, = axs.plot(xre,yre,'ro')
    fig.canvas.draw()
fig.canvas.draw()
ans = input('hit return to exit')
