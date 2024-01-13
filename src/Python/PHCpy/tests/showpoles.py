"""
Produces the same plots as showpaths.py,
but now computed with the apriori step size control
and with the plots of the location of the closest poles.
"""
import matplotlib.pyplot as plt
from phcpy.dimension import set_seed, get_seed
from phcpy.solutions import strsol2dict, clear_double_solutions
from phcpy.starters import total_degree_start_system
from phcpy.curves import set_default_parameters, write_parameters
from phcpy.curves import initialize_double_artificial_homotopy
from phcpy.curves import set_double_solution, get_double_solution
from phcpy.curves import double_predict_correct
from phcpy.curves import double_t_value, double_closest_pole

set_seed(12871)
print('the seed :', get_seed())
p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
print('constructing a total degree start system ...')
q, qsols = total_degree_start_system(p)
clear_double_solutions()
print('number of start solutions :', len(qsols))
set_default_parameters()
write_parameters()
initialize_double_artificial_homotopy(p, q, False)
plt.ion()
fig1 = plt.figure()
allpoles = []
for k in range(len(qsols)):
    if(k == 0):
        axs = fig1.add_subplot(221)
    elif(k == 1):
        axs = fig1.add_subplot(222)
    elif(k == 2):
        axs = fig1.add_subplot(223)
    elif(k == 3):
        axs = fig1.add_subplot(224)
    startsol = qsols[k]
    set_double_solution(len(p), startsol)
    dictsol = strsol2dict(startsol)
    xpoints =  [dictsol['x']]
    ypoints =  [dictsol['y']]
    poles = []
    for k in range(100):
        ns = get_double_solution()
        dictsol = strsol2dict(ns)
        xpoints.append(dictsol['x'])
        ypoints.append(dictsol['y'])
        tval = dictsol['t'].real
        if(tval == 1.0):
            break
        double_predict_correct()
        pole = double_closest_pole()
        tval = double_t_value()
        locp = (tval+pole[0], pole[1])
        poles.append(locp)
    print(ns)
    xre = [point.real for point in xpoints]
    yre = [point.real for point in ypoints]
    axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
    axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
    dots, = axs.plot(xre,yre,'b-')
    dots, = axs.plot(xre,yre,'bo')
    fig1.canvas.draw()
    allpoles.append(poles)
fig1.canvas.draw()
ans = input('hit return for the next figure')
fig2 = plt.figure()
for k in range(len(qsols)):
    if(k == 0):
        axs = fig2.add_subplot(221)
    elif(k == 1):
        axs = fig2.add_subplot(222)
    elif(k == 2):
        axs = fig2.add_subplot(223)
    elif(k == 3):
        axs = fig2.add_subplot(224)
    poles = allpoles[k]
    pl0 = [pole[0] for pole in poles]
    pl1 = [pole[1] for pole in poles]
    axs.set_xlim(-0.2, 1.2)
    axs.set_ylim(-0.5, 0.5)
    dots, = axs.plot(pl0,pl1,'r+')
    fig2.canvas.draw()
fig2.canvas.draw()
ans = input('hit return to exit')
