
require(["nbextensions/snippets_menu/main"], function (snippets_menu) {
    snippets_menu.options['menus'] = [{
        'name' : 'PHCpy',
        'sub-menu-direction' : 'left',
        'sub-menu' : [
    {
    'name' : 'blackbox solver',
    'sub-menu' : [
        {
        'name' : "solving trinomials",
        'sub-menu' : [
            {
            'name' : "solving a random case",
            'snippet' : ["from phcpy.solver import random_trinomials", "f = random_trinomials()", "for pol in f: print pol", "from phcpy.solver import solve", "sols = solve(f, verbose=False)", "for sol in sols: print sol", "print len(sols), \"solutions found\""],
            },
            {
            'name' : "solving a specific case",
            'snippet' : ["f = ['x^2*y^2 + 2*x - 1;', 'x^2*y^2 - 3*y + 1;']", "from phcpy.solver import solve", "sols = solve(f)", "for sol in sols: print sol"],
            }],
        },
        {
        'name' : "representations of isolated solutions",
        'sub-menu' : [
            {
            'name' : "from string to dictionary",
            'snippet' : ["p = ['x + y - 1;', '2*x - 3*y + 1;']", "from phcpy.solver import solve", "sols = solve(p)", "print sols[0]", "from phcpy.solutions import strsol2dict", "dsol = strsol2dict(sols[0])", "print dsol.keys()", "for key in dsol.keys(): print 'the value for', key, 'is', dsol[key]"],
            },
            {
            'name' : "verify by evaluation",
            'snippet' : ["p = ['x + y - 1;', '2*x - 3*y + 1;']", "from phcpy.solver import solve", "sols = solve(p)", "from phcpy.solutions import strsol2dict, evaluate", "dsol = strsol2dict(sols[0])", "eva = evaluate(p, dsol)", "for val in eva: print val"],
            },
            {
            'name' : "making a solution",
            'snippet' : ["from phcpy.solutions import make_solution", "s0 = make_solution(['x', 'y'], [float(3.14), complex(0, 2.7)])", "print s0", "s1 = make_solution(['x', 'y'], [int(2), int(3)])", "print s1"],
            },
            {
            'name' : "filtering solution lists",
            'snippet' : ["from phcpy.solutions import make_solution, is_real, filter_real",  "s0 = make_solution(['x', 'y'], [float(3.14), complex(0, 2.7)])", "print is_real(s0, 1.0e-8)", "s1 = make_solution(['x', 'y'], [int(2), int(3)])", "print is_real(s1, 1.0e-8)", "realsols = filter_real([s0, s1], 1.0e-8, 'select')", "for sol in realsols: print sol"],
            },
            {
            'name' : "coordinates, names and values",
            'snippet' : ["from phcpy.solver import solve", "p = ['x^2*y^2 + x + 1;', 'x^2*y^2 + y + 1;']", "s = solve(p)", "print(s[0])", "from phcpy.solutions import coordinates, make_solution", "(names, values) = coordinates(s[0])", "print names", "print values", "s0 = make_solution(names, values)", "print s0"],
            }],
        },
        {
        'name' : "reproducible runs with fixed seeds",
        'sub-menu' : [
            {
            'name' : "fixing and retrieving the seed",
            'snippet' : ["from phcpy.phcpy2c2 import py2c_set_seed", "py2c_set_seed(2013)", "from phcpy.phcpy2c2 import py2c_get_seed", "print py2c_get_seed()"],
            }],
        },
        {
        'name' : "shared memory parallelism",
        'sub-menu' : [
            {
            'name' : "solving with 4 tasks",
            'snippet' : ["from phcpy.solver import solve", "from phcpy.families import cyclic", "nbrt = 4 # number of tasks", "pols = cyclic(6)", "print 'solving the cyclic 6-roots problem :'", "for pol in pols: print pol", "from time import time", "starttime = time()", "sols = solve(pols, verbose=False)", "stoptime = time()", "elapsed = stoptime - starttime", "print 'solving with no multitasking :'", "print 'elapsed time : %.2f seconds' % elapsed", "starttime = time()", "sols = solve(pols, verbose=False, tasks=nbrt)", "stoptime = time()", "elapsed = stoptime - starttime", "print 'solving with %d tasks :' % nbrt", "print 'elapsed time : %.2f seconds' % elapsed" ],
            }],
        },
        {
        'name' : "root counting methods",
        'sub-menu' : [
            {
            'name' : "four different root counts",
            'snippet' : ["f = ['x^3*y^2 + x*y^2 + x^2;', 'x^5 + x^2*y^3 + y^2;']", "from phcpy.solver import total_degree", "print 'the total degree :', total_degree(f)", "from phcpy.solver import m_homogeneous_bezout_number as mbz", "(bz, part) = mbz(f)", "print 'a multihomogeneous Bezout number :', bz", "from phcpy.solver import linear_product_root_count as lrc", "print 'a linear-product root count :', lrc(f, verbose=False)", "from phcpy.solver import mixed_volume", "(mv, smv) = mixed_volume(f, stable=True)", "print 'the mixed volume :', mv", "print 'the stable mixed volume :', smv"],
            }],
        },
        {
        'name' : "Newton's method and deflation",
        'sub-menu' : [
            {
            'name' : "the Griewank-Osborne example",
            'snippet' : ["p = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']", "from phcpy.solutions import make_solution", "s = make_solution(['x', 'y'],[float(1.0e-6), float(1.0e-6)])", "print s", "from phcpy.solver import newton_step", "s2 = newton_step(p, [s])", "print s2[0]", "s3 = newton_step(p, s2)", "print s3[0]", "from phcpy.solver import standard_deflate", "sd = standard_deflate(p, [s])", "print sd[0]"],
            },
            {
            'name' : "deflating an overconstrained system",
            'snippet' : ["from phcpy.solutions import make_solution", "from phcpy.solver import standard_deflate", "sol = make_solution(['x', 'y'], [float(1.0e-6), float(1.0e-6)])", "print sol", "pols = ['x**2;', 'x*y;', 'y**2;']", "sols = standard_deflate(pols, [sol], tolrnk=1.0e-8)", "print sols[0]", "sols = standard_deflate(pols, [sol], tolrnk=1.0e-4)", "print sols[0]"],
            }],
        },
        {
        'name' : "equation and variable scaling",
        'sub-menu' : [
            {
            'name' : "solving without scaling",
            'snippet' : ["from phcpy.solver import solve", "p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']", "psols = solve(p, verbose=False)", "print psols[0]"],
            },
            {
            'name' : "solving after scaling",
            'snippet' : ["p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']", "from phcpy.solver import standard_scale_system as scalesys", "from phcpy.solver import standard_scale_solutions as scalesols", "(q, c) = scalesys(p)", "for pol in q: print pol", "qsols = solve(q, verbose=False)", "ssols = scalesols(len(q), qsols, c)", "for sol in ssols: print sol"],
            }],
        }],
    },

    {
    'name' : 'path trackers',
    'sub-menu' : [
        {
        'name' : "a simple example",
        'sub-menu' : [
            {
            'name' : "a total degree start system",
            'snippet' : ["from phcpy.solver import total_degree", "from phcpy.solver import total_degree_start_system", "from phcpy.trackers import track", "p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']", "d = total_degree(p)", "print 'the total degree :', d", "(q, qsols) = total_degree_start_system(p)", "print 'the number of start solutions :', len(qsols)", "print 'the start system :', q", "s = track(p, q, qsols)", "print 'the number of solutions :', len(s)", "for sol in s: print sol"],
            },
            {
            'name' : "track one solution path",
            'snippet' : ["from phcpy.solver import total_degree_start_system", "from phcpy.trackers import track", "p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']", "(q,qsols) = total_degree_start_system(p)", "s1 = track(p, q, [qsols[2]])", "print s1[0]", "s2 = track(p, q,[qsols[2]])", "print s2[0]"],
            }],
        },
        {
        'name' : "fixing the gamma constant",
        'sub-menu' : [
            {
            'name' : "specifying the gamma parameter",
            'snippet' : ["from phcpy.solver import total_degree_start_system", "from phcpy.trackers import track", "p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']", "(q, qsols) = total_degree_start_system(p)", "s3 = track(p, q, [qsols[2]], gamma=complex(0.824372806319,0.56604723848934))", "print 'the solution at the end:'", "print s3[0]"],
            }],
        },
        {
        'name' : "give the next solution on a path",
        'sub-menu' : [
            {
            'name' : "tracking with a generator",
            'snippet' : ["from phcpy.solver import total_degree_start_system", "p = ['x**2 + 4*x**2 - 4;', '2*y**2 - x;']", "(q, s) = total_degree_start_system(p)", "from phcpy.trackers import initialize_standard_tracker", "from phcpy.trackers import initialize_standard_solution", "from phcpy.trackers import next_standard_solution", "initialize_standard_tracker(p, q)", "initialize_standard_solution(len(p), s[0])", "s1 = next_standard_solution()", "print 'the next point on the solution path :'", "print s1", "print next_standard_solution()", "print next_standard_solution()", "initialize_standard_solution(len(p), s[1])", "points = [next_standard_solution() for i in range(11)]", "from phcpy.solutions import strsol2dict", "dicpts = [strsol2dict(sol) for sol in points]", "xvals = [sol['x'] for sol in dicpts]", "for x in xvals: print(x)"],
            },
            {
            'name' : "plotting trajectories",
            'snippet' : ["p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']", "print 'constructing a total degree start system ...'", "from phcpy.solver import total_degree_start_system as tds", "q, qsols = tds(p)", "print 'number of start solutions :', len(qsols)", "from phcpy.trackers import initialize_standard_tracker", "from phcpy.trackers import initialize_standard_solution", "from phcpy.trackers import next_standard_solution", "initialize_standard_tracker(p, q, False)", "from phcpy.solutions import strsol2dict", "import matplotlib.pyplot as plt", "plt.ion()", "fig = plt.figure()", "for k in range(len(qsols)):", "    if(k == 0):", "        axs = fig.add_subplot(221)", "    elif(k == 1):", "        axs = fig.add_subplot(222)", "    elif(k == 2):", "        axs = fig.add_subplot(223)", "    elif(k == 3):", "        axs = fig.add_subplot(224)", "    startsol = qsols[k]", "    initialize_standard_solution(len(p),startsol)", "    dictsol = strsol2dict(startsol)", "    xpoints =  [dictsol['x']]", "    ypoints =  [dictsol['y']]", "    for k in range(300):", "        ns = next_standard_solution()", "        dictsol = strsol2dict(ns)", "        xpoints.append(dictsol['x'])", "        ypoints.append(dictsol['y'])", "        tval = eval(dictsol['t'].lstrip().split(' ')[0])", "        if(tval == 1.0):", "            break", "    print(ns)", "    xre = [point.real for point in xpoints]", "    yre = [point.real for point in ypoints]", "    axs.set_xlim(min(xre)-0.3, max(xre)+0.3)", "    axs.set_ylim(min(yre)-0.3, max(yre)+0.3)", "    dots, = axs.plot(xre,yre,'r-')", "    fig.canvas.draw()", "fig.canvas.draw()"],
            }],
        },
        {
        'name' : "solving with polyhedral homotopies",
        'sub-menu' : [
            {
            'name' : "solving a random coefficient system",
            'snippet' : ["p = ['x^3*y^2 - 3*x^3 + 7;','x*y^3 + 6*y^3 - 9;']", "from phcpy.solver import mixed_volume", "print 'the mixed volume :', mixed_volume(p)", "from phcpy.solver import random_coefficient_system", "(q,qsols) = random_coefficient_system(verbose=False)", "print 'the number of start solutions :', len(qsols)", "from phcpy.trackers import track", "psols = track(p, q, qsols)", "print 'the number of solutions at the end :', len(psols)", "for sol in psols: print sol"],
            }],
        },
        {
        'name' : "Newton's method at higher precision",
        'sub-menu' : [
            {
            'name' : "using a linear-product start system",
            'snippet' : ["p = ['x*y^3 + y - 2;', 'x^3*y + x - 8;']", "from phcpy.solver import linear_product_root_count", "r = linear_product_root_count(p)", "from phcpy.solver import random_linear_product_system", "(q, qsols) = random_linear_product_system(p)", "print 'the number of start solutions :', len(qsols)", "from phcpy.trackers import track", "psols = track(p,q,qsols)", "print 'the number of end solutions :', len(psols)", "from phcpy.solver import newton_step", "psols_dd = newton_step(p,psols,precision='dd')", "print 'the solutions in double double precision :'", "for sol in psols_dd: print sol"],
            }],
        },
        {
        'name' : "multitasked path tracking",
        'sub-menu' : [
            {
            'name' : "tracking with 4 tasks",
            'snippet' : ["from phcpy.solver import random_linear_product_system as rlps", "from phcpy.families import noon", "from phcpy.trackers import track", "nbrt = 4 # number of tasks", "pols = noon(5)", "print 'solving the 5-variable Noonburg system :'", "for pol in pols: print pol", "(startpols, startsols) = rlps(pols)", "print 'number of paths :', len(startsols)", "from time import time", "starttime = time()", "sols = track(pols, startpols, startsols)", "stoptime = time()", "elapsed = stoptime - starttime", "print 'elapsed time with no multitasking : %.2f seconds' % elapsed", "starttime = time()", "sols = track(pols, startpols, startsols, tasks=nbrt)", "stoptime = time()", "elapsed = stoptime - starttime", "print 'elapsed time with %d tasks : %.2f seconds' % (nbrt, elapsed)"],
            }],
        },
        {
        'name' : "sweep homotopies",
        'sub-menu' : [
            {
            'name' : "towards a quadratic turning point",
            'snippet' : ["circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']",
            "from phcpy.solutions import make_solution as makesol", "first = makesol(['x', 'y', 's'], [float(1), float(0), float(0)])", "second = makesol(['x', 'y', 's'], [float(-1), float(0), float(0)])", "startsols = [first, second]", "from phcpy.sweepers import standard_real_sweep as sweep", "newsols = sweep(circle, startsols)", "print newsols[0]"],
            }],
        },
        {
        'name' : "real versus complex sweeps",
        'sub-menu' : [
            {
            'name' : "complex parameter homotopy continuation",
            'snippet' : ["circle = ['x^2 + y^2 - 1;']", "from phcpy.solutions import make_solution as makesol", "first = makesol(['x', 'y'], [float(1), float(0)])", "second = makesol(['x', 'y'], [float(-1), float(0)])", "startsols = [first, second]", "par = ['y']", "start = [0, 0] ", "target = [2, 0]", "from phcpy.sweepers import standard_complex_sweep as sweep", "newsols = sweep(circle, startsols, 2, par, start, target)", "print newsols[0]"],
            }],
        },
        {
        'name' : "a polyhedral end game",
        'sub-menu' : [
            {
            'name' : "numerical tropism computation",
            'snippet' : ["from phcpy.tuning import order_endgame_extrapolator_set as set", "set(4)", "from phcpy.tuning import order_endgame_extrapolator_get as get", "get()", "f = ['x + y^3 - 1;', 'x + y^3 + 1;']", "from phcpy.solver import mixed_volume as mv", "from phcpy.solver import random_coefficient_system as rcs", "print 'the mixed volume :', mv(f)", "(g, gsols) = rcs(f)", "print 'the number of start solutions :', len(gsols)", "from phcpy.trackers import standard_double_track as track", "sols = track(f, g, gsols)", "from phcpy.tropisms import standard_retrieve as retrieve", "(w, d, e) = retrieve(len(sols), len(f))", "print 'the numerical direction :', d", "print 'the error :', e", "print w"],
            }],
        }],
    },

    {
    'name' : 'solution sets',
    'sub-menu' : [
        {
        'name' : "witness sets",
        'sub-menu' : [
            {
            'name' : "embedding the twisted cubic",
            'snippet' : ["twisted = ['x^2 - y;', 'x^3 - z;']", "from phcpy.sets import embed", "e = embed(3,1,twisted)", "for pol in e: print pol"],
            },
            {
            'name' : "a witness set for the twisted cubic",
            'snippet' : ["twisted = ['x^2 - y;', 'x^3 - z;']", "from phcpy.sets import embed", "e = embed(3,1,twisted)", "from phcpy.solver import solve", "s = solve(e, verbose=False)", "print 'number of generic points :', len(s)", "for sol in s: print sol"],
            }],
        },
        {
        'name' : "homotopy membership test",
        'sub-menu' : [
            {
            'name' : "cyclic 4-roots on coordinates",
            'snippet' : ["from phcpy.families import cyclic", "c4 = cyclic(4)", "from phcpy.sets import embed", "c4e1 = embed(4, 1, c4)", "from phcpy.solver import solve", "sols = solve(c4e1)", "from phcpy.solutions import filter_zero_coordinates as filter", "genpts = filter(sols, 'zz1', 1.0e-8, 'select')", "for sol in genpts: print sol", "point0 = [1, 0, -1, 0, 1, 0, -1, 0]", "from phcpy.sets import membertest", "print 'point0 :', point0", "print 'Is point0 a member ?', membertest(c4e1, genpts, 1, point0)", "point1 = [1, 0, 1, 0, -1, 0, -1, 0]", "print 'point1 :', point1", "print 'Is point1 a member ?', membertest(c4e1, genpts, 1, point1)"],
            },
            {
            'name' : "cyclic 4-roots on solutions",
            'snippet' : ["from phcpy.families import cyclic", "c4 = cyclic(4)", "from phcpy.sets import embed", "c4e1 = embed(4, 1, c4)", "from phcpy.solver import solve", "sols = solve(c4e1)", "from phcpy.solutions import filter_zero_coordinates as filter", "genpts = filter(sols, 'zz1', 1.0e-8, 'select')", "for sol in genpts: print sol", "names = ['x0', 'x1', 'x2', 'x3']", "coord0 = [complex(1, 0), complex(-1, 0), complex(1, 0), complex(-1, 0)]", "from phcpy.solutions import make_solution", "point0 = make_solution(names, coord0)", "from phcpy.sets import is_member", "print 'point0 :'", "print point0", "print 'Is point0 a member ?', is_member(c4e1, genpts, 1, point0, verbose=False)", "coord1 = [complex(1, 0), complex(1, 0), complex(-1, 0), complex(-1, 0)]", "point1 = make_solution(names, coord1)", "print 'point1 :'", "print point1", "print 'Is point1 a member ?', is_member(c4e1, genpts, 1, point1, verbose=False)"],
            }],
        },
        {
        'name' : "cascade of homotopies",
        'sub-menu' : [
            {
            'name' : "an illustrative example",
            'snippet' : ["pol1 = '(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);'", "pol2 = '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);'", "pol3 = '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);'", "pols = [pol1, pol2, pol3]", "from phcpy.cascades import run_cascade", "otp = run_cascade(3, 2, pols)", "dims = otp.keys()", "dims.sort(reverse=True)", "for dim in dims: print 'number of solutions at dimension', dim, ':', len(otp[dim][1])"]
            },
            {
            'name' : "a Laurent system",
            'snippet' : ["pol1 = '0.710358341606049*t1 + 0.46*t2 - 0.41*t3 + 0.240761300555115 + 1.07248215701824*I;'", "pol2 = 't2*(-0.11 + 0.49*I) + 0.41*t3 - 0.502195181179589*t4 + 0.41*t5;'", "pol3 = '0.502195181179589*t4 + t5*(-0.0980434782608696 + 0.436739130434783*I) - 0.775518556663656*t6 - 1.2;'", "pol4 = '0.710358341606049*t1**(-1) + 0.46*t2**(-1) - 0.41*t3**(-1) + 0.240761300555115 - 1.07248215701824*I;'", "pol5 = 't2**(-1)*(-0.11 - 0.49*I) + 0.41*t3**(-1) - 0.502195181179589*t4**(-1) + 0.41*t5**(-1);'", "pol6 = '0.502195181179589*t4**(-1) + t5**(-1)*(-0.0980434782608696 - 0.436739130434783*I) - 0.775518556663656*t6**(-1) - 1.2;'", "pols = [pol1, pol2, pol3, pol4, pol5, pol6]", "from phcpy.cascades import run_cascade", "otp = run_cascade(6, 1, pols, islaurent=True)", "(epols, esols) = otp[1]", "print 'the generic points at a 1-dimensional curve :'", "for sol in esols: print sol"]
            }],
        },
        {
        'name' : "factoring into irreducibles",
        'sub-menu' : [
            {
            'name' : "factoring a cubic polynomial",
            'snippet' : ["p = '(x+1)*(x^2 + y^2 + 1);'", "from phcpy.sets import witness_set_of_hypersurface as wh", "(w, s) = wh(2, p)", "print 'number of witness points :', len(s)", "from phcpy.factor import factor", "f = factor(1, w, s)", "for fact in f: print fact"],
            }],
        },
        {
        'name' : "numerical irreducible decomposition",
        'sub-menu' : [
            {
            'name' : "an example",
            'snippet' : [ "pol0 = '(x1-1)*(x1-2)*(x1-3)*(x1-4);'", "pol1 = '(x1-1)*(x2-1)*(x2-2)*(x2-3);'", "pol2 = '(x1-1)*(x1-2)*(x3-1)*(x3-2);'", "pol3 = '(x1-1)*(x2-1)*(x3-1)*(x4-1);'", "pols = [pol0, pol1, pol2, pol3]", "from phcpy.factor import solve, write_decomposition", "deco = solve(4, 3, pols, verbose=False)", "write_decomposition(deco)"],
            }],
        },
        {
        'name' : "diagonal homotopies",
        'sub-menu' : [
            {
            'name' : "sphere intersected with a cylinder",
            'snippet' : ["sph = 'x^2 + y^2 + z^2 - 1;'", "cyl = 'x^2 + y - y + (z - 0.5)^2 - 1;'", "from phcpy.sets import witness_set_of_hypersurface as witsurf", "sphwit = witsurf(3, sph)", "spheqs, sphpts = sphwit", "cylwit = witsurf(3, cyl)", "cyleqs, cylpts = cylwit", "from phcpy.diagonal import diagonal_solver as diagsolve", "quawit = diagsolve(3, 2, spheqs, sphpts, 2, cyleqs, cylpts, verbose=False)", "quaeqs, quapts = quawit", "for pol in quaeqs: print pol", "for sol in quapts: print sol"],
            }],
        }],
    },
    {
    'name' : 'families of systems',
    'sub-menu' : [
        {
        'name' : "systems in a paper by Noonburg",
        'sub-menu' : [
            {
            'name' : "for linear-product start systems",
            'snippet' : ["from phcpy.examples import noon3", "f = noon3()", "for p in f: print p", "from phcpy.solver import solve", "sols = solve(f)", "print 'the number of solutions :', len(sols)"],
            }],
        },
        {
        'name' : "the cyclic n-roots problem",
        'sub-menu' : [
            {
            'name' : "for polyhedral homotopies",
            'snippet' : ["from phcpy.families import cyclic", "c5 = cyclic(5)", "for p in c5: print p", "from phcpy.solver import solve", "sols = solve(c5)", "print 'the number of solutions :', len(sols)"],
            }],
        }],
    },
    {
    'name' : 'Schubert calculus',
    'sub-menu' : [
        {
        'name' : "Pieri homotopies",
        'sub-menu' : [
            {
            'name' : "curves in the Grassmannian",
            'snippet' : ["from phcpy.schubert import pieri_root_count, random_complex_matrix", "(m,p,q) = (2,2,1)", "n = m*p + q*(m+p)", "r = pieri_root_count(m,p,q)", "L = [random_complex_matrix(m+p,m) for k in range(n)]", "points = random_complex_matrix(n,1)", "from phcpy.schubert import run_pieri_homotopies", "(f,fsols) = run_pieri_homotopies(m,p,q,L,points)", "print 'number of solutions :', len(fsols)"],
            }],
        },
        {
        'name' : "Littlewood-Richardson homotopies",
        'sub-menu' : [
            {
            'name' : "resolving Schubert conditions",
            'snippet' : ["from phcpy.schubert import resolve_schubert_conditions as rsc", "brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]", "rsc(6, 3, brackets)"],
            },
            {
            'name' : "solving a generic instance",
            'snippet' : ["brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]", "from phcpy.schubert import littlewood_richardson_homotopies as lrh", "(count, flags, sys, sols) = lrh(6, 3, brackets, verbose=False)", "print 'the root count :', count", "for sol in sols: print(sol)", "print 'the number of solutions :', len(sols)"],
            }],
        }],
    },
    {
    'name' : 'Newton polytopes',
    'sub-menu' : [
        {
        'name' : "convex hulls of lattice polytopes",
        'sub-menu' : [
            {
            'name' : "vertices and edge normals",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "from phcpy.polytopes import planar_convex_hull as pch", "points = rp(2, 7, -9, 9)", "for point in points: print point", "(vertices, normals) = pch(points)", "print 'the vertex points :', vertices", "print 'the edge normals :', normals"],
            },
            {
            'name' : "facets in 3-space",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "points = rp(3, 10, -9, 9)", "for point in points: print point", "from phcpy.polytopes import convex_hull as ch", "facets = ch(3, points)", "for facet in facets: print facet"],
            }],
        },
        {
        'name' : "mixed volumes",
        'sub-menu' : [
            {
            'name' : "volume of one random polytope",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "from phcpy.polytopes import mixed_volume as mv", "p1 = rp(3, 5, -9, 9)", "print p1", "mv([3], [p1])"],
            },
            {
            'name' : "mixed volume of two random polytopes",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "from phcpy.polytopes import mixed_volume as mv", "p1 = rp(3, 5, -9, 9); p2 = rp(3, 5, -9, 9)", "mv([2, 1],[p1, p2])", "mv([1, 2],[p1, p2])"],
            }],
        },
        {
        'name' : "solving binomial systems",
        'sub-menu' : [
            {
            'name' : "solution curves are maps",
            'snippet' : ["f = [ 'x**2*y - z*x;', 'x**2*z - y**2*x;' ]", "from phcpy.maps import solve_binomials", "maps = solve_binomials(3, f)", "for map in maps: print map"],
            }],
        },
        {
        'name' : "power series solutions",
        'sub-menu' : [
            {
            'name' : "intersecting the Viviani curve",
            'snippet' : ["plane = '(1-s)*y + s*(y-1);'", "vp0 = 'x^2 + y^2 + z^2 - 4;'", "vp1 = '(x-1)^2 + y^2 - 1;'", "vivplane = [plane, vp0, vp1]", "vivs0 = vivplane + ['s;']", "from phcpy.solver import solve", "sols = solve(vivs0, verbose=False)", "for sol in sols: print sol", "from phcpy.series import standard_newton_series", "sersols = standard_newton_series(vivplane, sols, verbose=False)", "for srs in sersols: print srs"],
            }],
        }],
    },
    {
    'name' : 'the extension module',
    'sub-menu' : [
        {
        'name' : "the module interface",
        'sub-menu' : [
            {
            'name' : "storing and loading a system",
            'snippet' : ["from phcpy.interface import store_standard_system, load_standard_system", "store_standard_system(['x^2 - 1/3;'])", "load_standard_system()"],
            }],
        },
        {
        'name' : "wrappers to the C interface",
        'sub-menu' : [
            {
            'name' : "reading and writing a system",
            'snippet' : ["from phcpy.phcpy2c2 import py2c_syscon_read_standard_system as readsys", "readsys()", "from phcpy.phcpy2c2 import py2c_syscon_write_standard_system as writesys", "writesys()"],
            }],
        }],
    }],
    }];
});

