
require(["nbextensions/snippets_menu/main"], function (snippets_menu) {
    snippets_menu.options['menus'] = [{
        'name' : 'PHCpy',
        'sub-menu-direction' : 'left',
        'sub-menu' : [
    {
    'name' : 'blackbox',
    'sub-menu' : [
        {
        'name' : "solving trinomials",
        'sub-menu' : [
            {
            'name' : "solving a random case",
            'snippet' : ["from phcpy.solver import random_trinomials", "f = random_trinomials()", "for pol in f: print pol", "from phcpy.solver import solve", "sols = solve(f, silent=True)", "for sol in sols: print sol", "print len(sols), \"solutions found\""],
            },
            {
            'name' : "solving a specific case",
            'snippet' : ["f = ['x*y^2 + y - 3;', 'x^3 - y + 1;']", "from phcpy.solver import solve", "sols = solve(f)", "for sol in sols: print sol"],
            }],
        },
        {
        'name' : "representations of solutions of polynomial systems ",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["print(s[0])", "from phcpy.solutions import strsol2dict", "d = strsol2dict(s[0])", "d.keys()", "d['x1']"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.solutions import evaluate", "e = evaluate(f,d)", "for x in e: print(x)"],
            },
            {
            'name' : "2",
            'snippet' : ["from phcpy.solutions import make_solution", "s0 = make_solution(['x', 'y'], [1, complex(0, 1)])", "print(s0)", "s1 = make_solution(['x', 'y'], [2, 3])", "print(s1)"],
            },
            {
            'name' : "3",
            'snippet' : ["from phcpy.solutions import is_real, filter_real", "is_real(s0, 1.0e-8)", "is_real(s1, 1.0e-8)", "realsols = filter_real([s0, s1], 1.0e-8, 'select')", "for sol in realsols: print(sol)"],
            },
            {
            'name' : "4",
            'snippet' : ["from phcpy.solver import solve", "p = ['x**2 - 3*y + 1;', 'x*y - 3;']", "s = solve(p, silent=True)", "print(s[0])"],
            },
            {
            'name' : "5",
            'snippet' : ["from phcpy.solutions import coordinates, make_solution", "(names, values) = coordinates(s[0])", "names", "values", "s0 = make_solution(names, values)", "print(s0)"],
            }],
        },
        {
        'name' : "reproducible runs with fixed seeds",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.phcpy2c import py2c_set_seed", "py2c_set_seed(2013)"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.phcpy2c import py2c_get_seed", "py2c_get_seed()"],
            }],
        },
        {
        'name' : "root counting methods",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["f = ['x^3*y^2 + x*y^2 + x^2;', 'x^5 + x^2*y^3 + y^2;']", "from phcpy.solver import total_degree", "total_degree(f)", "from phcpy.solver import m_homogeneous_bezout_number as mbz", "mbz(f)", "from phcpy.solver import linear_product_root_count as lrc", "lrc(f)", "from phcpy.solver import mixed_volume", "mixed_volume(f, stable=True)"],
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
            'name' : "0",
            'snippet' : ["from phcpy.solver import solve", "p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']", "psols = solve(p, silent=True)", "print(psols[0])"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.solver import standard_scale_system as scalesys", "from phcpy.solver import standard_scale_solutions as scalesols", "(q, c) = scalesys(p)", "q[0]", "q[1]"],
            },
            {
            'name' : "2",
            'snippet' : ["qsols = solve(q, silent=True)", "ssols = scalesols(len(q), qsols, c)", "for sol in ssols: print(sol)"],
            }],
        }],
    },

    {
    'name' : 'pathtrack',
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
            'name' : "0",
            'snippet' : ["s3 = track(p, q, [qsols[2]], gamma=complex(0.824372806319,0.56604723848934))", "print(s3[0])"],
            }],
        },
        {
        'name' : "give the next solution on a path",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.solver import total_degree_start_system", "p = ['x**2 + 4*x**2 - 4;', '2*y**2 - x;']", "(q, s) = total_degree_start_system(p)", "from phcpy.trackers import initialize_standard_tracker", "from phcpy.trackers import initialize_standard_solution", "from phcpy.trackers import next_standard_solution", "initialize_standard_tracker(p, q)", "initialize_standard_solution(len(p), s[0])", "s1 = next_standard_solution()", "print(s1)", "print(next_standard_solution())", "print(next_standard_solution())"],
            },
            {
            'name' : "1",
            'snippet' : ["initialize_standard_solution(len(p), s[1])", "points = [next_standard_solution() for i in range(11)]", "from phcpy.solutions import strsol2dict", "dicpts = [strsol2dict(sol) for sol in points]", "xvals = [sol['x'] for sol in dicpts]", "for x in xvals: print(x)"],
            },
            {
            'name' : "2",
            'snippet' : ["p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']", "print('constructing a total degree start system ...')", "from phcpy.solver import total_degree_start_system as tds", "q, qsols = tds(p)", "print('number of start solutions :', len(qsols))", "from phcpy.trackers import initialize_standard_tracker", "from phcpy.trackers import initialize_standard_solution", "from phcpy.trackers import next_standard_solution", "initialize_standard_tracker(p, q, False)", "from phcpy.solutions import strsol2dict", "import matplotlib.pyplot as plt", "plt.ion()", "fig = plt.figure()", "for k in range(len(qsols)):", "    if(k == 0):", "        axs = fig.add_subplot(221)", "    elif(k == 1):", "        axs = fig.add_subplot(222)", "    elif(k == 2):", "        axs = fig.add_subplot(223)", "    elif(k == 3):", "        axs = fig.add_subplot(224)", "    startsol = qsols[k]", "    initialize_standard_solution(len(p),startsol)", "    dictsol = strsol2dict(startsol)", "    xpoints =  [dictsol['x']]", "    ypoints =  [dictsol['y']]", "    for k in range(300):", "        ns = next_standard_solution()", "        dictsol = strsol2dict(ns)", "        xpoints.append(dictsol['x'])", "        ypoints.append(dictsol['y'])", "        tval = eval(dictsol['t'].lstrip().split(' ')[0])", "        if(tval == 1.0):", "            break", "    print(ns)", "    xre = [point.real for point in xpoints]", "    yre = [point.real for point in ypoints]", "    axs.set_xlim(min(xre)-0.3, max(xre)+0.3)", "    axs.set_ylim(min(yre)-0.3, max(yre)+0.3)", "    dots, = axs.plot(xre,yre,'r-')", "    fig.canvas.draw()", "fig.canvas.draw()", "ans = raw_input('hit return to exit')"],
            }],
        },
        {
        'name' : "solving with polyhedral homotopies",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["p = ['x^3*y^2 - 3*x^3 + 7;','x*y^3 + 6*y^3 - 9;']", "from phcpy.solver import mixed_volume", "mixed_volume(p)", "from phcpy.solver import random_coefficient_system", "(q,qsols) = random_coefficient_system(silent=True)", "len(qsols)", "from phcpy.trackers import track", "psols = track(p,q,qsols)", "len(psols)", "print(psols[4])"],
            }],
        },
        {
        'name' : "Newton's method at higher precision",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["psols_dd = newton_step(p,psols,precision='dd')"],
            },
            {
            'name' : "1",
            'snippet' : ["p = ['x*y^3 + y - 2;', 'x^3*y + x - 8;']", "from phcpy.solver import linear_product_root_count", "r = linear_product_root_count(p)", "from phcpy.solver import random_linear_product_system", "(q,qsols) = random_linear_product_system(p)", "len(qsols)", "from phcpy.trackers import track", "psols = track(p,q,qsols)", "len(psols)", "from phcpy.solver import newton_step", "psols_dd = newton_step(p,psols,precision='dd')"],
            }],
        },
        {
        'name' : "multitasked path tracking",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["$ time python trackcyclic7.py", "number of start solutions : 924", "starting the path tracking with 1 task(s) ...", "tracked 924 solution paths"],
            },
            {
            'name' : "1",
            'snippet' : ["from sys import argv", "if(len(argv) == 1):", "    nbtasks = 1", "else:", "    nbtasks = eval(argv[1])", "from phcpy.phcpy2c import py2c_read_standard_target_system_from_file", "from phcpy.phcpy2c import py2c_read_standard_start_system_from_file", "from phcpy.phcpy2c import py2c_copy_target_system_to_container", "from phcpy.phcpy2c import py2c_copy_start_system_to_container", "from phcpy.phcpy2c import py2c_copy_start_solutions_to_container", "from phcpy.phcpy2c import py2c_solcon_number_of_solutions", "from phcpy.solver import load_standard_system, load_standard_solutions", "from phcpy.trackers import standard_double_track", "cyclic7 = '/Users/jan/PHCv2/Demo/cyclic7'", "cyclic7q = '/Users/jan/PHCv2/Demo/cyclic7q'", "fail = py2c_read_standard_target_system_from_file(len(cyclic7),cyclic7)", "fail = py2c_copy_target_system_to_container()", "target = load_standard_system()", "fail = py2c_read_standard_start_system_from_file(len(cyclic7q),cyclic7q)", "fail = py2c_copy_start_system_to_container()", "start = load_standard_system()", "fail = py2c_copy_start_solutions_to_container()", "sols = load_standard_solutions()", "print('number of start solutions :', py2c_solcon_number_of_solutions())", "print('starting the path tracking with', nbtasks, 'task(s) ...')", "endsols = standard_double_track(target, start, sols, 0, nbtasks)", "print('tracked', len(endsols), 'solution paths')"],
            }],
        },
        {
        'name' : "GPU accelerated path tracking",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["GPU = 1 # use the GPU", "DIR = '/home/jan/Problems/GPUdata/MultiPath' # location of systems", "from phcpy.phcpy2c \\", "import py2c_read_standard_target_system_from_file as read_target", "from phcpy.phcpy2c \\", "import py2c_read_standard_start_system_from_file as read_start", "cyc10tarfile = DIR + '/cyclic10.target'", "cyc10stafile = DIR + '/cyclic10.start'", "fail = read_target(len(cyc10tarfile), cyc10tarfile)", "from phcpy.interface import load_standard_system as loadsys", "from phcpy.interface import load_standard_solutions as loadsols", "cyc10 = loadsys()", "print('the cyclic 10-roots problem :')", "for pol in cyc10:", "    print(pol)", "fail = read_start(len(cyc10stafile), cyc10stafile)", "cyc10q = loadsys()", "print('a start system for the cyclic 10-roots problem :')", "for pol in cyc10q:", "    print(pol)", "cyc10qsols = loadsols()", "print('number of start solutions :', len(cyc10qsols))", "print('the first solution :')", "print(cyc10qsols[0])", "print('calling the path tracker...')", "if(GPU == 0):", "    from phcpy.trackers import ade_double_track", "    cyc10sols = ade_double_track(cyc10,cyc10q,cyc10qsols,verbose=0)", "else:", "    from phcpy.trackers import gpu_double_track", "    cyc10sols = gpu_double_track(cyc10,cyc10q,cyc10qsols,verbose=0)", "print('number of solutions :', len(cyc10sols))", "for sol in cyc10sols:", "    print(sol)"],
            }],
        },
        {
        'name' : "sweep homotopies",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.solutions import make_solution as makesol", "first = makesol(['x', 'y', 's'], [1, 0, 0])", "second = makesol(['x', 'y', 's'], [-1, 0, 0])", "startsols = [first, second]", "from phcpy.sweepers import standard_real_sweep as sweep", "newsols = sweep(circle, startsols)", "print(newsols[0])"],
            },
            {
            'name' : "2",
            'snippet' : ["t :  0.00000000000000E+00   0.00000000000000E+00", "m : 1", "the solution for t :", " x : -2.46519032881566E-32   0.00000000000000E+00", " y :  1.00000000000000E+00   0.00000000000000E+00", " s :  5.00000000000000E-01   0.00000000000000E+00", "== err :  0.000E+00 = rco :  1.000E+00 = res :  0.000E+00 ="],
            }],
        },
        {
        'name' : "real versus complex sweeps",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["circle = ['x^2 + y^2 - 1;']", "from phcpy.solutions import make_solution as makesol", "first = makesol(['x', 'y'], [1, 0])", "second = makesol(['x', 'y'], [-1, 0])", "startsols = [first, second]", "par = ['y']", "start = [0, 0] ", "target = [2, 0]", "from phcpy.sweepers import standard_complex_sweep as sweep", "newsols = sweep(circle, startsols, 2, par, start, target)"],
            },
            {
            'name' : "1",
            'snippet' : ["print(newsols[0])"],
            }],
        },
        {
        'name' : "tuning parameters, settings, and tolerances",
        'sub-menu' : [],
        },
        {
        'name' : "a polyhedral end game",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.tuning import order_endgame_extrapolator_set as set", "set(4)"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.tuning import order_endgame_extrapolator_get as get", "get()"],
            },
            {
            'name' : "2",
            'snippet' : ["f = ['x + y^3 - 1;', 'x + y^3 + 1;']", "from phcpy.solver import mixed_volume as mv", "from phcpy.solver import random_coefficient_system as rcs", "mv(f)", "(g, gsols) = rcs(f)", "len(gsols)"],
            },
            {
            'name' : "3",
            'snippet' : ["from phcpy.trackers import standard_double_track as track", "sols = track(f, g, gsols)", "from phcpy.tropisms import standard_retrieve as retrieve", "(w, d, e) = retrieve(len(sols), len(f))", "w"],
            }],
        }],
    },

    {
    'name' : 'posdimsols',
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
            'snippet' : ["twisted = ['x^2 - y;', 'x^3 - z;']", "from phcpy.sets import embed", "e = embed(3,1,twisted)", "from phcpy.solver import solve", "s = solve(e, silent=True)", "print 'number of generic points :', len(s)", "for sol in s: print sol"],
            }],
        },
        {
        'name' : "homotopy membership test",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.families import cyclic", "c4 = cyclic(4)", "from phcpy.sets import embed", "c4e1 = embed(4, 1, c4)", "from phcpy.solver import solve", "sols = solve(c4e1)", "from phcpy.solutions import filter_zero_coordinates as filter", "genpts = filter(sols, 'zz1', 1.0e-8, 'select')", "for sol in genpts:"],
            },
            {
            'name' : "1",
            'snippet' : ["point = [1, 0, -1, 0, 1, 0, -1, 0]", "from phcpy.sets import membertest", "membertest(c4e1, genpts, 1, point)"],
            },
            {
            'name' : "2",
            'snippet' : ["point = [-1, 0, -1, 0, 1, 0, 1, 0]", "membertest(c4e1, genpts, 1, point)"],
            },
            {
            'name' : "3",
            'snippet' : ["ddpoint = [-1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]"],
            },
            {
            'name' : "4",
            'snippet' : ["membertest(c4e1, genpts, 1, ddpoint, memtol=1.e-12, precision='dd')"],
            }],
        },
        {
        'name' : "cascade of homotopies",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["pols = ['(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);', \\"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.sets import embed", "topemb = embed(3, 2, pols)", "from phcpy.solver import solve", "topsols = solve(topemb, silent=True)"],
            },
            {
            'name' : "2",
            'snippet' : ["from phcpy.solutions import filter_zero_coordinates as filter", "topsols0 = filter(topsols, 'zz2', 1.0e-8, 'select')", "topsols1 = filter(topsols, 'zz2', 1.0e-8, 'remove')", "print('generic points on the two dimensional surface :')", "for sol in topsols0:"],
            },
            {
            'name' : "3",
            'snippet' : ["from phcpy.sets import cascade_step", "lvl1sols = cascade_step(topemb, topsols1)"],
            },
            {
            'name' : "4",
            'snippet' : ["from phcpy.sets import drop_variable_from_polynomials as drop1poly", "from phcpy.sets import drop_coordinate_from_solutions as drop1sols", "lvl1emb = drop1poly(topemb, 'zz2')", "lvl1emb = lvl1emb[:-1]  # dropping the last polynomial", "lvl1solsdrop = drop1sols(lvl1sols, len(topemb), 'zz2')", "lvl1sols0 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'select') ", "lvl1sols1 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'remove') "],
            },
            {
            'name' : "5",
            'snippet' : ["from phcpy.solutions import filter_regular as regfilt", "reglvl1sols0 = regfilt(lvl1sols0, 1.0e-8, 'select')", "for sol in reglvl1sols0:"],
            },
            {
            'name' : "6",
            'snippet' : ["lvl2sols = cascade_step(lvl1emb, lvl1sols1)", "lvl2solsdrop = drop1sols(lvl2sols, len(lvl1emb), 'zz1')", "for sol in reglvl2solsdrop:"],
            }],
        },
        {
        'name' : "factoring into irreducibles",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["p = '(x+1)*(x^2 + y^2 + 1);'"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.sets import witness_set_of_hypersurface as wh", "(w, s) = wh(2, p)", "len(s)"],
            },
            {
            'name' : "2",
            'snippet' : ["from phcpy.sets import factor", "f = factor(1, w, s)", "f"],
            },
            {
            'name' : "3",
            'snippet' : ["[([1, 2], 8.537360146292391e-15), ([3], 2.1316282072803006e-14)]"],
            },
            {
            'name' : "4",
            'snippet' : ["f = factor(1, w, s, precision='dd')", "f = factor(1, w, s, precision='qd')"],
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
            'name' : "0",
            'snippet' : ["sph = 'x^2 + y^2 + z^2 - 1;'", "cyl = 'x^2 + y - y + (z - 0.5)^2 - 1;'"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.sets import witness_set_of_hypersurface as witsurf", "sphwit = witsurf(3, sph)", "spheqs, sphpts = sphwit", "cylwit = witsurf(3, cyl)", "cyleqs, cylpts = cylwit"],
            },
            {
            'name' : "2",
            'snippet' : ["from phcpy.sets import diagonal_solver as diagsolve", "quawit = diagsolve(3, 2, spheqs, sphpts, 2, cyleqs, cylpts)", "quaeqs, quapts = quawit", "for pol in quaeqs:", "for sol in quapts:"],
            }],
        }],
    },

    {
    'name' : 'examfams',
    'sub-menu' : [
        {
        'name' : "interactive regression testing",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.examples import noon3", "f = noon3()", "for p in f: print(p)"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.solver import solve", "s = solve(f,silent=True)", "len(s)", "print(s[0])"],
            }],
        },
        {
        'name' : "the cyclic n-roots problem",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.families import cyclic", "c4 = cyclic(4)", "for p in c4: print(p)"],
            }],
        }],
    },

    {
    'name' : 'numschub',
    'sub-menu' : [
        {
        'name' : "Pieri homotopies",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.schubert import *", "(m,p,q) = (2,2,1)", "n = m*p + q*(m+p)", "r = Pieri_root_count(m,p,q)", "L = [random_complex_matrix(m+p,m) for k in range(n)]", "points = random_complex_matrix(n,1)", "(f,fsols) = run_Pieri_homotopies(m,p,q,L,points)"],
            }],
        },
        {
        'name' : "Littlewood-Richardson homotopies ",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.schubert import resolve_schubert_conditions as rsc", "brackets = [[2,4,6],[2,4,6],[2,4,6]]", "rsc(6,3,brackets)"],
            },
            {
            'name' : "1",
            'snippet' : ["from phcpy.schubert import littlewood_richardson_homotopies as lrh", "(count, flags, sys, sols) = lrh(6, 3, brackets, verbose=False)", "count", "for sol in sols: print(sol)", " len(sys)"],
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
            'name' : "0",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "from phcpy.polytopes import planar_convex_hull as pch", "points = rp(2,7,-9,9)", "points", "(vertices, normals) = pch(points)", "vertices", "normals"],
            },
            {
            'name' : "1",
            'snippet' : [" 9 - 2 \\times 8 = 5 - 2 \\times 6 = -7 "],
            },
            {
            'name' : "2",
            'snippet' : [" x_1 - 2 x_2 \\geq -7"],
            },
            {
            'name' : "3",
            'snippet' : ["from phcpy.polytopes import random_points as rp", "points = rp(3,10,-9,9)", "for point in points: print(point)", "from phcpy.polytopes import convex_hull as ch", "facets = ch(3, points)", "for facet in facets: print(facet)"],
            },
            {
            'name' : "4",
            'snippet' : ["90 x_1 - 65 x_2 - 6 x_3 \\geq -597"],
            },
            {
            'name' : "5",
            'snippet' : ["vertices = []", "for facet in facets:", "vertices", "len(vertices)"],
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
            'name' : "0",
            'snippet' : ["vivplane = ['(1-s)*y + s*(y-1);',", "vivs0 = vivplane + ['s;']", "from phcpy.solver import solve", "sols = solve(vivs0, silent=True)", "print(sols[0])", ""],
            },
            {
            'name' : "1",
            'snippet' : ["from series import standard_newton_series", "sersols = standard_newton_series(vivplane, sols, verbose=False)", "sersols[0]", ""],
            }],
        }],
    },

    {
    'name' : 'modphcpy2c3',
    'sub-menu' : [
        {
        'name' : "design of the Python to C interface",
        'sub-menu' : [],
        },
        {
        'name' : "the interface to PHCpack",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.interface import store_standard_system, load_standard_system", "store_standard_system(['x^2 - 1/3;'])", "load_standard_system()"],
            }],
        },
        {
        'name' : "wrappers to the C interface to PHCpack",
        'sub-menu' : [
            {
            'name' : "0",
            'snippet' : ["from phcpy.phcpy2c3 import py2c_syscon_read_standard_system as readsys", "from phcpy.phcpy2c3 import py2c_syscon_write_standard_system as writesys", "readsys()"],
            }],
        }],
    }],
    }];
});

