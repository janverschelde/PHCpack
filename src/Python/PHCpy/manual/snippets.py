"""
The functions in this module export the code snippets of the menus
of the Jupyter notebook extension of PHCpy.
Every function is self contained and illustrates one particular feature.
"""
# I. snippets on the blackbox solver

def solve_random_trinomials():
    """
    Illustrates the solution of random trinomials.
    """
    from phcpy.solver import random_trinomials
    f = random_trinomials()
    for pol in f: print(pol)
    from phcpy.solver import solve
    sols = solve(f)
    for sol in sols: print(sol)
    print(len(sols), "solutions found")

def solve_specific_trinomials():
    """
    Illustrates the solution of specific trinomials.
    """
    f = ['x^2*y^2 + 2*x - 1;', 'x^2*y^2 - 3*y + 1;']
    from phcpy.solver import solve
    sols = solve(f)
    for sol in sols: print(sol)

def solution_from_string_to_dictionary():
    """
    Illustrates the dictionary format of a solution.
    """
    p = ['x + y - 1;', '2*x - 3*y + 1;']
    from phcpy.solver import solve
    sols = solve(p)
    print(sols[0])
    from phcpy.solutions import strsol2dict
    dsol = strsol2dict(sols[0])
    print(dsol.keys())
    for key in dsol.keys(): print('the value for', key, 'is', dsol[key])

def verify_by_evaluation():
    """
    Illustrates the verification of solutions by evaluation.
    """
    p = ['x + y - 1;', '2*x - 3*y + 1;']
    from phcpy.solver import solve
    sols = solve(p)
    from phcpy.solutions import strsol2dict, evaluate
    dsol = strsol2dict(sols[0])
    eva = evaluate(p, dsol)
    for val in eva: print(val)

def make_a_solution():
    """
    Illustrates the making of a solution.
    """
    from phcpy.solutions import make_solution
    s0 = make_solution(['x', 'y'], [float(3.14), complex(0, 2.7)])
    print(s0)
    s1 = make_solution(['x', 'y'], [int(2), int(3)])
    print(s1)

def filter_solution_lists():
    """
    Illustrates the filtering of solution lists.
    """
    from phcpy.solutions import make_solution, is_real, filter_real
    s0 = make_solution(['x', 'y'], [float(3.14), complex(0, 2.7)])
    print(is_real(s0, 1.0e-8))
    s1 = make_solution(['x', 'y'], [int(2), int(3)])
    print(is_real(s1, 1.0e-8))
    realsols = filter_real([s0, s1], 1.0e-8, 'select')
    for sol in realsols: print(sol)

def coordinates_names_values():
    """
    Illustrates coordinates, names, values of a solution.
    """
    from phcpy.solver import solve
    p = ['x^2*y^2 + x + 1;', 'x^2*y^2 + y + 1;']
    s = solve(p)
    print(s[0])
    from phcpy.solutions import coordinates, make_solution
    (names, values) = coordinates(s[0])
    print(names)
    print(values)
    s0 = make_solution(names, values)
    print(s0)

def fix_retrieve_seed():
    """
    Illustrates fixing and retrieving the seed.
    """
    from phcpy.dimension import set_seed, get_seed
    set_seed(2024)
    print(get_seed())


def solve_with_four_threads():
    """
    Ilustrates solving with four threads.
    """
    from phcpy.solver import solve
    from phcpy.families import cyclic
    nbrt = 4 # number of tasks
    pols = cyclic(6)
    print('solving the cyclic 6-roots problem :')
    for pol in pols: print(pol)
    from datetime import datetime
    starttime = datetime.now()
    sols = solve(pols)
    stoptime = datetime.now()
    elapsed = stoptime - starttime
    print(f'elapsed time with no multithreading :')
    print(elapsed)
    starttime = datetime.now()
    sols = solve(pols, tasks=nbrt)
    stoptime = datetime.now()
    elapsed = stoptime - starttime
    print(f'elapsed time with {nbrt} threads :')
    print(elapsed)

def four_root_counts():
    """
    Illustrates four root counts.
    """
    f = ['x^3*y^2 + x*y^2 + x^2;', 'x^5 + x^2*y^3 + y^2;']
    from phcpy.starters import total_degree
    print('the total degree :', total_degree(f))
    from phcpy.starters import m_homogeneous_bezout_number as mbz
    (bz, part) = mbz(f)
    print('a multihomogeneous Bezout number :', bz)
    from phcpy.starters import linear_product_root_count as lrc
    (rc, ssrc) = lrc(f)
    print('a linear-product root count :', rc)
    from phcpy.volumes import stable_mixed_volume
    (mv, smv) = stable_mixed_volume(f)
    print('the mixed volume :', mv)
    print('the stable mixed volume :', smv)

def newton_and_deflation():
    """
    Illustrates Newton's method and deflation.
    """
    p = ['(29/16)*x^3 - 2*x*y;', 'x^2 - y;']
    from phcpy.solutions import make_solution
    s = make_solution(['x', 'y'],[float(1.0e-6), float(1.0e-6)])
    print(s)
    from phcpy.deflation import double_newton_step
    s2 = double_newton_step(p, [s])
    print(s2[0])
    s3 = double_newton_step(p, s2)
    print(s3[0])
    from phcpy.deflation import double_deflate
    sd = double_deflate(p, [s])
    print(sd[0])

def overconstrained_deflation():
    """
    Illustrates deflation of an overconstrained system.
    """
    from phcpy.solutions import make_solution
    from phcpy.deflation import double_deflate
    sol = make_solution(['x', 'y'], [float(1.0e-6), float(1.0e-6)])
    print(sol)
    pols = ['x**2;', 'x*y;', 'y**2;']
    sols = double_deflate(pols, [sol], tolrnk=1.0e-8)
    print(sols[0])
    sols = double_deflate(pols, [sol], tolrnk=1.0e-4)
    print(sols[0])

def equation_and_variable_scaling():
    """
    Illustrates equation and variable scaling.
    """
    print('solving without scaling ...')
    from phcpy.solver import solve
    p = ['0.000001*x^2 + 0.000004*y^2 - 4;', '0.000002*y^2 - 0.001*x;']
    psols = solve(p)
    for sol in psols: print(sol)
    print('solving after scaling ...')
    from phcpy.scaling import double_scale_system as scalesys
    from phcpy.scaling import double_scale_solutions as scalesols
    (q, c) = scalesys(p)
    print('the scaled polynomial system :')
    for pol in q: print(pol)
    qsols = solve(q)
    ssols = scalesols(len(q), qsols, c)
    for sol in ssols: print(sol)

def blackbox_solver():
    """
    Runs the code snippets on the blackbox solver.
    """
    print("I. blackbox solver")
    print("I.1. solving trinomials")
    print("I.1.1 solving a random case")
    solve_random_trinomials()
    print("I.1.2. solving a specific case")
    solve_specific_trinomials()
    print("I.2. representations of isolated solutions")
    print("I.2.1 from string to dictionary")
    solution_from_string_to_dictionary()
    print("I.2.2 verify by evaluation")
    verify_by_evaluation()
    print("I.2.3 making a solution")
    make_a_solution()
    print("I.2.4 filtering solution lists")
    filter_solution_lists()
    print("I.2.5 coordinates, names, and values")
    coordinates_names_values()
    print("I.3 reproducible runs with fixed seed")
    print("I.3.1 fixing and retrieving the seed")
    fix_retrieve_seed()
    print("I.4 shared memory parallelism")
    print("I.4.1 solving with four threads")
    solve_with_four_threads()
    print("I.5 root counting methods")
    print("I.5.1 four different root counts")
    four_root_counts()
    print("I.6 Newton's method and deflation")
    print("I.6.1 the Griewank-Osborne example")
    newton_and_deflation()
    print("I.6.2 deflating an overconstrained system")
    overconstrained_deflation()
    print("I.7 equation and variable scaling")
    equation_and_variable_scaling()

# snippets on the path trackers

def total_degree_start_system():
    """
    The product of the degrees of the polynomials (the total degree)
    provide an upper bound on the number of solutions.
    A total degree start system is a simple system that has
    as many solutions as the product of the degrees.
    """
    from phcpy.starters import total_degree
    from phcpy.starters import total_degree_start_system
    from phcpy.trackers import double_track
    p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    d = total_degree(p)
    print('the total degree :', d)
    (q, qsols) = total_degree_start_system(p)
    print('the number of start solutions :', len(qsols))
    print('the start system :', q)
    s = double_track(p, q, qsols)
    print('the number of solutions :', len(s))
    for sol in s: print(sol)

def track_one_path():
    """
    Illustrates the tracking of one solution path.
    The order in which the solutions appear at the end
    depends on the gamma constant.
    """
    from phcpy.starters import total_degree_start_system
    from phcpy.trackers import double_track
    p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    (q, qsols) = total_degree_start_system(p)
    g1, s1 = double_track(p, q, [qsols[2]])
    print('first gamma :', g1)
    print(s1[0])
    g2, s2 = double_track(p, q, [qsols[2]])
    print('second gamma :', g2)
    print(s2[0])

def track_with_fixed_gamma():
    """
    Fixing the gamma in the homotopies fixes the order
    in which the solutions appear at the end of the paths.
    """
    from phcpy.starters import total_degree_start_system
    from phcpy.trackers import double_track
    p = ['x^2 + 4*y^2 - 4;', '2*y^2 - x;']
    (q, qsols) = total_degree_start_system(p)
    g3, s3 = double_track(p, q, [qsols[2]], \
        gamma=complex(0.824372806319,0.56604723848934))
    print('gamma :', g3)
    print('the solution at the end:')
    print(s3[0])

def get_next_point_on_path():
    """
    A step-by-step path tracker gives control to the user
    who can ask for the next point on a path.
    """
    from phcpy.starters import total_degree_start_system
    p = ['x**2 + 4*x**2 - 4;', '2*y**2 - x;']
    (q, s) = total_degree_start_system(p)
    from phcpy.trackers import initialize_double_tracker
    from phcpy.trackers import initialize_double_solution
    from phcpy.trackers import next_double_solution
    initialize_double_tracker(p, q)
    initialize_double_solution(len(p), s[0])
    s1 = next_double_solution()
    print('the next point on the solution path :')
    print(s1)
    print(next_double_solution())
    print(next_double_solution())
    initialize_double_solution(len(p), s[1])
    points = [next_double_solution() for i in range(11)]
    from phcpy.solutions import strsol2dict
    dicpts = [strsol2dict(sol) for sol in points]
    xvals = [sol['x'] for sol in dicpts]
    print('the x-coordinates on the path :')
    for x in xvals: print(x)

def plot_trajectories():
    """
    The step-by-step path tracker is applied to plot the
    trajectories of the solutions using matplotlib.
    """
    import matplotlib.pyplot as plt
    p = ['x^2 + y - 3;', 'x + 0.125*y^2 - 1.5;']
    print('constructing a total degree start system ...')
    from phcpy.starters import total_degree_start_system
    q, qsols = total_degree_start_system(p)
    print('number of start solutions :', len(qsols))
    from phcpy.trackers import initialize_double_tracker
    from phcpy.trackers import initialize_double_solution
    from phcpy.trackers import next_double_solution
    initialize_double_tracker(p, q, False)
    from phcpy.solutions import strsol2dict
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
            if(tval >= 1.0):
                break
        print(ns)
        xre = [point.real for point in xpoints]
        yre = [point.real for point in ypoints]
        axs.set_xlim(min(xre)-0.3, max(xre)+0.3)
        axs.set_ylim(min(yre)-0.3, max(yre)+0.3)
        dots, = axs.plot(xre,yre,'r-')
        fig.canvas.draw()
    fig.canvas.draw()
    ans = input('hit return to continue')

def polyhedral_homotopies():
    """
    Polyhedral homotopies solve random coefficient systems
    tracking an optimal number of solution paths,
    that is equal to the mixed volume.
    """
    from phcpy.volumes import mixed_volume
    from phcpy.volumes import double_polyhedral_homotopies
    from phcpy.trackers import double_track
    p = ['x^3*y^2 - 3*x^3 + 7;','x*y^3 + 6*y^3 - 9;']
    print('the mixed volume :', mixed_volume(p))
    (q, qsols) = double_polyhedral_homotopies()
    print('the number of start solutions :', len(qsols))
    gamma, psols = double_track(p, q, qsols)
    print('the number of solutions at the end :', len(psols))
    for sol in psols: print(sol)

def path_trackers():
    """
    Runs the snippets on the path trackers.
    """
    print("II path trackers")
    print("II.1 total degree start system")
    total_degree_start_system()
    print("II.2 track one path")
    track_one_path()
    print("II.3 track with fixed gamma")
    track_with_fixed_gamma()
    print("II.4 get next point on path")
    get_next_point_on_path()
    print("II.5 plot trajectories")
    plot_trajectories()
    print("II.6 polyhedral homotopies")
    polyhedral_homotopies()

# snippets on the sweep homotopies

def quadratic_turning_point():
    """
    Arc length parameter continuation is applied in a real sweep
    which ends at a quadratic turning point.
    """
    from phcpy.sweepers import double_real_sweep
    from phcpy.solutions import make_solution
    circle = ['x^2 + y^2 - 1;', 'y*(1-s) + (y-2)*s;']
    first = make_solution(['x', 'y', 's'], [1.0, 0.0, 0.0])
    second = make_solution(['x', 'y', 's'], [-1.0, 0.0, 0.0])
    startsols = [first, second]
    newsols = double_real_sweep(circle, startsols)
    for sol in newsols: print(sol)

def complex_parameter_homotopy_continuation():
    """
    Sweeps the parameter space with a convex linear combination
    of the parameters.  By a random gamma constant, no singularities
    are encountered during this complex sweep.
    """
    from phcpy.sweepers import double_complex_sweep
    from phcpy.solutions import make_solution
    circle = ['x^2 + y^2 - 1;']
    first = make_solution(['x', 'y'], [1.0, 0.0])
    second = make_solution(['x', 'y'], [-1.0, 0.0])
    startsols = [first, second]
    par = ['y']
    start = [0, 0]
    target = [2, 0]
    newsols = double_complex_sweep(circle, startsols, 2, par, start, target)
    for sol in newsols: print(sol)

def sweep_homotopies():
    """
    We distinguish between real and complex sweeps
    through the parameter space.
    """
    print("III sweep homotopies")
    print("III.1 computing a quadratic turning point")
    quadratic_turning_point()
    print("III.2 convex parameter homotopy continuation")
    complex_parameter_homotopy_continuation()

# snippets on schubert calculus

def lines_meeting_four_lines():
    """
    Applies Pieri homotopies to compute all lines meeting
    four given lines in 3-space.
    """
    from phcpy.schubert import pieri_root_count
    from phcpy.schubert import random_complex_matrix, run_pieri_homotopies
    (m, p, q) = (2, 2, 0)
    dim = m*p + q*(m+p)
    roco = pieri_root_count(m, p, q)
    print('the root count :', roco)
    L = [random_complex_matrix(m+p, m) for _ in range(dim)]
    (f, fsols) = run_pieri_homotopies(m, p, q, L)
    for sol in fsols: print(sol)
    print('number of solutions :', len(fsols))

def line_producing_interpolating_curves():
    """
    Applies Pieri homotopies to compute line producing curves,
    interpolating at given lines in 3-space.
    """
    from phcpy.schubert import pieri_root_count
    from phcpy.schubert import random_complex_matrix, run_pieri_homotopies
    (m, p, q) = (2, 2, 1)
    dim = m*p + q*(m+p)
    roco = pieri_root_count(m, p, q)
    print('the root count :', roco)
    L = [random_complex_matrix(m+p, m) for _ in range(dim)]
    points = random_complex_matrix(dim, 1)
    (f, fsols) = run_pieri_homotopies(m, p, q, L, 0, points)
    print('number of solutions :', len(fsols))

def resolve_some_schubert_conditions():
    """
    Resolves an example of a Schubert condition,
    applying the Littlewood-Richardson rule.
    """
    from phcpy.schubert import resolve_schubert_conditions
    brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    roco = resolve_schubert_conditions(6, 3, brackets)
    print('number of solutions :', roco)

def solve_generic_schubert_problem():
    """
    Runs the Littlewood-Richardson homotopies to solve
    a generic instance of a Schubert problem.
    """
    brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
    from phcpy.schubert import double_littlewood_richardson_homotopies as lrh
    (count, flags, sys, sols) = lrh(6, 3, brackets, verbose=False)
    print('the root count :', count)
    for sol in sols: print(sol)
    print('the number of solutions :', len(sols))

def schubert_calculus():
    """
    Runs the code snippets on the schubert module.
    """
    print("IV. schubert calculus")
    print("IV.1 Pieri homotopies")
    print("IV.1.1 lines meeting four given lines")
    lines_meeting_four_lines()
    print("IV.1.2 line producing interpolating curves")
    line_producing_interpolating_curves()
    print("IV.2 Littlewood-Richardson homotopies")
    print("IV.2.1 resolving Schubert conditions")
    resolve_some_schubert_conditions()
    print("IV.2.2 solving a generic Schubert problem")
    solve_generic_schubert_problem()

def example4pade():
    """
    The function f(z) = ((1 + 1/2*z)/(1 + 2*z))^(1/2) is
    a solution x(s) of (1-s)*(x^2 - 1) + s*(3*x^2 - 3/2) = 0
    """
    from phcpy.solutions import make_solution
    from phcpy.series import double_newton_at_point
    from phcpy.series import double_pade_approximants
    pol = ['(x^2 - 1)*(1-s) + (3*x^2 - 3/2)*s;']
    variables = ['x', 's']
    sol1 = make_solution(variables, [1, 0])
    sol2 = make_solution(variables, [-1, 0])
    sols = [sol1, sol2]
    print('start solutions :')
    for sol in sols: print(sol)
    srs = double_newton_at_point(pol, sols, idx=2)
    print('The series :')
    for ser in srs: print(ser)
    pad = double_pade_approximants(pol, sols, idx=2)
    print('the Pade approximants :')
    for app in pad: print(app)

def viviani_expansion(vrblvl=0):
    """
    Computes the power series expansion for the Viviani curve,
    from a natural parameter perspective, in double precision.
    """
    from phcpy.series import double_newton_at_series
    pols = [ '2*t^2 - x;', \
             'x^2 + y^2 + z^2 - 4;' , \
             '(x-1)^2 + y^2 - 1;']
    lser = [ '2*t^2;', '2*t;', '2;']
    nser = double_newton_at_series(pols, lser, maxdeg=12, nbr=8)
    variables = ['x', 'y', 'z']
    for (var, pol) in zip(variables, nser): print(var, '=', pol)

def apollonius_expansions():
    """
    Compare the series expansions at two solutions
    for the problem of Apollonius.
    """
    from phcpy.series import double_newton_at_series
    pols = [ 'x1^2 + 3*x2^2 - r^2 - 2*r - 1;', \
             'x1^2 + 3*x2^2 - r^2 - 4*x1 - 2*r + 3;', \
       '3*t^2 + x1^2 - 6*t*x2 + 3*x2^2 - r^2 + 6*t - 2*x1 - 6*x2 + 2*r + 3;']
    lser1 = ['1;', '1 + 0.536*t;', '1 + 0.904*t;']
    lser2 = ['1;', '1 + 7.464*t;', '1 + 11.196*t;']
    nser1 = double_newton_at_series(pols, lser1, idx=4, nbr=7)
    nser2 = double_newton_at_series(pols, lser2, idx=4, nbr=7)
    variables = ['x', 'y', 'z']
    print('the first solution series :')
    for (var, pol) in zip(variables, nser1): print(var, '=', pol)
    print('the second solution series :')
    for (var, pol) in zip(variables, nser2): print(var, '=', pol)

def power_series():
    """
    Runs the code snippets on the series module.
    """
    print("V. power series")
    print("V.1 example for Pade approximants")
    example4pade()
    print("V.2 series expansion for the Viviani curve")
    viviani_expansion()
    print("V.3 series expansions for the problem of Apollonius")
    apollonius_expansions()

def twisted_cubic_witness_set():
    """
    Makes a witness set for the twisted cubic.
    """
    twisted = ['x^2 - y;', 'x^3 - z;']
    from phcpy.sets import double_embed
    embtwist = double_embed(3, 1, twisted)
    print('embedded system augmented with a random hyperplane :')
    for pol in embtwist:
        print(pol)
    from phcpy.solver import solve
    sols = solve(embtwist)
    print('number of generic points :', len(sols))
    for sol in sols: print(sol)

def homotopy_membership_test():
    """
    Homotopies to decide whether a point belongs to a positive
    dimensional solution set.
    """
    from phcpy.families import cyclic
    c4 = cyclic(4)
    from phcpy.sets import double_embed
    c4e1 = double_embed(4, 1, c4)
    from phcpy.solver import solve
    sols = solve(c4e1)
    from phcpy.solutions import filter_zero_coordinates as filter
    genpts = filter(sols, 'zz1', 1.0e-8, 'select')
    print('generic points on the cyclic 4-roots set :')
    for sol in genpts: print(sol)
    from phcpy.sets import double_membertest
    pt0 = [1, 0, -1, 0, 1, 0, -1, 0]
    ismbr = double_membertest(c4e1, sols, 1, pt0)
    print('Is', pt0, 'a member?', ismbr)
    pt1 = [1, 0, 1, 0, -1, 0, -1, 0]
    ismbr = double_membertest(c4e1, sols, 1, pt1)
    print('Is', pt1, 'a member?', ismbr)

def factor_a_cubic():
    """
    Illustrates the numerical factorization of a cubic.
    """
    cubic = '(x+1)*(x^2 + y^2 + 1);'
    from phcpy.sets import double_hypersurface_set
    (wit, pts) = double_hypersurface_set(2, cubic)
    for pol in wit:
        print(pol)
    print('number of witness points :', len(pts))
    for (idx, sol) in enumerate(pts):
        print('Solution', idx+1, ':')
        print(sol)
    from phcpy.factor import double_monodromy_breakup, write_factorization
    deco = double_monodromy_breakup(wit, pts, dim=1)
    write_factorization(deco)

def cascades_of_homotopies():
    """
    A cascade of homotopies computes generic points on all positive
    components of the solution set, for all dimensions.
    """
    pol1 = '(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);'
    pol2 = '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);'
    pol3 = '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);'
    pols = [pol1, pol2, pol3]
    from phcpy.cascades import double_top_cascade, double_cascade_filter
    (embpols, sols0, sols1) = double_top_cascade(3, 2, pols)
    print('at dimension 2, degree :', len(sols0))
    (wp1, ws0, ws1) = double_cascade_filter(2, embpols, sols1, tol=1.0e-8)
    print('at dimension 1, candidate generic points :', len(ws0))
    (wp0, ws0, ws1) = double_cascade_filter(1, wp1, ws1, tol=1.0e-8)
    print('candidate isolated points :', len(ws0))

def numerical_irreducible_decomposition():
    """
    Runs a test example on a numerical irreducible decomposition.
    """
    pol0 = '(x1-1)*(x1-2)*(x1-3)*(x1-4);'
    pol1 = '(x1-1)*(x2-1)*(x2-2)*(x2-3);'
    pol2 = '(x1-1)*(x1-2)*(x3-1)*(x3-2);'
    pol3 = '(x1-1)*(x2-1)*(x3-1)*(x4-1);'
    pols = [pol0, pol1, pol2, pol3]
    from phcpy.decomposition import solve, write_decomposition
    deco = solve(pols, verbose=False)
    write_decomposition(deco)

def sphere_cylinder_intersection():
    """
    A sphere intersected by a cylinder defines a quartic curve.
    """
    sphere = 'X^2 + Y^2 + Z^2 - 1;'
    cylinder = 'X^2 + 1.0e-14*Y^2 + (Z - 0.5)^2 - 1;'
    from phcpy.sets import double_hypersurface_set
    (spheqs, sphpts) = double_hypersurface_set(3, sphere)
    (cyleqs, cylpts) = double_hypersurface_set(3, cylinder)
    from phcpy.diagonal import double_diagonal_solve
    quaeqs, quapts = double_diagonal_solve\
        (3, 2, spheqs, sphpts, 2, cyleqs, cylpts)
    for pol in quaeqs: print(pol)
    for sol in quapts: print(sol)

def positive_dimensional_sets():
    """
    Runs the code snippets on the methods to compute
    a numerical irreducible decomposition.
    """
    print("VI. positive dimensional solution sets")
    print("VI.1 witnes set for twisted cubic")
    twisted_cubic_witness_set()
    print("VI.2 homotopy membership test")
    homotopy_membership_test()
    print("VI.3 cascades of homotopies")
    cascades_of_homotopies()
    print("VI.4 factor a cubic polynomial")
    factor_a_cubic()
    print("VI.5 numerical irreducible decomposition")
    numerical_irreducible_decomposition()
    print("VI.6 diagonal homotopies")
    sphere_cylinder_intersection()

def main():
    """
    Runs all code snippets.
    """
    blackbox_solver()
    path_trackers()
    sweep_homotopies()
    schubert_calculus()
    power_series()
    positive_dimensional_sets()

if __name__=='__main__':
    main()
