Code Snippets
=============

The functions are based on the code snippets of the menus
of the Jupyter notebook extension of PHCpy.
Every function is self contained and illustrates one particular feature.

the blackbox solver
-------------------

::

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

::

    def solve_specific_trinomials():
        """
        Illustrates the solution of specific trinomials.
        """
        f = ['x^2*y^2 + 2*x - 1;', 'x^2*y^2 - 3*y + 1;']
        from phcpy.solver import solve
        sols = solve(f)
        for sol in sols: print(sol)

::

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

::

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

::

    def make_a_solution():
        """
        Illustrates the making of a solution.
        """
        from phcpy.solutions import make_solution
        s0 = make_solution(['x', 'y'], [float(3.14), complex(0, 2.7)])
        print(s0)
        s1 = make_solution(['x', 'y'], [int(2), int(3)])
        print(s1)

::

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

::

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

::

    def fix_retrieve_seed():
        """
        Illustrates fixing and retrieving the seed.
        """
        from phcpy.dimension import set_seed, get_seed
        set_seed(2024)
        print(get_seed())

::

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

::

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

::

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

::

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

::

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

path trackers
-------------

::

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

::

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

::

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

::

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

::

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

::

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

sweep homotopies
----------------

::

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

::

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

schubert calculus
-----------------

::

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

::

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

::

    def resolve_some_schubert_conditions():
        """
        Resolves an example of a Schubert condition,
        applying the Littlewood-Richardson rule.
        """
        from phcpy.schubert import resolve_schubert_conditions
        brackets = [[2, 4, 6], [2, 4, 6], [2, 4, 6]]
        roco = resolve_schubert_conditions(6, 3, brackets)
        print('number of solutions :', roco)

::

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
