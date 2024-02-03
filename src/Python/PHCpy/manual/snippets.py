"""
The functions in this module export the code snippets of the menus
of the Jupyter notebook extension of PHCpy.
Every function is self contained and illustrates one particular feature.
"""
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
    print("II. schubert calculus")
    print("II.1 Pieri homotopies")
    print("II.1.1 lines meeting four given lines")
    lines_meeting_four_lines()
    print("II.1.2 line producing interpolating curves")
    line_producing_interpolating_curves()
    print("II.2 Littlewood-Richardson homotopies")
    print("II.2.1 resolving Schubert conditions")
    resolve_some_schubert_conditions()
    print("II.2.2 solving a generic Schubert problem")
    solve_generic_schubert_problem()

def main():
    """
    Runs all code snippets.
    """
    blackbox_solver()
    schubert_calculus()

if __name__=='__main__':
    main()
