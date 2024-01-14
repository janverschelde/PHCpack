"""
The 2-by-2 minors of a 2-by-n matrix of indeterminates
has an interesting numerical irreducible decomposition.
"""
from phcpy.families import adjacent_minors
from phcpy.dimension import get_core_count
from phcpy.solver import solve
from phcpy.sets import double_embed
from phcpy.factor import double_monodromy_breakup
from phcpy.factor import write_decomposition

for col in range(3, 6):
    nvr = 2*col
    neq = col - 1
    dim = nvr - neq
    print('  number of columns :', col)
    print('number of variables :', nvr)
    print('number of equations :', neq)
    print('dimension of solution set :', dim)
    minors = adjacent_minors(2, col)
    print(f'the adjacent 2-by-2 minors of a 2-by-{col} matrix :')
    for pol in minors:
        print(pol)
    eminors = double_embed(nvr, dim, minors)
    esols = solve(eminors)
    print('the degree :', len(esols))
    deco = double_monodromy_breakup(eminors, esols, nvr-neq)
    print('the decomposition :')
    write_decomposition(deco)
    print('the number of factors :', len(deco))
    
