"""
Illustrative example for a numerical irreducible decomposition.
The script illustrates a two-stage cascade to compute candidate
generic points on all components, on all dimensions of the solution set.
"""
pols = ['(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);', \
        '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);', \
        '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);']
from phcpy.cascades import top_cascade
(topemb, topsols0, topsols1) = top_cascade(3, 2, pols, 1.0e-8)
from phcpy.solutions import filter_zero_coordinates as filter
print 'generic points on the two dimensional surface :'
for sol in topsols0:
    print sol
raw_input('hit enter to continue')
"""
Run the cascade with the nonsolutions.
"""
from phcpy.cascades import cascade_step
print 'running a cascade with %d paths ...' % len(topsols1)
lvl1sols = cascade_step(2, topemb, topsols1)
"""
After the filtering, we must drop variables, coordinates,
and hyperplane for the next level in the cascade.
"""
from phcpy.sets import drop_variable_from_standard_polynomials as drop1poly
from phcpy.sets import drop_coordinate_from_standard_solutions as drop1sols
lvl1emb = drop1poly(topemb, 'zz2')
lvl1solsdrop = drop1sols(lvl1sols, len(topemb), 'zz2')
lvl1sols0 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'select') 
lvl1sols1 = filter(lvl1solsdrop, 'zz1', 1.0e-8, 'remove') 
from phcpy.solutions import filter_regular as regfilt
reglvl1sols0 = regfilt(lvl1sols0, 1.0e-8, 'select')
for sol in reglvl1sols0:
    print sol
raw_input('hit enter to continue')
lvl1emb = lvl1emb[:-1]  # dropping the last polynomial
print 'running a cascade with %d paths ...' % len(lvl1sols1)
lvl2sols = cascade_step(1, lvl1emb, lvl1sols1)
print 'finished the cascade'
lvl2solsdrop = drop1sols(lvl2sols, len(lvl1emb), 'zz1')
from phcpy.solutions import filter_regular as regfilt
reglvl2sols = regfilt(lvl2solsdrop, 1.0e-8, 'select')
print 'the solutions at the end of the cascade :'
for sol in reglvl2sols:
    print sol
