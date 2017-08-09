"""
Illustrative example for a numerical irreducible decomposition.
This python3 script illustrates a two-stage cascade to compute candidate
generic points on all components, on all dimensions of the solution set.
"""
pols = ['(x^2 + y^2 + z^2 - 1)*(y - x^2)*(x - 0.5);', \
        '(x^2 + y^2 + z^2 - 1)*(z - x^3)*(y - 0.5);', \
        '(x^2 + y^2 + z^2 - 1)*(z - x*y)*(z - 0.5);']
from phcpy.cascades import top_cascade, cascade_filter
(topemb, topsols0, topsols1) = top_cascade(3, 2, pols, 1.0e-8)
print('generic points on the two dimensional surface :')
for sol in topsols0:
    print(sol)
input('hit enter to continue')
(lvl1emb, lvl1sols0, lvl1sols1) = cascade_filter(2, topemb, topsols1, 1.0e-8)
print('candidate generic points at level 1 :')
for sol in lvl1sols0:
    print(sol)
input('hit enter to continue')
(lvl0emb, lvl2sols) = cascade_filter(1, lvl1emb, lvl1sols1, 1.0e-8)
print('finished the cascade')
print('the solutions at the end of the cascade :')
for sol in lvl2sols:
    print(sol)
