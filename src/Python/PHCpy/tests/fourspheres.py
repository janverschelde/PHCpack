"""
Computes all tangent lines to four mutually touching spheres.

The original formulation as polynomial system came from
Cassiano Durand, then at the CS department in Purdue.
The positioning of the centers of the spheres, each with radius
0.5 at the vertices of a tetrahedron came from Thorsten Theobald,
then at TU Muenich.  The centers of the four spheres are

 c1 = (0, 0, 0)
 c2 = (1, 0, 0);
 c3 = (1/2, sqrt(3)/2, 0)
 c4 = (1/2, sqrt(3)/6, sqrt(6)/3);   

Tangent vector t = (x0,x1,x2) and moment vector m = (x3,x4,x5).
The first equation is ||t||=1, the second m.t = 0,
the other equations are ||m - c_i x t ||^2 - r^2 = 0, where
the radius r = 1/2.
"""
from sympy import var, sqrt
from sympy.vector import CoordSys3D, Vector
from phcpy.solver import solve

N = CoordSys3D('N')
x0, x1, x2 = var('x0, x1, x2')
# tangent vector
vt = Vector.zero + x0*N.i + x1*N.j + x2*N.k
normt = vt.dot(vt) - 1
x3, x4, x5 = var('x3, x4, x5')
# moment vector
vm = Vector.zero + x3*N.i + x4*N.j + x5*N.k
momvt = vt.dot(vm)
# centers and radii
ctr1 = (0, 0, 0)
ctr2 = (1, 0, 0)
ctr3 = (0.5, sqrt(3.0)/2, 0)
ctr4 = (0.5, sqrt(3.0)/6, sqrt(6.0)/3)
radius = 0.5
centers = [ctr1, ctr2, ctr3, ctr4]
radii = [radius for _ in range(4)]
eqs = [normt, momvt]
for (ctr, rad) in zip(centers, radii):
    vc = Vector.zero + ctr[0]*N.i + ctr[1]*N.j + ctr[2]*N.k
    left = vm - vc.cross(vt)
    equ = left.dot(left) - rad**2
    eqs.append(equ)
fourspheres = []
print('the polynomial system :')
for pol in eqs:
    print(pol.expand())
    fourspheres.append(str(pol.expand()) + ';')
print('calling the blackbox solver ..')
sols = solve(fourspheres)
print('the solutions :')
for (idx, sol) in enumerate(sols):
    print('Solution', idx+1, ':')
    print(sol)
