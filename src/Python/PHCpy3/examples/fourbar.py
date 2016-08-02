"""
This scripts sets up the equations to design a 4-bar mechanism,
using python 3 and sympy.
The system appears in a paper by A.P. Morgan and C.W. Wampler on
Solving a Planar Four-Bar Design Using Continuation, published in
the Journal of Mechanical Design, volume 112, pages 544-550, 1990.
"""
from sympy import var
from sympy.matrices import Matrix

def polynomials(d0,d1,d2,d3,d4,a):
    """
    Given in d0, d1, d2, d3, d4 are the coordinates of
    the precision points, given as Matrix objects.
    Also the coordinates of the pivot in a is stored in a Matrix.
    Returns the system of polynomials to design the 4-bar
    mechanism with a coupler passing through the precision points.
    """
    # the four rotation matrices
    c1, s1, c2, s2 = var('c1, s1, c2, s2')
    c3, s3, c4, s4 = var('c3, s3, c4, s4')
    R1 = Matrix([[c1, -s1], [s1, c1]])
    R2 = Matrix([[c2, -s2], [s2, c2]])
    R3 = Matrix([[c3, -s3], [s3, c3]])
    R4 = Matrix([[c4, -s4], [s4, c4]])
    # the first four equations reflecting cos^2(t) + sin^(t) = 1
    p1, p2 = 'c1^2 + s1^2 - 1;', 'c2^2 + s2^2 - 1;'
    p3, p4 = 'c3^2 + s3^2 - 1;', 'c4^2 + s4^2 - 1;'
    # the second four equations on X
    x1, x2 = var('x1, x2')
    X = Matrix([[x1], [x2]])
    c1x = 0.5*(d1.transpose()*d1 - d0.transpose()*d0)
    c2x = 0.5*(d2.transpose()*d2 - d0.transpose()*d0)
    c3x = 0.5*(d3.transpose()*d3 - d0.transpose()*d0)
    c4x = 0.5*(d3.transpose()*d4 - d0.transpose()*d0)
    e1x = (d1.transpose()*R1 - d0.transpose())*X + c1x
    e2x = (d2.transpose()*R2 - d0.transpose())*X + c2x
    e3x = (d3.transpose()*R3 - d0.transpose())*X + c3x
    e4x = (d4.transpose()*R4 - d0.transpose())*X + c4x
    s1, s2 = str(e1x[0]) + ';', str(e2x[0]) + ';'
    s3, s4 = str(e3x[0]) + ';', str(e4x[0]) + ';'
    # the third group of equations on Y
    y1, y2 = var('y1, y2')
    Y = Matrix([[y1], [y2]])
    c1y = c1x - a.transpose()*(d1 - d0)
    c2y = c2x - a.transpose()*(d2 - d0)
    c3y = c3x - a.transpose()*(d3 - d0)
    c4y = c4x - a.transpose()*(d4 - d0)
    e1y = ((d1.transpose() - a.transpose())*R1 \
         - (d0.transpose() - a.transpose()))*Y + c1y
    e2y = ((d2.transpose() - a.transpose())*R2 \
         - (d0.transpose() - a.transpose()))*Y + c2y
    e3y = ((d3.transpose() - a.transpose())*R3 \
         - (d0.transpose() - a.transpose()))*Y + c3y
    e4y = ((d4.transpose() - a.transpose())*R4 \
         - (d0.transpose() - a.transpose()))*Y + c4y
    s5, s6 = str(e1y[0]) + ';', str(e2y[0]) + ';'
    s7, s8 = str(e3y[0]) + ';', str(e4y[0]) + ';'
    return [p1, p2, p3, p4, s1, s2, s3, s4, s5, s6, s7, s8]

def main():
    """
    Defines the pivot and generates random coordinates,
    uniformly distributed in [-1, +1], for the five points
    through which the coupler must pass.
    Calls the blackbox solver to solve the system.
    """
    # the five precision points are random points
    from random import uniform as u
    pt0 = Matrix(2, 1, lambda i,j: u(-1,+1))
    pt1 = Matrix(2, 1, lambda i,j: u(-1,+1))
    pt2 = Matrix(2, 1, lambda i,j: u(-1,+1))
    pt3 = Matrix(2, 1, lambda i,j: u(-1,+1))
    pt4 = Matrix(2, 1, lambda i,j: u(-1,+1))
    # the pivot is a
    piv = Matrix([[1], [0]])
    equ = polynomials(pt0,pt1,pt2,pt3,pt4,piv)
    print('the polynomial system :')
    for pol in equ:
        print(pol)
    from phcpy.solver import solve
    sols = solve(equ)
    print('the solutions :')
    for sol in sols:
        print(sol)
    print('computed', len(sols), 'solutions')

if __name__ == "__main__":
    main()
