"""
Python script to compute with sympy the sines and cosines
of multiples of pi/1024 with 160 decimal places.
The output is formatted as Ada assignments of constant strings.
"""
from sympy import pi, sin, cos, evalf
S = []
C = []
for k in range(1, 257):
    p = k*pi/1024
    s = sin(p).evalf(160)
    S.append(s)
    c = cos(p).evalf(160)
    C.append(c)
print('sine table :')
for k in range(len(S)):
    print('    sint%03d : constant string' % k)
    s = str(S[k])
    L = len(s)//3
    (s0, s1, s2) = (s[0:L], s[L:2*L], s[2*L:])
    print('            := \"%s\"' % s0)
    print('             & \"%s\"' % s1)
    print('             & \"%s\";' % s2)
print('cosine table :')
for k in range(len(C)):
    print('    cost%03d : constant string' % k)
    c = str(C[k])
    L = len(s)//3
    (c0, c1, c2) = (c[0:L], c[L:2*L], c[2*L:])
    print('            := \"%s\"' % c0)
    print('             & \"%s\"' % c1)
    print('             & \"%s\";' % c2)
