"""SU(2) stuff."""

from sympy import simplify, symbols, ask, Q
from sympy import eye, sqrt, exp, sin, cos, atan2, re, im, I, Matrix, Rational
from sympy.algebras.quaternion import Quaternion
from sympy.physics.matrices import msigma

def f(a, b, c, d):
    """Isomorphism from unitary quaternions to SU(2)."""

    _atan2 = lambda y, x: atan2(y, x) if not x==y==0 else 0

    theta = _atan2(b, a)
    phi = _atan2(d, c)
    cosu = sqrt(a**2+b**2)
    sinu = sqrt(1-cosu**2)

    return Matrix([
        [exp(I*theta)*cosu, exp(I*phi)*sinu],
        [-exp(-I*phi)*sinu, exp(-I*theta)*cosu]
        ])

def finv(u):
    """Inverse isomorphism from SU(2) to unitary quaternions."""

    a, b = re(u[0]), im(u[0])
    c, d = re(u[1]), im(u[1])

    return a, b, c, d

# Sanity checks
assert finv(f(1,0,0,0)) == (1, 0, 0, 0)
assert finv(f(0,1,0,0)) == (0, 1, 0, 0)
assert finv(f(0,0,1,0)) == (0, 0, 1, 0)
assert finv(f(0,0,0,1)) == (0, 0, 0, 1)

# su(2)
Id, sx, sy, sz = eye(2), msigma(1), msigma(2), msigma(3)

assert f(1,0,0,0) == Id
assert f(0,1,0,0) == I*sz
assert f(0,0,1,0) == I*sy
assert f(0,0,0,1) == I*sx

# Division
def div(q1, q2):
    """Divide two quaternions as q1/q2 = f^{-1}[f(q1)*f(q2)^{-1}]."""

    return finv(f(*q1)*(f(*q2).inv()))

q1 = Rational(1,5), Rational(2,5), Rational(2,5), Rational(4,5)
q2 = Rational(2,3), Rational(2,3), Rational(1,3), 0
qd = simplify(div(q1,q2))

# Killing form
def killing(x, y):
    return -Rational(1,2)*(x*y).trace()

# Generic matrices in su(2)
a1, b1, c1 = symbols('a1 b1 c1', real=True)
a2, b2, c2 = symbols('a2 b2 c2', real=True)

x = I*(a1*sx + b1*sy + c1*sz)
y = I*(a2*sx + b2*sy + c2*sz)

# Killing form is an inner product - trace is linear and symmetric
assert ask(Q.zero(simplify(killing(x,y)-simplify(killing(x,y))))) # symmetric
assert ask(Q.nonnegative(simplify(killing(x,x)))) # positive definite - sum of squares => 0 iff all are 0

# Orthongonality of Pauli matrices
assert ask(Q.zero(killing(sx,sy)))
assert ask(Q.zero(killing(sx,sz)))
assert ask(Q.zero(killing(sy,sz)))

# Structure constants
e1, e2, e3 = I*sx, I*sy, I*sz
comm = lambda x, y: x*y-y*x
f = [[comm(u, v) for u in (e1, e2, e3)] for v in (e1, e2, e3)]
