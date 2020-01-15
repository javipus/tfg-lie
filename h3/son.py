"""SO(N) stuff."""

import numpy as np

from sympy import simplify, symbols, ask, Q
from sympy import eye, sqrt, exp, sin, cos, atan2, re, im, I, Matrix, flatten, Rational, solve, linsolve
#from sympy.algebras.quaternion import Quaternion
#from sympy.physics.matrices import msigma

def killing(x, y):
    """Killing inner product in Lie algebras."""

    return -Rational(1,2)*(x*y).trace()

def expand_in_basis(m, basis, real=True):
    """Expand matrix in basis."""

    cs = symbols('c0:{}'.format(len(basis)), real=real)
    eqs = flatten(m-sum(map(lambda c, e: c*e, cs, basis), np.full(m.shape, 0)))
    sol = linsolve(eqs, cs)

    if len(sol.args)==1:
        return sol.args[0]
    else:
        raise Exception('No valid solution: {}'.format(sol))

class SO:
    """Lie algebra so(N)."""

    def __init__(self, N):
        """Create algebra given natural N."""

        self.N = N
        self.dim = int(N*(N-1)/2)

    @property
    def generators(self):
        """Standard basis."""

        gens = []

        for i in range(1,self.N+1):
            for j in range(1,i):
                m = Matrix([[int(k==i and l==j) for k in range(1,self.N+1)] for l in range(1,self.N+1)])
                m = (-1)**(i+j)*(m-m.T)
                gens += [m]

        return gens

    @property
    def structure(self):
        """Structure constants."""

        f = np.full((self.dim, self.dim, self.dim), None)

        for i,g in enumerate(self.generators):
            for j,h in enumerate(self.generators):
                f[i,j,:] = expand_in_basis(g*h-h*g, self.generators, real=True)

        return f.tolist()

    def __repr__(self):
        return 'so({})'.format(self.N)


# DEPRECATED #
# FAILED #
# MISERY #
# HORROR #

delta = lambda i,j: Q.zero(i-j)

class eij(dict):

    def __init__(self, i, j):
        self[(i,j)] = 1

    def __mul__(self, other):
        new = eij()

        for (i,j),u in self.items():
            for (k,l),v in other.items():
                new[(i,l)] = delta(j,k)

        return new

    def __add__(self, other):
        return {k: self.get(k,0)+other.get(k,0) for k in set(list(self.keys())+list(other.keys()))}

    def __neg__(self):
        return {k: -v for k, v in self.items()}

    def __sub__(self, other):
        return self.__add__(other.__neg__())

Eij = lambda i,j: eij(i,j)-eij(j,i)
comm = lambda u,v: u*v-v*u
