"""SU(N) stuff."""

import numpy as np

from sympy import simplify, symbols, ask, Q
from sympy import zeros, eye, sqrt, exp, sin, cos, atan2, re, im, diag, I, Poly, Matrix, flatten, Rational, solve, linsolve
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

class LieAlgebra:
    """Base class, for subclassing only."""

    def __init__(self):
        self._generators = None

    @property
    def generators(self):
        return self._generators

    @generators.setter
    def generators(self, new):
        self._generators = new

    @property
    def structure(self):
        """Structure constants."""

        f = np.full((self.dim, self.dim, self.dim), None)

        for i,g in enumerate(self.generators):
            for j,h in enumerate(self.generators):
                f[i,j,:] = expand_in_basis(g*h-h*g, self.generators, real=True)

        return f.astype(np.complex)

    @property
    def killing(self):
        """Killing form in coordinates given by generators."""

        fs = np.array(self.structure) #.astype(np.float)

        return np.einsum('kil,ljk->ij', fs, fs)

    def adj(self, m):
        """Find adjoint representation of matrix m."""

        comm = lambda x,y: x*y-y*x
        cs = symbols('l0:{}'.format(self.dim), real=True)

        x = zeros(self.N)
        for c, e in zip(cs, self.generators):
            x += c*e

        y = comm(m,x)
        print(y)
        coords = expand_in_basis(y, self.generators, real=True)

        return Matrix([[coords[i].subs({l: 1 if l==c else 0 for l in cs}) for c in cs] for i, _ in enumerate(cs)])

    @property
    def roots(self):
        """Root system."""

        hs = [self.adj(g) for g in self.generators if g.is_diagonal()]

        return [Matrix(v) for v in filter(any, Matrix([[h.tolist()[i][i] for i in range(self.dim)] for h in hs]).T.tolist())]

    @property
    def _roots(self):
        """Get root system."""

        # TODO I'm using a dirty trick - computing commutators of generators
        # What I should do in general is get the adjoint representation
        # of the elements in the cartan subalgebra and diagonalise

        root = []
        gs = self.generators
        hs = [g for g in gs if g.is_diagonal()]
        xs = [g for g in gs if not g.is_diagonal()]

        for i, x in enumerate(gs):

            if x in hs:
                continue

            r = []

            for h in hs:
                comm = expand_in_basis(h*x-x*h, gs, real=False)
                print(comm)
                rr = [c for j,c in enumerate(comm) if c!=0 or j==i]
                if len(rr)==1:
                    r += [rr[0]]
                else:
                    break

            root += [Matrix(r)] if any(r) else []

        return root

class SL(LieAlgebra):
    """Special linear Lie algebra."""

    def __init__(self, N):
        """SL(N)"""

        super(SL, self).__init__()

        self.N = N
        self.dim = self.N**2-1

        self._generators = self._getDefaultGenerators()

    def _getDefaultGenerators(self):
        """Standard basis."""

        gens = [diag(*[int((-1)**(i-j)*int(i==j or i==j+1)) for i in range(self.N)]) for j in range(self.N-1)]

        for i in range(self.N):
            for j in range(self.N):
                if i!=j:
                    eij = Matrix([[int(k==i and l==j) for k in range(self.N)] for l in range(self.N)])
                    gens += [eij]

        return gens

    def __repr__(self):
        return 'sl({})'.format(self.N)

class SU(LieAlgebra):
    """Special unitary Lie algebra."""

    def __init__(self, N):
        """SU(N)"""

        super(SU, self).__init__()

        self.N = N
        self.dim = self.N**2-1

        self._generators = self._getDefaultGenerators()

    def _getDefaultGenerators(self):
        """Standard basis."""

        gens = [I*diag(*[int((-1)**(i-j)*int(i==j or i==j+1)) for i in range(self.N)]) for j in range(self.N-1)]

        for i in range(self.N):
            for j in range(i,self.N):
                if i!=j:
                    eij = Matrix([[int(k==i and l==j) for k in range(self.N)] for l in range(self.N)])
                    gens += [eij.T-eij]
                    gens += [I*(eij+eij.T)]

        return gens

    def __repr__(self):
        return 'su({})'.format(self.N)

def GellMann(anti=True):
    """The (anti-)Hermitian Gell-Mann matrices."""

    l1 = Matrix([
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 0]
        ])

    l2 = Matrix([
        [0, -I, 0],
        [I, 0, 0],
        [0, 0, 0]
        ])

    l3 = Matrix([
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, 0]
        ])

    l4 = Matrix([
        [0, 0, 1],
        [0, 0, 0],
        [1, 0, 0]
        ])

    l5 = Matrix([
        [0, 0, -I],
        [0, 0, 0],
        [I, 0, 0]
        ])

    l6 = Matrix([
        [0, 0, 0],
        [0, 0, 1],
        [0, 1, 0]
        ])

    l7 = Matrix([
        [0, 0, 0],
        [0, 0, -I],
        [0, I, 0]
        ])

    l8 = Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, -2]
        ])

    ls = [l1, l2, l3, l4, l5, l6, l7, sqrt(3)*Rational(1,3)*l8]

    ls = [sqrt(12)*Rational(1,12)*l for l in ls]

    if anti:
        ls = [I*l for l in ls]

    return ls
