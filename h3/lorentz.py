"""
Finite-dimensional (m,n)-representations of SO(1,3) - connected component of the Lorentz group containing the identity.
"""

from sympy import Rational, Matrix, eye, sqrt, exp, symbols, simplify
from sympy.physics.quantum import TensorProduct as TP
from sympy.printing import pretty

class Pi:

    def __init__(self, m, n, theta=None, beta=None):
        """(m,n)-irrep of dimension d=(2m+1)(2n+1)."""

        self.m = m
        self.n = n

        self.theta = theta or symbols('theta', real=True)
        self.beta = beta or symbols('beta', real=True)

        self.dim = (2*self.m+1)*(2*self.n+1)

        d1, d2 = int(2*self.m+1), int(2*self.n+1)

        # Infinitesimal generators - Lie algebra
        self._gen_j = {'J{}'.format(i): TP(eye(d1), self._j(n, i)) + TP(self._j(m, i), eye(d2)) for i in (1,2,3)}
        self._gen_k = {'K{}'.format(i): 1j*(TP(eye(d1), self._j(n, i)) - TP(self._j(m, i), eye(d2))) for i in (1,2,3)}
        self.generators = {**self._gen_j, **self._gen_k}

        # Transformations - Lie group
        self._g_j = {ji: simplify(exp(self.theta * j)) for ji, j in self._gen_j.items()}
        self._g_k = {ki: simplify(exp(self.beta * k)) for ki, k in self._gen_k.items()}
        self._g = {**self._g_j, **self._g_k}

    def _j(self, n, i):
        """Auxiliary spin matrices aka ladder operators."""

        rng = range(-n,n+1) if int(n)==n else [i + Rational(1/2) for i in range(int(-n-Rational(1/2)), int(n+Rational(1/2)))]

        upper = lambda n: [[sqrt((n-j)*(n+j+1)) if i==j+1 else 0 for i in rng] for j in rng] or [0]
        lower = lambda n: [[sqrt((n+j)*(n-j+1)) if i==j-1 else 0 for i in rng] for j in rng] or [0]
        diag = lambda n: [[i if i==j else 0 for i in rng] for j in rng] or [0]

        if i==1:
            J = Rational(1,2) * (Matrix(upper(n)) + Matrix(lower(n)))
        elif i==2:
            J = -1j*Rational(1,2) * (Matrix(upper(n)) - Matrix(lower(n)))
        elif i==3:
            J = Matrix(diag(n))
        else:
            raise ValueError('Index i must be <=3!')

        return J

    def __add__(self, other):
        """Direct sum."""
        # TODO: need different class to represent reducible representations

        return Pi

    def __repr__(self):
        return 'Rep({}, {})'.format(self.m, self.n)

class DirectSum:
    """Direct sum of several irreps."""

    def __init__(self, *args):
        """Create direct sum."""

        self.

if __name__=='__main__':

    # Scalar field - e.g. Higgs
    sc = Pi(0, 0)

    # Weyl fermions s=1/2 - e.g. left and right handed electron/positron
    lweyl = Pi(Rational(1,2), 0)
    rweyl = Pi(0, Rational(1,2))

    # Gauge bosons - e.g. photon
    gauge = Pi(Rational(1,2), Rational(1,2))

    # Traceless symmetric tensor - e.g. stress-energy T^mu_nu aka graviton
    grav = Pi(1,1)
