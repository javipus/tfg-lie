"""
Irreducible representations of su(2).
"""

from sympy import pi, sqrt, exp, diag, I, Matrix, Rational

class Irrep:

    def __init__(self, n):
        """Create irreducible representation of dimension dim."""

        self.n = n
        self._dim = self.n+1
        self._weights = range(-self.n, self.n+1, 2)

        self.z = diag(*[m for m in self._weights])

        self._raise = Rational(1,2)*Matrix([[sqrt(self.n*(self._dim)-i*(i+1)) if i==j+1 else 0 for j in range(self._dim)] for i in range(self._dim)])
        self._lower = Rational(1,2)*Matrix([[sqrt(self.n*(self._dim)-i*(i-1)) if i==j-1 else 0 for j in range(self._dim)] for i in range(self._dim)])

        self.x = Rational(1,2) * (self._raise + self._lower)
        self.y = -I*Rational(1,2) * (self._raise - self._lower)

    def __repr__(self):
        return 'pi_{}'.format(self.n)
