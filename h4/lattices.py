"""
Lattices, root systems, Dynkin diagrams or whatever.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sympy import sqrt, Matrix, Rational

# How to get the Dynkin diagram associated with a Lie group?
#
# 1. Get the generators of the Lie algebra
# 2. Calculate the Cartan subalgebra, i.e. the tangent space of the maximal torus in G - that is, the matrices must be diagonal
# 3. Find al the nonzero a_ij satisfying [H_i, X_j] = a_ij*X_j with H_i generating the Cartan subalgebra and X_j generating the whole algebra
# 4. Now (a_ij) is a root system, find its simple roots somehow - they are a basis for the whole space
# 5. Use the rules to generate a Dynkin diagram from the simple roots

def _inv(m):
    try:
        return m.inv()
    except AttributeError:
        return np.linalg.inv(m)

def _array(vects):
    if isinstance(vects[0], np.ndarray):
        return np.array(vects)
    else:
        return Matrix([v.T for v in vects]).T

def hexlat(n=3, dual=False):
    """Plot hexagonal lattice in 2D with coefficients up to n."""

    if dual:
        e1 = np.array([1, np.sqrt(3)/3])
        e2 = np.array([0, 2*np.sqrt(3)/3])
    else:
        e1 = np.array([1, 0])
        e2 = np.array([-.5, np.sqrt(3)/2])

    points = np.array([k1*e1+k2*e2 for k1 in range(-n, n+1) for k2 in range(-n, n+1)])

    return points.T

def dual(*lattice_vectors):
    return _inv(_array(lattice_vectors)).T


class Lattice:

    def __init__(self, dots):
        """
        Create lattice from set of pairwise dot products.

        @param dots: Dictionary with pairs (ei, ej) as keys and inner products <ei,ej> as values.
        """

        self._dots = dots
        self._names = [k[0] for k, _ in self._dots.items()]

        # Put inner products in d-by-d dataframe
        self.inners = pd.DataFrame(0, index=self._names, columns=self._names)

        for ei in self._names:
            for ej in self._names:
                if (ei,ej) in self._dots:
                    self.inners.loc[ei,ej] = self._dots[(ei,ej)]
                else:
                    self.inners.loc[ei,ej] = self._dots[(ej,ei)]

        # Put angles in d-by-d dataframe
        self.angles = pd.DataFrame(0, index=self._names, columns=self._names)

        for ei in self._names:
            for ej in self._names:
                self.angles.loc[ei,ej] = self.inners.loc[ei,ej]/(self.inners.loc[ei,ei]*self.inners.loc[ej,ej])**(.5)

        # Diagonal contains squared norms
        self.norms = dict(zip(self._names, self._data.values.diag()**(.5)))

        # Find vectors

    def draw(self, n=3):
        """Draw lattice."""


