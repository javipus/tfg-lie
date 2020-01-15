"""
Character of some irreducible representations of small symmetric gropus, S3, S4.
"""

from copy import deepcopy
from itertools import permutations

import numpy as np
factorial = np.math.factorial
binom = lambda n,m: int(factorial(n)/(factorial(m)*factorial(n-m)))

from sympy import sin, cos, pi, Matrix
array = Matrix

class Permutation(list):

    """
    Permutation algebra.
    """

    def __init__(self, perm):
        """
        Create permutation from list [sigma(0), ..., sigma(n-1)].
        """

        if sorted(perm) != list(range(len(perm))):
            raise ValueError('{} is not a valid permutation!'.format(perm))
        else:
            super(Permutation, self).__init__(perm)

        self._n = len(self)
        self._unit = list(range(self._n))

    @property
    def sign(self):
        return (-1)**(sum([self[i]>self[j] for i in range(self._n) for j in range(i, self._n)]))

    @property
    def order(self):
        """
        Order of the permutation in the group.
        """
        o = 1
        g = self

        while not g == self._unit:
            g *= self
            o += 1

        return o

    def as_cycle(self):
        """
        Represent permutation in cycle notation.
        """
        c = [[0]]
        idx = [0]

        while len(c) < self._n:
            while self[c[-1][-1]] != c[-1][0] or len(c[-1]) == 1:
                c[-1] += [self[c[-1][-1]]]
                idx += [self[c[-1][-1]]]
            try:
                c += [[min(set(range(self._n))-set(idx))]]
            except ValueError:
                break

        c = ''.join(map(lambda x: str(tuple(map(lambda y: y+1, x))), filter(lambda x: x[0]!=x[1] if len(x)>1 else False, c)))

        return c

    def _inverse(self):
        """
        Inverse element.
        """
        return Permutation([self.index(i) for i in range(self._n)])

    def __mul__(self, other):
        """
        Permutation composition.
        """
        return Permutation([other[v] for v in self])

    def __pow__(self, n):
        """
        Permutation powers = composition + inverses.
        """
        sgn = int(n/abs(n))
        an = abs(n)

        if an>1:
            return self.__pow__(sgn)*self.__pow__(abs(n)-1)
        else:
            if an==1:
                return self
            elif an==-1:
                return self._inverse()
            elif an==0:
                return self._unit

    def __call__(self, x): # TODO type invariance
        """
        Apply permutation to iterable x.
        """
        if len(x) == self._n:
            return [x[i] for i in self]
        else:
            raise ValueError('Set cardinality is {} but permutation length is {}'.format(len(x), self._n))

    def __hash__(self):
        """
        Make class hashable.
        """
        return hash(tuple(self))

class S:
    """
    Symmetric group.
    """

    def __init__(self, n):
        """
        Create symmetric group of order n.
        """

        if n==int(n):
            self.n = n
        else:
            raise ValueError('n must be a positive integer!')

        self._elements = [Permutation(sigma) for sigma in permutations(range(n),n)]

    def conjugates(self, perm):
        """
        Calculate conjugacy class of element.

        @param perm: Permutation contained in group.
        """

        if perm not in self._elements:
            raise TypeError('Permutation needs to be contained in group!')

        conj = [perm]

        for sigma in self._elements:
            if sigma not in conj:
                for tau in self._elements:
                    if tau*sigma*tau._inverse()==perm:
                        conj += [sigma]
                        break

        return frozenset(conj)

    def conjugacy(self):
        """
        Return set of all conjugacy classes.
        """
        return set([self.conjugates(s) for s in self._elements])

    def rho_natural(self):
        """
        Natural representation of S_n over R^n.
        """
        return [array([[int(sigma[j]==i) for i in range(self.n)] for j in range(self.n)]) for sigma in self._elements]

    def rho_trivial(self):
        """
        Trivial representation of S_n over R.
        """
        return [1 for sigma in self._elements]

    def rho_sign(self):
        """
        Sign representation of S_n over R.
        """
        return [sigma.sign for sigma in self._elements]

    def __repr__(self):
        return '\n'.join(map(str,self._elements))

    def __str__(self):
        return '\n'.join(map(str,self._elements))

class A(S):
    """
    Alternating group.
    """

    def __init__(self, n):
        super(A, self).__init__(n)

        self._elements = [sigma for sigma in self._elements if sigma.sign==1]

class S3(S):
    """
    Symmetric group S3.
    """

    def __init__(self):
        super(S3, self).__init__(3)

    def rho_D3(self):
        """
        Irrep as dihedral group D3.
        """

        # Identity
        id2 = Matrix([
            [1, 0],
            [0, 1]
            ])

        # Rotation of angle 2pi/3
        r = Matrix([
            [cos(2*pi/3), -sin(2*pi/3)],
            [sin(2*pi/3), cos(2*pi/3)]
            ])

        # Reflection about OY
        s = Matrix([
            [-1, 0],
            [0, 1]
            ])

        # Rotation of angle 4pi/3
        r2 = r*r

        # Rotate 2pi/3 then reflect
        sr = s*r

        # Rotate 4pi/3 then reflect
        sr2 = s*r2

        # Translation
        _map = {
                Permutation([0,1,2]): id2, # identity
                Permutation([2,0,1]): r, # 3-cycle
                Permutation([1,2,0]): r2, # 3-cycle
                Permutation([2,1,0]): s, # reflect about OY
                Permutation([1,0,2]): sr, # rotate then reflect
                Permutation([0,2,1]): sr2, # rotate, rotate, reflect
                }

        return [_map[k] for k in self._elements]

def subfactorial(n):
    """
    Number of derangements (permutations with no fixed points) in S_n.
    """
    return int(factorial(n) * sum([(-1)**j/factorial(j) for j in range(n+1)]))

def lmn(n,m):
    """
    Number of permutations in S_n fixing m elements.
    """
    return binom(n,m)*subfactorial(n-m)
