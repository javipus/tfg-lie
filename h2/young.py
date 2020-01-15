"""
Young diagrams and tableaux.
"""

from functools import reduce
from math import factorial

from permutations import S

class YoungDiagram:

    def __init__(self, l):
        """
        Creates Young diagram for a given integer partition.

        @param: Iterable of non-increasing positive integers.
        """

        # Partition
        self.l = l

        # Rows and columns
        self.m, self.n = len(self.l), max(self.l)

        # Order and symmetric group associated to it
        self.N = sum(self.l)
        self._sn = S(self.N)

        # Cells
        self.cells = [(i+1,j+1) for i in range(self.m) for j in range(self.l[i])]

    def hook(self, cell):
        """
        Get hook of cell in diagram.
        """

        assert cell in self.cells, 'Cell {} not in diagram!'.format(cell)

        right = lambda c: c[0]==cell[0] and c[1]>cell[1]
        down = lambda c: c[0]>cell[0] and c[1]==cell[1]

        return set([_cell for _cell in self.cells if right(_cell) or down(_cell)] + [cell])

    def hookLength(self, cell):
        """
        Calculate hook length of given cell.
        """

        return len(self.hook(cell))

    @property
    def dim(self):
        """
        Calculate dimension of irrep associated to diagram using the hook length formula.
        """

        d = factorial(self.N) / reduce(lambda x,y: x*y, [self.hookLength(cell) for cell in self.cells])

        assert d==int(d), 'Non-integer dimension = {}, huh?'.format(d)

        return int(d)

    @property
    def tableaux(self):
        """
        Tableaux associated to this diagram.
        """

        return [YoungTableau(self, sigma(range(1,self.N+1))) for sigma in self._sn._elements]

    def subspace(self):
        """
        Subspace associated with this representation in the regular representation.
        """

        pass

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, new):
        assert all([l>0 for l in new]), 'All integers must be positive!'
        assert all([l0>=l1 for l1, l0 in zip(new[1:], new[:-1])]), 'Integer sequence must be non-increasing!'
        self._l = new

    def __repr__(self):
        return '\n'.join(['\t'+'O'*k for k in self.l])


class YoungTableau:

    def __init__(self, diagram, labels):
        """
        Create Young tableau given shape and labeling.

        @param diagram: YoungDiagram object.
        @param labels: Diagram numbering, from left to right, top to bottom.
        """

        self._diagram = diagram
        self._labels = labels

        assert len(self._diagram.cells)==len(self._labels), 'Number of cells must equal number of labels!'

        self._cells = zip(self._diagram.cells, self._labels)

    def __repr__(self):
        _labels = self._labels.copy()
        return '\n'.join(['\t'+' '.join([str(_labels.pop(0)) for i in range(row)]) for row in self._diagram.l])
