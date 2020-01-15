"""
The alternating group A4 is isomorphic to the isometries of a regular tetrahedron. This module allows one to calculate the characters of that representation.
"""
# Except it doesn't because the tetrahedron needs to be centered at the origin and you've put it with *its face* centered at the origin, so FIX THAT!

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sympy import symbols, simplify
from sympy import pi, flatten, sqrt
from sympy import sin as _sin
from sympy import cos as _cos
sin = lambda x: _sin(pi*x/180)
cos = lambda x: _cos(pi*x/180)

#sin = lambda x: _sin(x)
#cos = lambda x: _cos(x)

from sympy import Matrix, Rational
array = Matrix
norm = lambda v: v.norm()
cross = lambda x,y: x.cross(y)
dot = lambda x,y: x.dot(y)

class Triangle:
    """
    Equilateral triangle.
    """

    def __init__(self, L=1, c=(0,0), theta=0):
        """
        Create triangle.

        @param L: Side length.
        @param c: Triangle's center coordinates.
        @param theta: Rotation about triangle's center, in degrees.
        """

        self.L = L
        self.c = c
        self.theta = theta

    def contains(self, point):
        """
        Check if point is inside triangle using the ray casting algorithm.

        @param point: x and y coordinates.
        """

        ymin, ymax = map(lambda f: f(map(lambda v: v[1], self.vertices.values())), (min, max))

        x0, y0 = point[0], point[1]
        cuts = 0

        for _, eqn in self.equations.items():
            m, n = eqn['m'], eqn['n']

            if m==0:
                continue
            else:
                if ymin<=y0<=ymax and (y0-n)/m>=x0:
                    cuts += 1

        return bool(cuts%2)

    def plot(self, nPoints=100, ax=None, **kwds):
        """
        Plot triangle.
        """

        if ax is None:
            _, ax = plt.subplots()

        plt.sca(ax)

        x = np.linspace(self.limits['xmin'], self.limits['xmax'], nPoints)
        xl = x[int(nPoints/2):]
        xr = x[:int(nPoints/2)]

        ya = self.equations['a']['m']*x+self.equations['a']['n']
        yb = self.equations['b']['m']*xl+self.equations['b']['n']
        yc = self.equations['c']['m']*xr+self.equations['c']['n']

        plt.plot(x, ya, **kwds)
        plt.plot(xl, yb, **kwds)
        plt.plot(xr, yc, **kwds)

        return ax

    @property
    def height(self):
        """
        Triangle's height.
        """
        return sqrt(self.L**2 - (self.L*Rational(1,2))**2)

    @property
    def outer_radius(self):
        """
        Distance from vertex to center.
        """
        return (self.L*Rational(1,2))* (1/cos(30))

    @property
    def inner_radius(self):
        """
        Distance from side to center.
        """
        return sqrt(self.outer_radius**2 - (self.L*Rational(1,2))**2)

    @property
    def vertices(self):
        """
        Triangle vertices in this order:

                    A
                   / \
                  B---C
        """

        R = self.outer_radius

        A = self.c + R*array([cos(self.theta+90), sin(self.theta+90)])
        B = self.c + R*array([cos(self.theta+210), sin(self.theta+210)])
        C = self.c + R*array([cos(self.theta+330), sin(self.theta+330)])

        return {'A': A, 'B': B, 'C': C}

    @property
    def equations(self):
        """
        Slope and intercept defining triangle's sides in the following order:

                    / \
                   c   b
                  /__a__\
        """

        A, B, C = map(lambda k: self.vertices[k], ('A', 'B', 'C'))

        a = {'m': (C[1]-B[1])*(1/(C[0]-B[0]))}
        a['n'] = B[1]-a['m']*B[0]

        b = {'m': (C[1]-A[1])*(1/(C[0]-A[0]))}
        b['n'] = A[1]-b['m']*A[0]

        c = {'m': (B[1]-A[1])*(1/(B[0]-A[0]))}
        c['n'] = A[1]-c['m']*A[0]

        return {'a': a, 'b': b, 'c': c}

    @property
    def limits(self):
        """
        Maximum and minimum x and y coordinates.
        """
        A, B, C = map(lambda k: self.vertices[k], ('A', 'B', 'C'))

        xmin = min(A[0], B[0], C[0])
        xmax = max(A[0], B[0], C[0])
        ymin = min(A[1], B[1], C[1])
        ymax = max(A[1], B[1], C[1])

        return {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax}

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, new):
        if new>0:
            self._L = new
        else:
            raise ValueError('Triangle\'s side length must be positive!')

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, new):
        if not hasattr(new, '__len__'):
            raise TypeError('Center must be an iterable containing two coordinates!')
        if len(new)==2:
            self._c = array(new)
        else:
            raise ValueError('Center must be an iterable containing two coordinates!')

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, new):
        self._theta = new%360

class Tetrahedron:
    """
    Regular tetrahedron.
    """

    def __init__(self, face=None):
        """
        Create regular tetrahedron from face.

        @param base: Triangle object.
        """

        self.face = face

    @property
    def height(self):
        """
        Tetrahedron's height.
        """
        return sqrt(self.face.height**2 - self.face.inner_radius**2)

    @property
    def vertices(self):
        """
        Tetrahedron's vertices.
        """

        self._vertices = {k: array(flatten(v)+[0]) for k,v in self.face.vertices.items()}
        self._vertices['D'] = array([0, 0, self.height])

        return self._vertices

    @property
    def edges(self):
        """
        Unit vectors defining edges.
        """

        self._edges = {}

        for k, v in self.vertices.items():
            for kk, vv in self.vertices.items():

                if k==kk:
                    continue

                kkk = ''.join(sorted(k+kk))

                if kkk not in self._edges:
                    self._edges[kkk] = (v-vv)*(1/norm(v-vv))

        return self._edges

    @property
    def normals(self):
        """
        Unit normals to the faces.
        """

        self._normals = {'ABC': array([0,0,1])}

        for k, v in self.edges.items():
            for kk, vv in self.edges.items():

                if k==kk:
                    continue

                if any([v[-1]==0, vv[-1]==0]):
                    continue

                kkk = ''.join(sorted(set(k+kk)))

                if kkk not in self._normals:
                    self._normals[kkk] = cross(v,vv)*(1/norm(cross(v,vv)))

        return self._normals

    @property
    def face(self):
        return self._face

    @face.setter
    def face(self, new):
        if new is None:
            self._face = Triangle()
        else:
            self._face = new

    def __repr__(self):
        return '\n'.join(['\t{}: {}'.format(k, val.T.tolist()[0]) for k, val in self.vertices.items()])

class Rotation:
    """
    Rotation matrix.
    """

    def __init__(self, axis, angle):
        """
        Create rotation matrix.

        @param axis: 3D unit length array with axis of rotation.
        @param angle: Rotation angle, in degrees.
        """

        self.axis = axis
        self.angle = angle

    def __call__(self, v):
        """
        Rotate vector v.
        """

        new = array([0,0,0])
        new += cos(self.angle)*v
        new += sin(self.angle)*cross(self.axis,v)
        new += (1-cos(self.angle))*dot(self.axis,v)*self.axis

        return simplify(new)

    def as_matrix(self):
        """
        Get matrix elements.
        """

        x, y, z = symbols('v_1 v_2 v_3')
        v = array([x, y, z])

        w = self(v)

        m = []

        for coord in w:
            row = []
            for v in (x,y,z):
                coeff = coord.as_poly(v).coeffs()[0]
                if any([v in coeff.free_symbols for v in (x,y,z)]):
                    coeff = 0

                row += [coeff]

            m += [row]

        m = Matrix(m)

        return m

    @property
    def axis(self):
        return self._axis

    @axis.setter
    def axis(self, new):
        self._axis = new.normalized()

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, new):
        self._angle = new #%360

def character(d):
    """
    Characters.
    """

    return {k: val.trace() for k,val in d.items()}

if __name__ == '__main__':

    te = Tetrahedron()

    r1 = {face+'1': Rotation(axis, 120).as_matrix() for face, axis in te.normals.items()}
    r2 = {face+'2': Rotation(axis, 240).as_matrix() for face, axis in te.normals.items()}

    rotations = {**r1, **r2}

    crot = character(rotations)
