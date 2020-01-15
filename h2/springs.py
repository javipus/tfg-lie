"""
Drawing springs and masses in 2D.
"""

import numpy as np
import matplotlib.pyplot as plt

class Spring:
    """Curly spring."""

    def __init__(self, x0, y0, theta, l, slant, r, w, nPoints=1000):
        """
        Create spring.

        @params x0, y0: Center coordinates.
        @params theta: Angle wrt OY, in radians.
        @params l: Length.
        @params r: Radius.
        @params w: Frequency.
        """

        self.x0, self.y0 = x0, y0
        self.theta = theta
        self.l = l
        self.slant = slant
        self.r = r
        self.w = w
        self.nPoints = nPoints

        R = np.array([
            [np.cos(self.theta), -np.sin(self.theta)],
            [np.sin(self.theta), np.cos(self.theta)]
            ])

        self.t = np.linspace(-1-np.pi/self.w, 1, self.nPoints)

        self.x = self.r * np.cos(self.w*self.t+np.pi/2)
        self.y = .5 * self.l * self.t + self.slant*np.sin(self.w*self.t+np.pi/2)
        self.z = self.r * np.sin(self.w*self.t+np.pi/2)

        self.x, self.y = np.dot(R, np.vstack([self.x, self.y]))

        self.x += self.x0
        self.y += self.y0

    def draw2D(self, ax=None, **kwds):
        """
        Draw spring's projection onto XY plane.
        """

        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(self.x, self.y, **kwds)

        return ax

if __name__=='__main__':

    # Spring parameters
    r = .25
    w = 4 * 2*np.pi
    slant = .15
    l = 2*np.sqrt(3)

    # Arrows
    arrNorm = 2.5
    arrWidth = .5

    mode0 = [
            arrNorm*np.array([-1/2, np.sqrt(3)/6]),
            arrNorm*np.array([0, -np.sqrt(3)/3]),
            arrNorm*np.array([1/2, np.sqrt(3)/6])
            ]

    mode1 = [
            arrNorm*np.array([-np.sqrt(3)/6, -1/2]),
            arrNorm*np.array([np.sqrt(3)/3,0]),
            arrNorm*np.array([-np.sqrt(3)/6, 1/2])
            ]

    fig, axs = plt.subplots(1,2,figsize=(20,10))

    for j, (ax, mode) in enumerate(zip(axs, (mode0, mode1))):

        # Lower
        Spring(x0=0, y0=-1, theta=-np.pi/2, l=l, slant=slant, w=w, r=r).draw2D(ax=ax, color='k', lw=2)

        # Left
        Spring(x0=-np.sqrt(3)/2, y0=1/2, theta=-np.pi/6, l=l, slant=slant, w=w, r=r).draw2D(ax=ax, color='k', lw=2)

        # Right
        Spring(x0=np.sqrt(3)/2, y0=1/2, theta=np.pi/6, l=l, slant=slant, w=w, r=r).draw2D(ax=ax, color='k', lw=2)

        # Masses
        top = plt.Circle((0,2), .5, edgecolor='k', facecolor='w', zorder=2, lw=2)
        left = plt.Circle((-np.sqrt(3),-1), .5, edgecolor='k', facecolor='w', zorder=2, lw=2)
        right = plt.Circle((np.sqrt(3),-1), .5, edgecolor='k', facecolor='w', zorder=2, lw=2)
        ax.add_artist(top)
        ax.add_artist(left)
        ax.add_artist(right)

        ax.add_artist(plt.Arrow(np.sqrt(3), -1, mode[0][0], mode[0][1], width=arrWidth, color='k'))
        ax.add_artist(plt.Arrow(0, 2, mode[1][0], mode[1][1], width=arrWidth, color='k'))
        ax.add_artist(plt.Arrow(-np.sqrt(3), -1, mode[2][0], mode[2][1], width=arrWidth, color='k'))

        # Text
        fsize = 40
        ax.text(0, 2, '2', color='k', size=fsize, horizontalalignment='center', verticalalignment='center')
        ax.text(-np.sqrt(3), -1, '3', color='k', size=fsize, horizontalalignment='center', verticalalignment='center')
        ax.text(np.sqrt(3), -1, '1', color='k', size=fsize, horizontalalignment='center', verticalalignment='center')

        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 3)

        ax.set_axis_off()

    plt.savefig('springs_nm_1.png')
    plt.close()
