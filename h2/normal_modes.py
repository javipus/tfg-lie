"""
Normal modes of three identical masses connected by identical springs on a frictionless table.

Following Georgi's Lie Algebras in Particle Physics 1.16-17 the system looks like this:

    2 - 1
     \ /
      3

and we use coordinates

    (r11, r12, r21, r22, r31, r32) in R6

where the first index labels the mass and the second index labels the physical coordinate - either x or y.
"""

from functools import reduce

from sympy.physics.quantum import TensorProduct
from sympy import latex

from permutations import S3

def eigen(m, return_multiplicity=False):
    if return_multiplicity:
        return {lmbd: (mul, basis) for lmbd, mul, basis in m.eigenvects()}
    else:
        return {lmbd: basis for lmbd, _, basis in m.eigenvects()}

# Symmetric group S3
s3 = S3()
N = len(s3._elements) # group order

# Representation as dihedral group D3
d3 = s3.rho_D3()

# Natural representation as permutation of basis elements of R3
r3 = s3.rho_natural()

# Regular 6-dimensional representation, with R6 ~= R3 \otimes R2
reg = [TensorProduct(x,y) for x,y in zip(r3,d3)]

# Projector onto trivial representation
p0 = reduce(lambda x,y: x+y, reg)/N

if p0.rank()==1:
    e0 = eigen(p0)[1][0].normalized()
else:
    e0 = eigen(p0)[1]

with open('p0.tex', 'w') as f:
    f.write(latex(p0))

with open('e0.tex', 'w') as f:
    f.write(latex(e0.T))

# Projector onto sign representation
p1 = reduce(lambda x,y: x+y, [s.sign*r for s,r in zip(s3._elements, reg)])/N

if p1.rank()==1:
    e1 = eigen(p1)[1][0].normalized()
else:
    e1 = eigen(p1)[1]

with open('p1.tex', 'w') as f:
    f.write(latex(p1))

with open('e1.tex', 'w') as f:
    f.write(latex(e1.T))

# Projector onto D3 representation
p2 = 2*reduce(lambda x,y: x+y, [s.trace()*r for s,r in zip(d3,reg)])/N

if p2.rank()==1:
    e2 = eigen(p2)[1][0].normalized()
else:
    e2 = eigen(p2)[1]

with open('p2.tex', 'w') as f:
    f.write(latex(p2))

#with open('e2.tex', 'w') as f:
#    f.write(latex([]))
