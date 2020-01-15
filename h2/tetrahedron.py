from sympy import sqrt, Rational, Matrix, diag

# A MAGICIAN DID THIS #
# Vertices
v1 = Matrix([sqrt(Rational(8,9)), 0, -Rational(1,3)])
v2 = Matrix([-sqrt(Rational(2,9)), sqrt(Rational(2,3)), -Rational(1,3)])
v3 = Matrix([-sqrt(Rational(2,9)), -sqrt(Rational(2,3)), -Rational(1,3)])
v4 = Matrix([0,0,1])

# Change of basis
c = Matrix([
    [sqrt(Rational(1,6)), sqrt(Rational(2,3)), -sqrt(Rational(1,6))],
    [sqrt(Rational(1,2)), 0, sqrt(Rational(1,2))],
    [-sqrt(Rational(1,3)), sqrt(Rational(1,3)), sqrt(Rational(1,3))]
    ])

# /END MAGICIAN

# 2-cycles
r1 = diag(1,-1,-1)
r2 = diag(-1,1,-1)
r3 = diag(-1,-1,1)

# Check that it works
print('All these should be zero:')
print(c*r1*c.inv()*v1-v2)
print(c*r1*c.inv()*v3-v4)
print(c*r2*c.inv()*v2-v3)
print(c*r2*c.inv()*v1-v4)
print(c*r3*c.inv()*v2-v4)
print(c*r3*c.inv()*v1-v3)
