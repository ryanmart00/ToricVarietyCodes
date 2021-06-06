from fields import *
from itertools import product
from math import log

"""
    Generates the set of all vectors of length n over the field F
    as an iterable
"""
def vectors(F, n):
    Fn = [F.set() for i in range(n)]
    return product(*Fn)

def nonzeroVectors(F, n):
    Fn = [F.nonzero() for i in range(n)]
    return product(*Fn)

def HammingDist(x,y):
    return len([1 for i in range(len(x)) if x[i] != y[i]])

def HammingWeight(x):
    return len([1 for i in range(len(x)) if x[i] != 0])

class Vector(tuple):
    def __add__(self, other):
        if other == 0:
            return self
        if len(self) != len(other):
            raise TypeError("Mismatched Vector length")
        return Vector([self[i] + other[i] for i in range(len(self))])

    def __radd__(self, other):
        return self + other

    def __mul__(self, scalar):
        return Vector([scalar * self[i] for i in range(len(self))])

    def __rmul__(self, scalar):
        return self * scalar

    def __neg__(self):
        return Vector([-x for x in self])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

class LinearCode:
    def __init__(self, F, basis):
        if len(basis) == 0:
            raise Exception("Cannot construct the empty code")
        self.n = len(basis[0])
        if not all(len(el) == self.n for el in basis):
            raise Exception("Cannot construct a code with mismatched dimension")
        if not all(all(isinstance(x, F) for x in el) for el in basis):
            raise Exception("Cannot construct a code with mismatched field")
        self.field = F
        self.basis = basis

    def k(self):
        return len(self.basis)

    def d(self):
        weight = self.n
        for x in Span(self.field, self.basis):
            w = HammingWeight(x)
            if w < weight and w != 0:
                weight = w
        return weight



    def minDist(self):
        return min(HammingDist(x,y) for x in self.elements for y in self.elements if x != y)

def Span(F, basis):
    vecs = []
    for x in vectors(F,len(basis)):
        y = sum(x[i] * basis[i] for i in range(len(basis)))
        if not y in vecs:
            vecs = vecs + [y]
    return vecs

def CodeFromPolynomialBasis(F, basis, n):
    codeBasis = [Vector([f(x) for x in nonzeroVectors(F,n)]) for f in basis]
    return LinearCode(F, codeBasis)

def CodeFromLatticePoints(F, points):
    if len(points) == 0:
        raise Exception("No points given...")
    n = len(points[0])
    if not all(len(p) == n for p in points):
        raise Exception("Not all points were of the same dimension")
    P = [PolynomialsOver(F, n)]
    for i in range(n-1):
        P = [P[0].R] + P
    basis = []
    for point in points:
        monomial = F(1)
        for i in range(n):
            monomial = P[i](point[i] * [F(0)] + [monomial])
        basis = basis + [monomial]
    return CodeFromPolynomialBasis(F, basis, n)
