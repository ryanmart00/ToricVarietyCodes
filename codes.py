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

"""
    for x of length n, ReducedHammingWeight(x) is 
    n - largest number of equal digits in x 

    The idea is, we can remove the vector [1,1,1,...,1] from our
    basis if and then the minimum distance will be the minimum of 
    this Reduced Distance. This is because the largest number of equal
    digits y in x corresponds to the Hamming Weight of x - y[1,1,1,...,1] 
    and will be the smallest Hamming Weight of that form. Finally, every vector 
    is of this form so we have it
"""
def ReducedHammingWeight(x):
    digits = {}
    maxParity = 0
    for el in x:
        if not str(el) in digits.keys():
            digits[str(el)] = 1
        else:
            digits[str(el)] = digits[str(el)] + 1
        if digits[str(el)] > maxParity:
            maxParity = digits[str(el)]
    return len(x) - maxParity

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
    """
        Set reduced to true if the element [1,1,1,1 ... , 1] is in the basis 
        but you've removed it for performance reasons
    """
    def __init__(self, F, basis, reduced = False):
        if len(basis) == 0:
            raise Exception("Cannot construct the empty code")
        self.n = len(basis[0])
        if not all(len(el) == self.n for el in basis):
            raise Exception("Cannot construct a code with mismatched dimension")
        if not all(all(isinstance(x, F) for x in el) for el in basis):
            raise Exception("Cannot construct a code with mismatched field")
        self.field = F
        self.basis = basis
        self.reduced = reduced

        
    def k(self):
        return len(self.basis)

    def d(self):
        if self.reduced:
            return self.d_reduced()

        weight = self.n
        for x in Span(self.field, self.basis):
            w = HammingWeight(x)
            if w < weight and w != 0:
                weight = w
        return weight

    def d_reduced(self):
        weight = self.n
        for x in Span(self.field, self.basis):
            w = ReducedHammingWeight(x)
            if w < weight and w != 0:
                weight = w
        return weight

    def d_poly(self):
        weight = self.n
        for x in Span(self.field, self.basis):
            w = ReducedHammingWeight(x)
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

def CodeFromPolynomialBasis(F, basis, n, reduced=False):
    codeBasis = [Vector([f(x) for x in nonzeroVectors(F,n)]) for f in basis]
    return LinearCode(F, codeBasis, reduced)

def CodeFromLatticePoints(F, points, reduced=False):
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
    return CodeFromPolynomialBasis(F, basis, n, reduced)

def PolyBasisFromLatticePoints(F, points):
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
    return (n,basis)

def MaxZeros(F, points):
    n, basis = PolyBasisFromLatticePoints(F,points) 
    polys = [F(1)]
    num = 0
    for p in Span(F, basis):
        if p == F(0):
            continue
        temp = 0
        for x in nonzeroVectors(F,n):
            if p(x) == F(0):
                temp += 1
        if temp > num:
            polys = [p]
            num = temp
        elif temp == num:
            polys = polys + [p]

    return (num,polys)


    

