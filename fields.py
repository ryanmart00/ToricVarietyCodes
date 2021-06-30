from itertools import product
from math import sqrt

def isPrime(p):
    if not isinstance(p, int) or p < 2:
        return False
    if p == 2:
        return True
    if p % 2 == 0:
        return False
    for i in range(3, int(sqrt(p))+1,2):
        if p % i == 0:
            return False
    return True

def primeFactors(n):
    if not isinstance(n, int) or n < 2:
        return False
    if n % 2 == 0:
        yield 2
    while n % 2 == 0:
        n = n // 2

    for i in range(3, int(sqrt(n))+1,2):
        if n % i == 0:
            yield i
        while n % i == 0:
            n = n // i
    if n > 1:
        yield n




"""
    This meta class gives correct type equality
"""
class PrimeMeta(type):

    def __repr__(cls):
        return str(cls)

    def __str__(cls):
        return "Z/%iZ" % cls.P

    def __eq__(cls, other):
        if type(cls) == type(other):
            return cls.P == other.P
        return False

"""
    Returns the class for the field Z/pZ for the given value of p. Note that the 
    return type is a class, not an object of that class. In particular, to create 
    specific elements of Z/pZ call
    R = Z(p)
    r = R(x) # then r = x mod p 
"""
def _Z(prime):
    if not isPrime(prime):
        if not isPrime(prime):
            raise TypeError("p must be prime!")

    class Zp(tuple,metaclass=PrimeMeta):
        P = prime
        N = 1

        @classmethod
        def set(cls):
            for i in range(0, prime):
                yield cls(i)

        @classmethod
        def nonzero(cls):
            for i in range(1, prime):
                yield cls(i)

        @classmethod
        def constant(cls, x):
            return cls(x)
        
        def __new__(cls, val):
            if not isinstance(val, int):
                raise TypeError("Cannot construct an element of %s from an element of %s"\
                        % (cls, type(val)))
            return tuple.__new__(cls, (val % prime,))

        @classmethod
        def id(cls):
            return cls(0)

        @property
        def val(self):
            return self[0]

        @classmethod
        def _convert(cls, other):
            if not isinstance(other, cls):
                other = cls.constant(other)
            return other

        def __add__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Zp(self.val + other.val)

        def __radd__(self, other):
            return self.__add__(other)

        def __sub__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Zp(self.val - other.val)

        def __rsub__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__sub__(self)

        def __neg__(self):
            return Zp(-self.val)

        def __mul__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented

            return Zp(self.val * other.val)

        def __rmul__(self, other):
            return self.__mul__(other)

        def __pow__(self, power):
            if not isinstance(power, int):
                raise TypeError("Cannot raise to non-integer power")
            return Zp(pow(self.val, power, prime)) #Using pow for speed

        def __div__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented

            if other.val == 0:
                raise ZeroDivisionError()
            #Since the multiplicative group is cyclic of size
            #p-1, a** (p-2) * a = a ** (p-1) = 1
            return self * (other ** (prime - 2)) 

        def __rdiv__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__div__(self)

        __floordiv__ = __div__
        __rfloordiv__ = __rdiv__

        __truediv__ = __div__
        __rtruediv__ = __rdiv__
        
        def __str__(self):
            return str(self.val) 

        def __repr__(self):
            return str(self)

        def __eq__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return self.val == other.val;

        def __req__(self, other):
            return self.__eq__(other)

        def __ne__(self, other):
            return not self.__eq__(other)

        def __rne__(self, other):
            return not self.__eq__(other)

    return Zp

class PolyMeta(type):
    def __repr__(cls):
        return str(cls)

    def __str__(cls):
        return "%s[%s]" % (cls.R, cls.s)

    def __eq__(cls, other):
        if type(cls) == type(other):
            return cls.R == other.R
        return False

def _PolyOver(Ring, symbol='x', useParen=True):
    if not isinstance(symbol, str):
        raise TypeError("The symbol %s must be a string!", str(symbol))
    
    class Rx(tuple, metaclass=PolyMeta):
        R=Ring
        Base=Ring
        s=symbol

        @classmethod
        def set(cls, deg):
            lists = [Rx.R.nonzero()]
            for i in range(deg):
                lists = [Rx.R.set()] + lists
            

        @classmethod
        def monic(cls, deg):
            lists = []
            for i in range(deg):
                lists = [Rx.R.set()] + lists
            
            for vec in product(*lists):
                l = [x for x in vec]
                l.reverse()
                yield Rx(l+ [Rx.R.constant(1)])

        @classmethod
        def constant(cls, x):
            if isinstance(x, cls.R):
                return cls([x])
            elif isinstance(x, cls.Base) or isinstance(x, int):
                return cls([cls.R.constant(x)])
            else:
                raise TypeError("Value %s has ambiguous typing" % str(x))
        
        def __new__(cls, vals):
            vals = [val if isinstance(val, cls.R) else cls.R.constant(val)\
                    for val in vals]
            #trim high order zeros
            while len(vals) > 1 and vals[-1] == 0:
                vals = vals[:-1]
            return tuple.__new__(cls, vals)

        @property
        def vals(self):
            return self

        def degree(self):
            if self == 0:
                return -1
            return len(self) - 1

        @classmethod
        def _convert(cls, other):
            if not isinstance(other, cls):
                other = cls.constant(other)
            return other

        def __add__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            p1 = list(self)+max(0,len(other) - len(self)) * [0] 
            p2 = list(other)+max(0,len(self) - len(other)) * [0] 
            return Rx([p1[i] + p2[i] for i in range(len(p1))])

        def __radd__(self, other):
            return self.__add__(other)

        def __sub__(self, other):
            return self + (-other)

        def __rsub__(self, other):
            return -self + other

        def __neg__(self):
            return Rx([-a for a in self])

        def __mul__(self, other):
            if isinstance(other, Rx.R):
                return Rx([other * s for s in self])
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Rx([sum(self[i-j]*other[j]\
                for j in range(max(0,i-len(self)+1),1+min(i,len(other)-1)))\
                for i in range(0,len(self) + len(other)-1)])

        def __rmul__(self, other):
            return self.__mul__(other)

        def pow(self, n, mod):
            if not isinstance(n, int) or n < 0:
                raise TypeError("Power must be a positive integer, not %s" % str(n))
            return Rx._pow_helper(self, n, mod)

        @classmethod
        def _pow_helper(cls, self, n, mod):
            if n == 1:
                return self % mod
            if n == 0:
                return Rx([1])
            val = Rx._pow_helper(self, n // 2, mod) 
            if n % 2 == 0:
                return (val * val) % mod
            return (val * val * self) % mod

        def _str_helper(self, i):
            if i == 0:
                return str(self[i])
            elif i == 1:
                if self[i] == Rx.R.constant(1):
                    return symbol
                if useParen:
                    return '(' + str(self[i]) + ')' + symbol
                else:
                    return str(self[i]) + symbol
            else:
                if self[i] == Rx.R.constant(1):
                    return symbol + '^' + str(i)
                if useParen:
                    return '(' + str(self[i]) + ')' + symbol + '^' + str(i)
                else:
                    return str(self[i]) + symbol + '^' + str(i)


        def __str__(self):
            if len(self) == 0 or self == 0:
                return '0'
            i = 0
            while self[i] == 0:
                i = i + 1
            #find the first non-zero entry
            s = self._str_helper(i)
            for j in range(i+1, len(self)):
                if self[j] != 0:
                    s += '+' + self._str_helper(j)
            return s

        def __repr__(self):
            return str(self)

        """
            evaluate this polynomial
        """
        def __call__(self, x):
            if not isinstance(x, Rx.R):
                x = Rx.R.constant(x)
            return sum([self[i] * (x ** i) for i in range(len(self))])


        """
            This uses the fact that polynomials over a field are a Euclidean Domain
            to find polynomials q, r such that 

            self = q * other + r

            where r.degree < other.degree

            We use long division
            The return value is (q,r)
        """
        def __truediv__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            r = self
            q = []
            while r.degree() >= other.degree(): 
                q = [r[-1] / other[-1]] + q
                r = r -\
                    Rx((r.degree() - other.degree())*[Rx.R.constant(0)] + list(q[0] * other)) 
                         
            return (Rx(q), r)

        def __rtruediv__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__truediv__(self)


        """
            Returns the r value from the Euclidean algorithm '/' above
        """
        def __mod__(self, other):
            q,r = self.__truediv__(other)
            return r

        def __rmod__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__mod__(self)

        def __eq__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            if len(self) != len(other):
                return False
            return all([self[i] == other[i] for i in range(len(self))])

        def __req__(self, other):
            return self.__eq__(other)

        def __ne__(self, other):
            return not self.__eq__(other)

        def __rne__(self, other):
            return not self.__eq__(other)
                

    return Rx 

"""
    This meta class gives correct type equality
"""
class FieldMeta(type):
    def __str__(cls):
        return "F%i" % (cls.P ** cls.N)

    def __repr__(cls):
        return str(cls)

    def __eq__(cls, other):
        if type(cls) == type(other):
            return cls.P == other.P and cls.N == other.N
        return False

"""
    Euclid's Algorithm
"""
def gcd(a,b):
    if b == 0:
        return a
    return gcd(b, a % b)

"""
    This computes an irreducible degree deg polynomial in Poly by brute force and
    it turns out is quite slow
"""
def _getIrreducible(Poly, deg):
    R = Poly.R
    #Store lower order irreducibles
    irreducibles = []
    if not isinstance(deg, int) or deg < 2:
        raise TypeError("No irreducible polynomials of degree %s" % deg)
    for i in range(2, deg-1):
        for f in Poly.monic(i):
            if all(f(x) != 0 for x in R.set()):
                # This polynomial has no zeros <=> Has no linear factors
                if all(f % g != 0 for g in irreducibles):
                    # f has no irreducible factors either
                    irreducibles += [f]

    for f in Poly.monic(deg):
        if all(f % g != 0 for g in irreducibles):
            # f has no irreducible factors either
            if all(f(x) != 0 for x in R.set()):
                # This polynomial has no zeros <=> Has no linear factors
                return f
    raise ValueError("We couldn't find an irreducible polynomial...")

"""
    This is adapted from ''Three Ways to Test Irreducibility'' by Richard P. Brent

    Brent claims that P(x) of degree d over Fp is irreducible if and only if 

    x^(p^d) = x mod P(x)

    and for all prime divisors m of d

    gcd(x^(p^(d/m)) - x, P(x)) = 1
    (This is not what Brent says. He only gives the case for p=2 but I think this is 
    the correct extension)
"""

def gcdCheck(P, primes, f):
    for m in primes:
        x = P([P.R.constant(0),P.R.constant(1)]).pow(P.R.P ** m, f) 
        g = gcd(f, x - P([P.R.constant(0),P.R.constant(1)]))
        if g.degree() != 0:
            return False
    return True


def getIrreducible(P, d):
    primes = [d//m for m in primeFactors(d)]
        
    for f in P.monic(d):
        #compute x^(p^d) by doubling
        x = P([P.R.constant(0),P.R.constant(1)]).pow(P.R.P ** d, f)
        if ((x - P([P.R.constant(0),P.R.constant(1)])) % f) == 0 and gcdCheck(P, primes, f):
            return f
    raise Exception("We couldn't find an irreducible polynomial...")

"""
    F(p,n, symbol='x') Generates the field of size p^n with symbol x as the ajoined symbol
    Set bruteForce=True if you want to compute the irreducible polynomial by brute force
    (useful in case you don't trust Brent)
"""
def F(prime, power, symbol='a', bruteForce=False):
    if power == 1:
        return _Z(prime)
    if power < 1 or not isinstance(power,int):
        raise TypeError("The power %s is not a natural number!" % str(power))
    #First we find an irreducible polynomial of degree n
    R = _Z(prime)
    Poly = _PolyOver(R,symbol,False)
    
    if bruteForce:
        ideal = _getIrreducible(Poly, power)
    else:
        ideal = getIrreducible(Poly, power)

    class Fpn(tuple,metaclass=FieldMeta):
        P = prime
        N = power
        S = symbol
        I = ideal

        @classmethod
        def set(cls):
            for i in range(cls.P ** cls.N):
                yield Fpn([R(i // (Fpn.P ** j)) for j in range(Fpn.N)]) 

        @classmethod
        def nonzero(cls):
            for i in range(1, cls.P ** cls.N):
                yield Fpn([R(i // (Fpn.P ** j)) for j in range(Fpn.N)]) 

        @classmethod
        def constant(cls, x):
            return cls(Poly.constant(x))
        
        def __new__(cls, val):
            if isinstance(val, int):
                val = [val]
            return tuple.__new__(cls, (Poly(val) % Fpn.I,))

        @classmethod
        def id(cls):
            return cls(cls.P(0))

        @property
        def val(self):
            return self[0]

        @classmethod
        def _convert(cls, other):
            if not isinstance(other, cls):
                other = cls.constant(other)
            return other

        def __add__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Fpn(self.val + other.val)

        def __radd__(self, other):
            return self.__add__(other)

        def __sub__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Fpn(self.val - other.val)

        def __rsub__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__sub__(self)

        def __neg__(self):
            return Fpn(-self.val)

        def __mul__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented

            return Fpn(self.val * other.val)

        def __rmul__(self, other):
            return self.__mul__(other)

        def __pow__(self, power):
            if not isinstance(power, int):
                raise TypeError("Cannot raise to non-integer power")
            power = power % (Fpn.P ** Fpn.N - 1)
            counter = Fpn([1])
            for i in range(0,power): #This could be improved by splitting power up in binary
                counter = counter * self
            return counter

        def __div__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented

            if other.val == 0:
                raise ZeroDivisionError()
            #Since the multiplicative group is cyclic of size
            #q-1, a** (q-2) * a = a ** (q-1) = 1
            return self * (other ** ((Fpn.P ** Fpn.N) - 2)) 

        def __rdiv__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__div__(self)

        __floordiv__ = __div__
        __rfloordiv__ = __rdiv__

        __truediv__ = __div__
        __rtruediv__ = __rdiv__
        
        def __str__(self):
            return str(self.val) 

        def __repr__(self):
            return str(self)

        def __eq__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return self.val == other.val;

        def __req__(self, other):
            return self.__eq__(other)

        def __ne__(self, other):
            return not self.__eq__(other)

        def __rne__(self, other):
            return not self.__eq__(other)

    return Fpn

def _multiCall(self, x):
        y = x[-1]
        z = [x[i] for i in range(len(x) - 1)]
        if len(z) == 1:
            z = z[0]
        return sum([self[i](z) * (y ** i) for i in range(len(self))])


def PolynomialsOver(field, n=1, symbol='x'):
    if n == 1:
        return _PolyOver(field, symbol)

    P = _PolyOver(field, symbol + str(1) )

    for i in range(1,n):
        P = _PolyOver(P, symbol + str(i+1))
        P.__call__ = _multiCall
        P.Base = field
    return P


