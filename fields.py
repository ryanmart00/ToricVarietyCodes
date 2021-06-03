def isPrime(p):
    if not isinstance(p, int):
        return False
    if p <= 1:
        return False
    from math import sqrt
    for i in range(2, int(sqrt(p))+1):
        if p % i == 0:
            return False
    return True

"""
    This meta class gives correct type equality
"""
class PrimeMeta(type):
    def __repr__(cls):
        return "Z/%iZ" % cls.p

    def __eq__(cls, other):
        if type(cls) == type(other):
            return cls.p == other.p
        return False

"""
    Returns the class for the field Z/pZ for the given value of p. Note that the 
    return type is a class, not an object of that class. In particular, to create 
    specific elements of Z/pZ call
    R = Z(p)
    r = R(x) # then r = x mod p 
"""
def Z(prime):
    if not isPrime(prime):
        if not isPrime(prime):
            raise TypeError("p must be prime!")

    class Zp(tuple,metaclass=PrimeMeta):
        p = prime

        @classmethod
        def set(cls):
            for i in range(0, prime):
                yield cls(i)
        
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

        def _convert(self, other):
            if type(self) == type(other):
                return other
            return type(self)(other)

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

    return Zp

class PolyMeta(type):
    def __repr__(cls):
        return "%s[%s]" % (cls.R, cls.s)

    def __eq__(cls, other):
        if type(cls) == type(other):
            return cls.R == other.R
        return False

def PolyOver(Ring, symbol='x'):
    
    class Rx(tuple, metaclass=PolyMeta):
        R=Ring
        s=symbol

        @classmethod
        def set(cls, deg):
            from itertools import product
            
        
        def __new__(cls, vals):
            converted = [val if isinstance(val, Ring) else Ring(val)\
                    for val in vals]
            #trim high order zeros
            while converted[-1] == 0:
                converted = converted[:-1]
            return tuple.__new__(cls, converted)

        @property
        def vals(self):
            return self

        def degree(self):
            return len(self) - 1

        def _convert(self, other):
            if not type(other) == type(self):
                try:
                    return type(self)(other)
                except:
                    pass
                try:
                    return type(self)([other])
                except:
                    pass
                raise TypeError("Cannot convert an element of %s to an element of %s" %\
                        (str(type(other)),str(type(self))))
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
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            p1 = list(self)+max(0,len(other) - len(self)) * [0] 
            p2 = list(other)+max(0,len(self) - len(other)) * [0] 
            return Rx([p1[i] - p2[i] for i in range(len(p1))])

        def __rsub__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return other.__sub__(self)

        def __neg__(self):
            return Rx([-a for a in self])

        def __mul__(self, other):
            try:
                other = self._convert(other)
            except:
                return NotImplemented
            return Rx([sum([self[i-j]*other[j]\
                for j in range(max(0,i-len(self)+1),1+min(i,len(other)-1))])\
                for i in range(0,len(self) + len(other)-1)])

        def __rmul__(self, other):
            return self.__mul__(other)

        def __str__(self):
            if len(self) == 0:
                return "0"
            s = str(self[0])
            if len(self) == 1:
                return s
            s += ' + ' + str(self[1]) + symbol
            for i in range(2, len(self)):
                s += " + " + str(self[i]) + symbol + "^" + str(i);
            return s

        def __repr__(self):
            return str(self)

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
                r = r - ((r.degree() - other.degree())*[0]\
                        + list(q[0] * other)) 
                         
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
                

    return Rx 

"""
    F(p,n, symbol='x') Generates the field of size p^n with symbol x as the ajoined symbol
"""
def F(p,n, symbol='x'):
    pass
