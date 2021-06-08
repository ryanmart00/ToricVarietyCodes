import codes
import fields
from math import log

def isPrimePower(x):
    if fields.isPrime(x):
        return True 
    if len(list(fields.primeFactors(x))) == 1:
        return True

def main(path):
    #first check which have already been computed
    try:
        f = open(path,'r',1)
        lines = f.readlines()
        f.close()
        if len(lines) == 1:
            # if we only have one line then we just have the header
            i = 1
        else:
            i = int(lines[-1].split(' ', 1)[0]) 
    except:
        f= open(path, 'w')
        f.write('q d\n')
        f.close()
        i = 1


    while(True):
        #First find the next i
        i = i+1
        while not isPrimePower(i):
            i = i + 1
        p = list(fields.primeFactors(i))[0]
        n = int(log(i, p))
        if p ** n != i:
            raise Exception("%i ** %i is not %i" % (p,n,i))
        f = open(path,'a',1)
        f.write('\n%i %i' % (i, codes.CodeFromLatticePoints(fields.F(p,n),\
                [(1,1),(2,1),(1,2)], True).d()))
        f.close()

if __name__ == '__main__':
    main('exceptional_triangle.dat')
