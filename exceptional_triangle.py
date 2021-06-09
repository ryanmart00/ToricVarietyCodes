import codes
import fields
from math import log
import multiprocessing as mp

def isPrimePower(x):
    if fields.isPrime(x):
        return True 
    if len(list(fields.primeFactors(x))) == 1:
        return True

def primePowerIterator(i):
    while True:
        if isPrimePower(i):
            yield i
        i = i+1

path = 'exceptional_triangle.dat'

def init(l):
    global lock
    lock = l

def main():
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
    l = mp.Lock()
    p = mp.Pool()
    for (j, d) in p.imap_unordered(process, primePowerIterator(i),1):
        f = open(path, 'a')
        f.write('\n%i %i' % (j, d))
        f.close()
        


def process(i):
    p = list(fields.primeFactors(i))[0]
    n = int(log(i, p))
    d = codes.CodeFromLatticePoints(fields.F(p,n),\
            [(1,1),(2,1),(1,2)], True).d()
    return (i,d)


if __name__ == '__main__':
    main()
