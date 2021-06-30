#!/usr/bin/env python3
from codes import *
from fields import *

points = [(0, 0, 5), (4, 2, 0), (2, 3, 0),\
    (2, 1, 0), (0, 0, 1), (2, 1, 1), (0, 0, 2), (2, 1, 2),\
    (0, 0, 3), (0, 0, 4), (1, 1, 0), (3, 2, 0), (1, 1, 1),\
    (1, 1, 2), (1, 1, 3), (2, 2, 0), (2, 2, 1) ]
f = F(7,1)
print("Computing min dist for %s over %s" % (str(points), str(f)), flush=True)
print(CodeFromLatticePoints(f, points, True).d(), flush=True)
