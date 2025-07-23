from sage.all import *
from random import randint
import numpy
import scipy.stats

for delta in range(760, 30, -1):
    # delta = 290
    n = 512 - delta
    d = 768 - delta

    Sigma = 1
    q = 7681
    gh = lambda Beta: sqrt(Beta / (2 * pi * exp(1)))
    Alpha = lambda Beta: gh(Beta) ** (2 / (Beta - 1))
    for Beta in range(d, 2, -1):
        LHS = sqrt(Beta) * Sigma
        RHS = sqrt(Alpha(Beta)) ** (2 * Beta - d - 1) * q ** (n / d)
        if LHS > RHS:
            print(delta, Beta + 1)
            break
    # print(Beta, LHS.n(), RHS.n(), Alpha(Beta).n(), gh(Beta).n())

# NTTRU
# delta beta
# 0     139
# 290   119


# s = []
# for i in range(1024):
#     r = 0
#     for _ in range(2):
#         r += randint(0,1) - randint(0,1)
#     s.append(r)
# print(s[:20])
# print(numpy.std(vector(s)))
