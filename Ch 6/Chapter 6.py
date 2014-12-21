__author__ = 'maurice'
import sys
import os
# import copy
# import math
# import re
# from collections import deque
# from operator import add
# import sqlite3
import numpy as np
# import random
from itertools import permutations, product

sys.path.append(os.path.relpath('../'))
from common import *

def plusAndMinusPermutations(items):
    for p in permutations(items):
        for signs in product([-1,1], repeat=len(items)):
            yield [a*sign for a,sign in zip(p,signs)]

def greedySorting(pIn):

    def po(x):
        st = '('
        for o in x:
            st += '{:+d} '.format(o)
        out.append(st.rstrip(' ')+')')

    out = []
    p = np.array(pIn)
    approxReverseDistance = 0
    for k in xrange(len(p)):
        m = k+1
        if m != abs(p[k]):
            s = [-p[k]]
            for i in xrange(k+1, len(p)):
                s.append(-p[i])
                if m == abs(p[i]):
                    break
            s.reverse()
            for i, j in enumerate(s):
                p[k+i] = j
            approxReverseDistance += 1
            po(p)
        if m == -p[k]:
            p[k] = -p[k]
            approxReverseDistance += 1
            po(p)

    with open('outf.txt', 'w') as fil:
        fil.writelines('\n'.join(out))

    return approxReverseDistance

# outp = 0
# for _ in plusAndMinusPermutations([1,2,3,4,5,6,7]):
#     outp += 1
outp = greedySorting(lReadB('strA.txt'))
print outp