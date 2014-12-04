__author__ = 'maurice'
import sys
import os
import copy
import math
# import random

sys.path.append(os.path.relpath('../'))
from common import *

def paths(m, n):
    return math.factorial(m+n)/(math.factorial(m)*math.factorial(n))

def coins(denom, n):
    m = n + 1
    ans = {i:0 for i in range(m)}
    vals = [i for i in denom if i < m]
    z = 1
    while vals != []:
        for i in vals:
            ans[i] = z
        new = set(j+c for j in vals for c in denom)
        vals = [p for p in new if p < m and ans[p] == 0]
        z += 1
    return ans

def manhattanTourist(n, m, grids):
    down = grids[:n]
    right = grids[n:]
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] + right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max([s[i-1][j]+down[i-1][j],s[i][j-1]+right[i][j-1]])
    return s[n][m]

def LCSBacktrack(v, w):
    def bactrack():
        s = [[0 for _ in range(len(w))] for _ in range(len(v))]
        b = copy.deepcopy(s)
        for i in range(1, len(v)):
            for j in range(1, len(w)):
                tes = [s[i-1][j], s[i][j-1]]
                eq = v[i] == w[j]
                if eq:
                    tes.append(s[i-1][j-1]+1)
                s[i][j] = max(tes)
                if s[i][j] == s[i-1][j]:
                    b[i][j] = 1
                if s[i][j] == s[i][j-1]:
                    b[i][j] = 2
                if eq and s[i][j] == s[i-1][j-1]+1:
                    b[i][j] = 3
        return b

    def outputLCS(i, j):
        if bt[i][j] == 1:
            ni = i-1
            nj = j
        elif bt[i][j] == 2:
            ni = i
            nj = j-1
        else:
            ni = i-1
            nj = j-1
            out.append(v[i])
        return ni, nj

    out = []
    bt = bactrack()
    x = len(v) - 1
    y = len(w) - 1
    while x >= 0 and y >= 0:
        x, y = outputLCS(x, y)
    out.reverse()
    return ''.join(out)

# print coins([9,5,3,1],  19163)
# print manhattanTourist(14, 11, gRead('grid.txt'))
print LCSBacktrack(fRead('strA.txt'), fRead('strB.txt'))