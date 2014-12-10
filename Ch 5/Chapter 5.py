__author__ = 'maurice'
import sys
import os
import copy
import math
import re
from collections import deque
import sqlite3
# import random

sys.path.append(os.path.relpath('../'))
from common import *

def blosum62Load(array):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    hdr = [s for s in array[0].split(' ') if len(s) > 0]
    for x in array[1:]:
        row = [s for s in x.split(' ') if len(s) > 0]
        for i, c in enumerate(row[1:]):
            px.execute("INSERT INTO blosum62 (aminoA, aminoB, coefficient) VALUES (?, ?, ?)", (row[0], hdr[i], int(c)))
    cx.commit()
    cx.close()

def pam250Load(array):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    hdr = [s for s in array[0].split(' ') if len(s) > 0]
    for x in array[1:]:
        row = [s for s in x.split(' ') if len(s) > 0]
        for i, c in enumerate(row[1:]):
            px.execute("INSERT INTO pam250 (aminoA, aminoB, coefficient) VALUES (?, ?, ?)", (row[0], hdr[i], int(c)))
    cx.commit()
    cx.close()

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

def topologicalOrdering(graph):

    def toList():
        nodes = {}
        for g in graph:
            m = re.findall(r'(\w+)', g)
            nodes[int(m[0])] = ([int(i) for i in m[1:]])
        return nodes

    def outE():
        inE = set()
        for m in copG.values():
            for n in m: #.keys():
                inE.add(n)
        return inE

    def noIn():
        return sorted([m for m in copG if outE().isdisjoint([m])])

    graphOrg = toList()
    copG = copy.deepcopy(graphOrg)
    outpt = []
    candidates = deque(noIn())
    while len(candidates) > 0:
        x = candidates.popleft()
        outpt.append(x)
        try:
            y = copG[x]
            del copG[x]
        except KeyError:
            y = {}
        rem = outE()
        for i in y:
            if rem.isdisjoint([i]):
                candidates.append(i)
    return outpt

def longestPathDAG(graph, s, e):

    def toList():
        nodes = {}
        for g in graph:
            m = re.findall(r'(\w+)', g)
            try:
                nodes[int(m[0])][int(m[1])] = int(m[2])
            except KeyError:
                nodes[int(m[0])] = {int(m[1]):int(m[2])}
        return nodes

    def trav(n, p, tot):
        if n == e:
            l[tot] = p[:]
            return
        try:
            for x, y in thisG[n].iteritems():
                trav(x, p + [x], tot + y)
        except KeyError:
            return

    l = {}
    thisG = toList()
    trav(s, [s], 0)
    if len(l) > 0:
        m = max(l.keys())
        return m, l[m]
    else:
        return 0, []

def blosum62(a, b):
    cx = sqlite3.connect("../Bioinformatics.db")
    px = cx.cursor()
    px.execute("SELECT coefficient FROM blosum62 WHERE aminoA = ? AND aminoB = ?", (a, b))
    out = px.fetchone()
    cx.close()
    if out is None:
        return ''
    return out[0]

def globalAlignment(v, w, sig = None, bs62 = {}):

    def blosum62DS():
        cx = sqlite3.connect("../Bioinformatics.db")
        px = cx.cursor()
        px.execute("SELECT * FROM blosum62")
        out = px.fetchall()
        cx.close()
        if out is None:
            return {}
        outD = {}
        for x in out:
            outD[(str(x[0]), str(x[1]))] = x[2]
        return outD

    def backtrack():
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        s1 = [0 for _ in xrange(cellsj)]
        s = [s1[:] for _ in xrange(cellsi)]
        b = [s1[:] for _ in xrange(cellsi)]

        for i in range(1, cellsi):
            s[i][0] = s[i-1][0] + sig
            b[i][0] = 1
        for j in range(1, cellsj):
            s[0][j] = s[0][j-1] + sig
            b[0][j] = 2
        for j in range(1, cellsj):
            for i in range(1, cellsi):
                tes = [s[i-1][j] + sig, s[i][j-1] + sig, s[i-1][j-1] + bs62[(v[i-1], w[j-1])]]
                s[i][j] = max(tes)
                if s[i][j] == tes[0]:
                    b[i][j] = 1
                if s[i][j] == tes[1]:
                    b[i][j] = 2
                if s[i][j] == tes[2]:
                    b[i][j] = 3
        score = s[len(v)][len(w)]
        return score, b

    def outputLCS(i, j):
        if bt[i][j] == 1:
            ni = i-1
            nj = j
            outA.append(v[i-1])
            outB.append('-')
        elif bt[i][j] == 2:
            ni = i
            nj = j-1
            outA.append('-')
            outB.append(w[j-1])
        else:
            ni = i-1
            nj = j-1
            outA.append(v[i-1])
            outB.append(w[j-1])
        return ni, nj

    outA = []
    outB = []
    if bs62 == {}:
        bs62 = blosum62DS()
    if sig == None:
        sig = -5
    sc, bt = backtrack()
    x = len(v)
    y = len(w)
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def localAlignment(v, w, sig = None):

    def pam250DS():
        cx = sqlite3.connect("../Bioinformatics.db")
        px = cx.cursor()
        px.execute("SELECT * FROM pam250")
        out = px.fetchall()
        cx.close()
        if out is None:
            return {}
        outD = {}
        for x in out:
            outD[(str(x[0]), str(x[1]))] = x[2]
        return outD

    def backtrack():
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        pam250 = pam250DS()
        score = (0,0,0)
        s1 = [0 for _ in xrange(cellsj)]
        s = [s1[:] for _ in xrange(cellsi)]
        b = [s1[:] for _ in xrange(cellsi)]

        for i in range(1, cellsi):
            s[i][0] = s[i-1][0] + sig
            b[i][0] = 1
        for j in range(1, cellsj):
            s[0][j] = s[0][j-1] + sig
            b[0][j] = 2
        for j in range(1, cellsj):
            for i in range(1, cellsi):
                tes = [0, s[i-1][j] + sig, s[i][j-1] + sig, s[i-1][j-1] + pam250[(v[i-1], w[j-1])]]
                s[i][j] = max(tes)
                if s[i][j] > score[2]:
                    score = (i, j, s[i][j])
                if s[i][j] == tes[1]:
                    b[i][j] = 1
                if s[i][j] == tes[2]:
                    b[i][j] = 2
                if s[i][j] == tes[3]:
                    b[i][j] = 3
        return score, b

    def outputLCS(i, j):
        if bt[i][j] == 1:
            ni = i-1
            nj = j
            outA.append(v[i-1])
            outB.append('-')
        elif bt[i][j] == 2:
            ni = i
            nj = j-1
            outA.append('-')
            outB.append(w[j-1])
        else:
            ni = i-1
            nj = j-1
            outA.append(v[i-1])
            outB.append(w[j-1])
        return ni, nj

    outA = []
    outB = []
    if sig == None:
        sig = -5
    (x, y, sc), bt = backtrack()
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def editDistance(a, b):
    (_, txtA, txtB) = globalAlignment(a, b, sig=0)
    sc = 0
    for i, c in enumerate(txtA):
        if c == txtB[i]:
            continue
        sc += 1
        print sc, c, txtB[i]
    return sc

# blosum62Load(aReadT('grid.txt'))
# pam250Load(aReadT('grid.txt'))
# print coins([9,5,3,1],  19163)
# print manhattanTourist(14, 11, gRead('grid.txt'))
# print LCSBacktrack(fRead('strA.txt'), fRead('strB.txt'))
# out = topologicalOrdering(aReadT('grid.txt'))
# x, out = longestPathDAG(aReadT('grid.txt'), 0, 44)
# out = globalAlignment(fRead('strA.txt'), fRead('strB.txt'))
# out = localAlignment(fRead('strA.txt'), fRead('strB.txt'))
print editDistance(fRead('strA.txt'), fRead('strB.txt'))
# print x
# for x in out:
#     print x
# print '->'.join([str(c) for c in out])