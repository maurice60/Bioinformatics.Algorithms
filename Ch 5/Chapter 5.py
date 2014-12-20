__author__ = 'maurice'
import sys
import os
import copy
import math
import re
from collections import deque
# from operator import add
import sqlite3
import numpy as np
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
    def backtrack():
        s = [[0 for _ in range(len(w)+1)] for _ in range(len(v)+1)]
        b = copy.deepcopy(s)
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        for i in range(1, cellsi):
            s[i][0] = s[i-1][0]
            b[i][0] = 1
        for j in range(1, cellsj):
            s[0][j] = s[0][j-1]
            b[0][j] = 2
        for i in range(1, cellsi):
            for j in range(1, cellsj):
                tes = [s[i-1][j], s[i][j-1]]
                eq = v[i-1] == w[j-1]
                if eq:
                    tes.append(s[i-1][j-1]+1)
                s[i][j] = max(tes)
                if s[i][j] == s[i-1][j]:
                    b[i][j] = 1
                if s[i][j] == s[i][j-1]:
                    b[i][j] = 2
                if eq and s[i][j] == s[i-1][j-1]+1:
                    b[i][j] = 3
        # for x in b:
        #     print x
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
            out.append(v[i-1])
        return ni, nj

    out = []
    bt = backtrack()
    x = len(v) # - 1
    y = len(w) # - 1
    while x > 0 and y > 0:
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
            # try:
            #     nodes[int(m[0])][int(m[1])] = int(m[2])
            # except KeyError:
            #     nodes[int(m[0])] = {int(m[1]):int(m[2])}
            try:
                nodes[m[0]][m[1]] = int(m[2])
            except KeyError:
                nodes[m[0]] = {m[1]:int(m[2])}
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

def globalAlignment(v, w, sig = None, bs62 = None):

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
                    continue
                if s[i][j] == tes[1]:
                    b[i][j] = 2
                    continue
                if s[i][j] == tes[2]:
                    b[i][j] = 3
        score = s[len(v)][len(w)]
        # for q in s:
        #     print q
        # print ''
        # for q in b:
        #     print q
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
    if bs62 is None:
        bs62 = blosum62DS()
    if sig is None:
        sig = -5
    sc, bt = backtrack()
    x = len(v)
    y = len(w)
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def localAlignment(v, w, sig = None, pam250 = None):

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
                    continue
                if s[i][j] == tes[2]:
                    b[i][j] = 2
                    continue
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
    if pam250 is None:
        pam250 = pam250DS()
    if sig is None:
        sig = -5
    (x, y, sc), bt = backtrack()
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def editDistance(a, b):

    def setupDict():
        cx = sqlite3.connect("../Bioinformatics.db")
        px = cx.cursor()
        px.execute("SELECT * FROM pam250")
        out = px.fetchall()
        cx.close()
        if out is None:
            return {}
        outD = {}
        for x in out:
            if x[0] == x[1]:
                outD[(str(x[0]), str(x[1]))] = 0
            else:
                outD[(str(x[0]), str(x[1]))] = -1
        return outD

    dic = setupDict()
    (_, txtA, txtB) = globalAlignment(a, b, sig=-1, bs62=dic)
    sc = 0
    for i, c in enumerate(txtA):
        if c == txtB[i]:
            continue
        sc += 1
    return sc

def fittingAlignment(v, w, sig = None, pam250 = None):
    def setupDict():
        cx = sqlite3.connect("../Bioinformatics.db")
        px = cx.cursor()
        px.execute("SELECT * FROM pam250")
        out = px.fetchall()
        cx.close()
        if out is None:
            return {}
        outD = {}
        for x in out:
            if x[0] == x[1]:
                outD[(str(x[0]), str(x[1]))] = 1
            else:
                outD[(str(x[0]), str(x[1]))] = sig
        return outD

    def backtrack():
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        score = (0,0,0)
        s1 = [0 for _ in xrange(cellsj)]
        s = [s1[:] for _ in xrange(cellsi)]
        b = [s1[:] for _ in xrange(cellsi)]

        # for i in range(1, cellsi):
        #     # s[i][0] = s[i-1][0] + sig
        #     b[i][0] = 0
        for j in range(1, cellsj):
            s[0][j] = s[0][j-1] + sig
            b[0][j] = 2
        for j in range(1, cellsj):
            for i in range(1, cellsi):
                tes = [s[i-1][j] + sig, s[i][j-1] + sig, s[i-1][j-1] + pam250[(v[i-1], w[j-1])]]
                s[i][j] = max(tes)
                if j == len(w) and s[i][j] > score[2]:
                    score = (i, j, s[i][j])
                if s[i][j] == tes[0]:
                    b[i][j] = 1
                    continue
                if s[i][j] == tes[1]:
                    b[i][j] = 2
                    continue
                if s[i][j] == tes[2]:
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
        elif bt[i][j] == 3:
            ni = i-1
            nj = j-1
            outA.append(v[i-1])
            outB.append(w[j-1])
        else:
            ni = 0
            nj = 0
        return ni, nj

    outA = []
    outB = []
    if sig is None:
        sig = -1
    if pam250 is None:
        pam250 = setupDict()
    (x, y, sc), bt = backtrack()
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def overlapAlignment(v, w):
    def setupDict():
        cx = sqlite3.connect("../Bioinformatics.db")
        px = cx.cursor()
        px.execute("SELECT * FROM pam250")
        out = px.fetchall()
        cx.close()
        if out is None:
            return {}
        outD = {}
        for x in out:
            if x[0] == x[1]:
                outD[(str(x[0]), str(x[1]))] = 1
            else:
                outD[(str(x[0]), str(x[1]))] = sig
        return outD

    def backtrack():
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        score = (0,0,0)
        s1 = [0 for _ in xrange(cellsj)]
        s = [s1[:] for _ in xrange(cellsi)]
        b = [s1[:] for _ in xrange(cellsi)]

        for j in range(1, cellsj):
            for i in range(1, cellsi):
                tes = [s[i-1][j] + sig, s[i][j-1] + sig, s[i-1][j-1] + dic[(v[i-1], w[j-1])]]
                s[i][j] = max(tes)
                if i == len(v) and s[i][j] > score[2]:
                    score = (i, j, s[i][j])
                if s[i][j] == tes[0]:
                    b[i][j] = 1
                    continue
                if s[i][j] == tes[1]:
                    b[i][j] = 2
                    continue
                if s[i][j] == tes[2]:
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
        elif bt[i][j] == 3:
            ni = i-1
            nj = j-1
            outA.append(v[i-1])
            outB.append(w[j-1])
        else:
            ni = 0
            nj = 0
        return ni, nj

    outA = []
    outB = []
    sig = -2
    dic = setupDict()
    (x, y, sc), bt = backtrack()
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def affineGlobalAlignment(v, w, sig = None, eps = None, bs62 = None):

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
        s = [0 for _ in xrange(cellsj)]
        u = [s[:] for _ in xrange(cellsi)]
        m = [s[:] for _ in xrange(cellsi)]
        l = [s[:] for _ in xrange(cellsi)]
        b = [s[:] for _ in xrange(cellsi)]
        for i in range(1, cellsi):
            # u[i][0] = u[i-1][0] + eps
            # m[i][0] = m[i-1][0] + sig
            # l[i][0] = l[i-1][0] + eps
            b[i][0] = 1
        for j in range(1, cellsj):
            # u[0][j] = u[0][j-1] + eps
            # m[0][j] = m[0][j-1] + sig
            # l[0][j] = l[0][j-1] + eps
            b[0][j] = 2
        for i in range(1, cellsi):
            for j in range(1, cellsj):
                tesl = [l[i-1][j] + eps, m[i-1][j] + sig]
                l[i][j] = max(tesl)
                tesu = [u[i][j-1] + eps, m[i][j-1] + sig]
                u[i][j] = max(tesu)
                tesm = [l[i][j], u[i][j], m[i-1][j-1] + bs62[(v[i-1], w[j-1])]]
                m[i][j] = max(tesm)
                if m[i][j] == tesm[0]:
                    b[i][j] = 1
                    continue
                if m[i][j] == tesm[1]:
                    b[i][j] = 2
                    continue
                if m[i][j] == tesm[2]:
                    b[i][j] = 3
        # for q in u:
        #     print q
        # print ''
        # for q in m:
        #     print q
        # print ''
        # for q in l:
        #     print q
        # print ''
        # for q in b:
        #     print q
        score = m[len(v)][len(w)]
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

    if bs62 is None:
        bs62 = blosum62DS()
    if sig is None:
        sig = -11
    if eps is None:
        eps = -1

    outA = []
    outB = []
    sc, bt = backtrack()
    x = len(v)
    y = len(w)
    while bt[x][y] != 0:
        x, y = outputLCS(x, y)
    outA.reverse()
    outB.reverse()
    return (sc, ''.join(outA), ''.join(outB))

def middleEdge(a, b):

    def backtrack(v, w, vec=None):
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        s = np.zeros((cellsi, cellsj))
        if vec is None:
            vec = [0]
            for i in range(1, cellsi):
                vec.append(vec[i-1] + sig)
            dn = []
            ptr = False
        else:
            dn = np.zeros((cellsi, cellsj))
            ptr = True
        for i in xrange(cellsi):
            s[i, 0] = vec[i]
        for j in range(1, cellsj):
            for i in range(1, cellsi):
                tes = [s[i-1, j] + sig, s[i, j-1] + sig, s[i-1, j-1] + dic[(v[i-1], w[j-1])]]
                t = max(tes)
                s[i, j] = t
                if ptr:
                    if t == tes[0]:
                        dn[i, j] = 1
                        continue
                    if t == tes[1]:
                        dn[i, j] = 2
                        continue
                    if t == tes[2]:
                        dn[i, j] = 3
        if ptr:
            lcol = []
        else:
            lcol = [val[-1] for val in s]
        return dn, lcol

    def track(i, j):
        if bt[i][j] == 1:
            ni = i-1
            nj = j
        elif bt[i][j] == 2:
            ni = i
            nj = j-1
        else:
            ni = i-1
            nj = j-1
        return ni, nj

    sig=-5
    p = (len(b) / 2)
    dic=blosum62DS()
    _, vec = backtrack(a, b[:p])
    bt, _ = backtrack(a, b[p:], vec)
    x = len(a)
    y = len(b[p:])
    x1 = x
    while bt[x][y] != 0:
        x, y = track(x, y)
        x1 = x

    return [(x, p), (x1, p+1)]

def findMid(a, b, dic):

    def backtrack(v, w):
        cellsi = len(v) + 1
        cellsj = len(w) + 1
        s = {0:0}
        for j in xrange(1, cellsj):
            s[j] = s[j-1] + sig
        for i in xrange(1, cellsi):
            t = s.copy()
            s = {0:sig*i}
            ch = v[i-1]
            for j in xrange(1, cellsj):
                s[j] = max([t[j] + sig, s[j-1] + sig, t[j-1] + dic[(ch, w[j-1])]])
        return s.values()

    sig=-5
    p = (len(a) / 2)
    vec1 = backtrack(a[:p], b)
    vec2 = backtrack(a[p:][::-1], b[::-1])
    maxV = 0
    xmax = 0
    for i, c in enumerate(vec1):
        sm = c + vec2[-i-1]
        if sm > maxV:
            xmax = i
            maxV = sm
    return (p, xmax)

def linearSpaceAlignment(a, b):
    dic=blosum62DS()
    # for x in dic:
    #     if x[0] == x[1]:
    #         dic[(str(x[0]), str(x[1]))] = 2
    #     else:
    #         dic[(str(x[0]), str(x[1]))] = sig
    (u, v) = findMid(a, b, dic)
    (s, t) = findMid(a[u:], b[v:], dic)
    (w, x) = findMid(a[:u], b[:v], dic)
    return (s,t), (u,v), (w,x)

def global3Way(u, v, w):

    def dic(a, b, c = None):
        if a == b:
            if b == c:
                return 1
        return 0

    def backtrack():
        cellsi = len(u) + 1
        cellsj = len(v) + 1
        cellsk = len(w) + 1
        s = np.zeros((cellsi, cellsj, cellsk)).astype('int32')
        b = np.zeros((cellsi, cellsj, cellsk)).astype('int32')

        for i in xrange(1, cellsi):
            b[i,0,0] = 1
        for j in xrange(1, cellsj):
            b[0,j,0] = 2
        for k in xrange(1, cellsk):
            b[0,0,k] = 3
        for i in xrange(1, cellsi):
            for j in xrange(1, cellsj):
                b[i,j,0] = 4
        for j in xrange(1, cellsj):
            for k in xrange(1, cellsk):
                b[0,j,k] = 5
        for i in xrange(1, cellsi):
            for k in xrange(1, cellsk):
                b[i,0,k] = 6
        for k in xrange(1, cellsk):
            for j in xrange(1, cellsj):
                for i in xrange(1, cellsi):
                    x1 = s[i-1,j,k]
                    y1 = s[i,j-1,k]
                    z1 = s[i,j,k-1]
                    x2 = s[i-1,j-1,k] #+ dic(u[i-1], v[j-1])
                    y2 = s[i,j-1,k-1] #+ dic(v[j-1], w[k-1])
                    z2 = s[i-1,j,k-1] #+ dic(u[i-1], w[k-1])
                    z3 = s[i-1,j-1,k-1] + dic(u[i-1], v[j-1], w[k-1])
                    tes = [x1,y1,z1,x2,y2,z2,z3]
                    mx = max(tes)
                    s[i,j,k] = mx
                    if mx == tes[6]:
                        b[i,j,k] = 7
                        continue
                    if mx == tes[3]:
                        b[i,j,k] = 4
                        continue
                    if mx == tes[4]:
                        b[i,j,k] = 5
                        continue
                    if mx == tes[5]:
                        b[i,j,k] = 6
                        continue
                    if mx == tes[0]:
                        b[i,j,k] = 1
                        continue
                    if mx == tes[1]:
                        b[i,j,k] = 2
                        continue
                    if mx == tes[2]:
                        b[i,j,k] = 3
                    else:
                        print '*'
        # for q in s:
        #     print q
        # print ''
        # for q in b:
        #     print q
        score = s[len(u),len(v),len(w)]
        return score, b

    def outputLCS(i, j, k):
        this = bt[i,j,k]
        # print i,j,k,this
        if this == 1:
            ni = i-1
            nj = j
            nk = k
            outA.append(u[i-1])
            outB.append('-')
            outC.append('-')
        elif this == 2:
            ni = i
            nj = j-1
            nk = k
            outA.append('-')
            outB.append(v[j-1])
            outC.append('-')
        elif this == 3:
            ni = i
            nj = j
            nk = k-1
            outA.append('-')
            outB.append('-')
            outC.append(w[k-1])
        elif this == 4:
            ni = i-1
            nj = j-1
            nk = k
            outA.append(u[i-1])
            outB.append(v[j-1])
            outC.append('-')
        elif this == 5:
            ni = i
            nj = j-1
            nk = k-1
            outA.append('-')
            outB.append(v[j-1])
            outC.append(w[k-1])
        elif this == 6:
            ni = i-1
            nj = j
            nk = k-1
            outA.append(u[i-1])
            outB.append('-')
            outC.append(w[k-1])
        else:
            ni = i-1
            nj = j-1
            nk = k-1
            outA.append(u[i-1])
            outB.append(v[j-1])
            outC.append(w[k-1])
        return ni, nj, nk

    outA = []
    outB = []
    outC = []
    sc, bt = backtrack()
    x = len(u)
    y = len(v)
    z = len(w)
    while bt[x,y,z] != 0:
        x, y, z = outputLCS(x, y, z)
    outA.reverse()
    outB.reverse()
    outC.reverse()
    return sc, (''.join(outA), ''.join(outB), ''.join(outC))

# blosum62Load(aReadT('grid.txt'))
# pam250Load(aReadT('grid.txt'))
# print coins([3,2],  24)
# print manhattanTourist(14, 11, gRead('grid.txt'))
# print LCSBacktrack(fRead('strA.txt'), fRead('strB.txt'))
# outp = topologicalOrdering(aReadT('grid.txt'))
# scr, outp = longestPathDAG(aReadT('grid.txt'), 'a', 'g')
# outp = globalAlignment(fRead('strA.txt'), fRead('strB.txt'))
# print ''
# outp = localAlignment(fRead('strA.txt'), fRead('strB.txt'))
# print editDistance(fRead('strA.txt'), fRead('strB.txt'))
# outp = fittingAlignment(fRead('strA.txt'), fRead('strB.txt'))
# outp = overlapAlignment(fRead('strA.txt'), fRead('strB.txt'))
# outp = affineGlobalAlignment(fRead('strA.txt'), fRead('strB.txt'))
# outp = middleEdge(fRead('strA.txt'), fRead('strB.txt'))
outp = linearSpaceAlignment(fRead('strA.txt'), fRead('strB.txt'))
# scr, outp = global3Way(fRead('strA.txt'), fRead('strB.txt'), fRead('grid.txt'))
# print outp
# print scr
for xoutp in outp:
    print xoutp
# print '->'.join([str(c) for c in outp])