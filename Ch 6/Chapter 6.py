__author__ = 'maurice'
import sys
import os
# import copy
# import math
from collections import deque
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

def numberOfBreakpoints(pIn):

    xo = 0
    ans = 0
    n1 = len(pIn) + 1
    for _, x in enumerate(pIn):
        if x - xo != 1 or -xo + x != 1:
            ans += 1
        xo = x
    if n1 - xo != 1 or -xo + n1 != 1:
        ans += 1
    return ans

def chromosomeToCycle(chromosome, printOut=True):
    nodes = []
    for c in chromosome:
        if c > 0:
            nodes.append(2*c-1)
            nodes.append(2*c)
        else:
            nodes.append(-2*c)
            nodes.append(-2*c-1)
    if printOut:
        st = '('
        for o in nodes:
            st += '{:d} '.format(o)
        st = st.rstrip(' ')+')'
        return st
    else:
        return nodes

def cycleToChromosome(nodes, printOut=True):
    cycle = deque(nodes)
    chromosome = []
    for _ in xrange(len(cycle)/2):
        n1 = cycle.popleft()
        n2 = cycle.popleft()
        if n1 < n2:
            chromosome.append(n2/2)
        else:
            chromosome.append(-n1/2)
    if printOut:
        st = '('
        for o in chromosome:
            st += '{:+d} '.format(o)
        st = st.rstrip(' ')+')'
        return st
    else:
        return chromosome

def colouredEdges(p, printOut=True):
    edges = set()
    for c in p:
        n = chromosomeToCycle(c, printOut=False)
        nodes = deque(n[1:] + [n[0]])
        for _ in xrange(len(c)):
            edges.add((nodes.popleft(), nodes.popleft()))
    if printOut:
        st = ''
        for o in edges:
            st += '{}, '.format(o)
        st = st.rstrip(', ')
        return st
    else:
        return edges

def blackEdges(p):
    edges = set()
    for c in p:
        n = deque(chromosomeToCycle(c, printOut=False))
        for _ in xrange(len(c)):
            edges.add((n.popleft(), n.popleft()))
    return edges

def graphToGenome(edges):

    try:
        assert isinstance(edges, set)
    except AssertionError:
        tx = edges
        edges = set()
        for w in re.findall('([0-9]+\W+[0-9]+)', tx):
            x, y = w.split(',')
            edges.add((int(x), int(y)))
    lc = {e[0]:e for e in edges}
    ld = {e[1]:e for e in edges}
    out = []
    outT = ''
    # drn = True
    while len(lc) > 0:
        ptn = []
        ix = min(lc.keys())
        if float(ix) % 2 == 0:
            fnsh = ix - 1
        else:
            fnsh = ix + 1
        while len(lc) > 0:
            try:
                this = lc.pop(ix)
                ld.pop(this[1])
                ptn.append(this[0])
                ptn.append(this[1])
                if this[1] == fnsh:
                    break
                if float(this[1]) % 2 == 0:
                    ix = this[1] - 1
                else:
                    ix = this[1] + 1
            except KeyError:
                this = ld.pop(ix)
                lc.pop(this[0])
                ptn.append(this[1])
                ptn.append(this[0])
                if this[0] == fnsh:
                    break
                if float(this[0]) % 2 == 0:
                    ix = this[0] - 1
                else:
                    ix = this[0] + 1
        ptn = [ptn[-1]] + ptn[:-1]
        out.append(ptn[:])
        outT += cycleToChromosome(ptn)
    return outT

def distance2Break(p, q):
    cycles = 0
    size = 0
    for sz in p:
        size += len(sz)
    pBlock = {}
    for b in colouredEdges(p, printOut=False):
        pBlock[b[0]] = b[1]
        pBlock[b[1]] = b[0]
    qBlock = {}
    for b in colouredEdges(q, printOut=False):
        qBlock[b[0]] = b[1]
        qBlock[b[1]] = b[0]
    st = 0
    ix = 0
    while len(pBlock) > 0:
        if st == 0:
            st = pBlock.keys()[0]
            ix = st
            cycles += 1
        this = pBlock.pop(ix)
        pBlock.pop(this)
        that = qBlock.pop(this)
        qBlock.pop(that)
        if that == st:
            st = 0
        else:
            ix = that
    return size - cycles

def onGenomeGraph2Break(genomeGraph, i, i1, j, j1):
    try:
        assert isinstance(genomeGraph, set)
    except AssertionError:
        tx = genomeGraph
        genomeGraph = set()
        for w in re.findall('([0-9]+\W+[0-9]+)', tx):
            x, y = w.split(',')
            genomeGraph.add((int(x), int(y)))
    out = set()
    for x in genomeGraph:
        if x[0] == i and x[1] == i1:
            out.add((i, j))
        elif x[0] == i1 and x[1] == i:
            out.add((j,i))
        elif x[0] == j and x[1] == j1:
            out.add((i1, j1))
        elif x[0] == j1 and x[1] == j:
            out.add((j1, i1))
        else:
            out.add(x)
    return out

def onGenome2Break(p, i, i1, j, j1):
    a = colouredEdges(p, printOut=False)
    # a = colouredEdges(p, printOut=False).union(blackEdges(p))
    b = onGenomeGraph2Break(a, i, i1, j, j1)
    return graphToGenome(b)

def shortestRearrangementScenario(p, q):
    def graph():
        pBlock = {}
        for b in red:
            pBlock[b[0]] = b[1]
            pBlock[b[1]] = b[0]
        qBlock = {}
        for b in blue:
            qBlock[b[0]] = b[1]
            qBlock[b[1]] = b[0]
        st = 0
        ix = 0
        lthis = 0
        cand = ()
        while len(pBlock) > 0:
            if st == 0:
                st = pBlock.keys()[0]
                ix = st
                # cycles += 1
                lthis = 0
            this = pBlock.pop(ix)
            pBlock.pop(this)
            that = qBlock.pop(this)
            toth = qBlock.pop(that)
            lthis += 1
            if that == st:
                st = 0
                if lthis > 1:
                    cand = (toth, that)
            else:
                ix = that
        return cand

    red = colouredEdges(p, printOut=False)
    blue = colouredEdges(q, printOut=False)
    print graphToGenome(red)
    cand = graph()
    while cand != ():
        i = 0
        i1 = 0
        for c in red:
            if c[0] == cand[0]:
                i = c[1]
            if c[1] == cand[0]:
                i = c[0]
            if c[0] == cand[1]:
                i1 = c[1]
            if c[1] == cand[1]:
                i1 = c[0]
        red = onGenomeGraph2Break(red, i, cand[0], i1,  cand[1])
        print graphToGenome(red)
        cand = graph()

def sharedKMers(a, b, k):
    aDic = {}
    for j in xrange(len(a) - k + 1):
        aDic[j] =  a[j:j+k]
    bDic = {}
    for j in xrange(len(b) - k + 1):
        pat = b[j:j+k]
        patR = reverseComplement(pat)
        try:
            bDic[pat].append(j)
        except KeyError:
            bDic[pat] = [j]
        try:
            bDic[patR].append(j)
        except KeyError:
            bDic[patR] = [j]
    out = []
    for i, t in aDic.iteritems():
        if bDic.has_key(t):
            for j in bDic[t]:
                out.append((i, j))
    return out

# outp = 0
# for _ in plusAndMinusPermutations([1,2,3,4,5,6,7]):
#     outp += 1
# outp = greedySorting(lReadB('strA.txt'))
# print numberOfBreakpoints(lReadB('strA.txt'))
# outp = chromosomeToCycle(lReadB('strA.txt'))
# print cycleToChromosome(lReadB('strA.txt'))
# outp = colouredEdges(lReadBA('strA.txt'))
# outp = blackEdges(lReadBA('strA.txt'))
# outp = graphToGenome(fRead('strA.txt'))
# outp = distance2Break(lReadBA('strA.txt'), lReadBA('strB.txt'))
# outp = onGenomeGraph2Break(fRead('strA.txt'),10, 12, 3, 5)
# sto = ''
# for o in outp:
#     sto += '{}, '.format(o)
# sto = sto.rstrip(' ').rstrip(',')
# print sto
# outp = onGenome2Break(lReadBA('strA.txt'), 105, 107, 74, 76)
# shortestRearrangementScenario(lReadBA('strA.txt'), lReadBA('strB.txt'))
outp = sharedKMers(fRead('../E-coli.txt'), tRead('../Salmonella_enterica.txt'), 30)
print len(outp)
# for xo in outp:
#     print xo
