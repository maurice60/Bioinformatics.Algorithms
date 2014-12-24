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
    out = []
    outT = ''
    while len(lc) > 0:
        ptn = []
        ix = min(lc.keys())
        while len(lc) > 0:
            this = lc.pop(ix)
            ptn.append(this[0])
            ptn.append(this[1])
            if this[0] > this[1]:
               break
            if float(this[1]) % 2 == 0:
                ix = this[1] - 1
            else:
                ix = this[1] + 1
        ptn = [ptn[-1]] + ptn[:-1]
        out.append(ptn[:])
        outT += cycleToChromosome(ptn)
    return outT

def distance2Break(p, q):
    cycles = 0
    pBlock = set()
    qBlock = set()
    for c in p:
        cycles += 1
        for b in c:
            pBlock.add(abs(b))
    for c in q:
        cycles += 1
        for b in c:
            qBlock.add(abs(b))
    return len(pBlock.intersection(qBlock)) - cycles

# outp = 0
# for _ in plusAndMinusPermutations([1,2,3,4,5,6,7]):
#     outp += 1
# outp = greedySorting(lReadB('strA.txt'))
# print numberOfBreakpoints(lReadB('strA.txt'))
# print chromosomeToCycle(lReadB('strA.txt'))
# print cycleToChromosome(lReadB('strA.txt'))
# outp = colouredEdges(lReadBA('strA.txt'))
outp = graphToGenome(fRead('strA.txt'))
# outp = distance2Break(lReadBA('strA.txt'), lReadBA('strB.txt'))
print outp
# for xo in outp:
#     print xo,
