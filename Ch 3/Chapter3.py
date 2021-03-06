__author__ = 'maurice'
from scipy.special import binom
import sys
import os
import math
import random

sys.path.append(os.path.relpath('../'))
from common import *

def occurs(n, a, k, t):
    return binom(n - t*(k - 1), t)/a**(t*k)

def motifEnumeration(dna, k, d):
    patterns = set()
    toTest = set()
    for str in dna:
        for i in range(len(str) - k+1):
            for j in neighbours(str[i:i+k], d):
                toTest.add(j)
    for w in toTest:
        ok = True
        for x in dna:
            if len([x[c:c+k] for c in range(len(x)-k+1) if hammingDistance(w, x[c:c+k]) <= d]) == 0:
                ok = False
                break
        if ok:
            patterns.add(w)
    return patterns

def hammingArray(dna, str):
    out = 0
    for this in dna:
        minT = len(this)
        for x in xrange(len(this)-len(str)):
            d = hammingDistance(this[x:x+len(str)], str)
            if d < minT:
                minT = d
        out += minT
    return out

def checkStrings(dna, k):
    out = set()
    for this in dna:
        for x in xrange(len(this)-k):
            out.add(this[x:x+k])
    return out

def medianString(dna, k):
    distance = len(dna)
    median = ''
    for pattern in checkStrings(dna, k):
        d = hammingArray(dna, pattern)
        if distance > d:
            distance = d
            median = pattern
    return median

def countMatrix(dna):
    global ADIC
    out = [[0 for _ in dna[0]] for _ in ADIC]
    for str in dna:
        for i, c in enumerate(str):
            out[ADIC[c]][i] += 1
    return out

def profile(matrix, pseudoCounts=True):
    if pseudoCounts:
        n = sum([x[0] + 4 for x in matrix])
        return [[(float(v)+1)/n for v in line] for line in matrix]
    n = sum([x[0] for x in matrix])
    return [[float(v)/n for v in line] for line in matrix]

def score(matrix):
    b = ['A', 'C', 'G', 'T']
    n = sum([x[0] for x in matrix])
    s = [n for _ in matrix[0]]
    a = ['' for _ in matrix[0]]
    for j, line in enumerate(matrix):
        for i, num in enumerate(line):
            if n - num < s[i]:
                s[i] = n - num
                a[i] = b[j]
    return sum(s), ''.join(a)

def entropy(profile):
    return sum([sum([-i*math.log(i, 2) for i in line if i > 0]) for line in profile])

def profileProbability(p, profile):
    global ADIC
    x = 1
    for i, c in enumerate(p):
        x *= profile[ADIC[c]][i]
    return x

def profileMostProbable(str, k, profile):
    maxP = 0
    maxS = str[0:k]
    for s in [str[i:i+k] for i in xrange(len(str) - k + 1)]:
        pro = profileProbability(s, profile)
        if maxP < pro:
            maxP = pro
            maxS = s
    return maxS

def greedyMotifSearch(dna, k, t=-1, pseudoCounts=True):
    if t == -1:
        t = len(dna)
    if t == 0:
        return []
    bestMotifs = [d[:k] for d in dna]
    motifs = bestMotifs[:]
    best, _ = score(countMatrix(bestMotifs))
    # trials = [dna[0][i:i+k] for i in xrange(len(dna[0]) - k + 1)]
    trials = composition(dna[0], k)
    for motif in trials:
        motifs[0] = motif
        for i in range(1, t):
            motifs[i] = profileMostProbable(dna[i], k, profile(countMatrix(motifs[0:i]), pseudoCounts))
        this, _ = score(countMatrix(motifs))
        if this < best:
            best = this
            bestMotifs = motifs[:]
    return bestMotifs

def randomisedMotifSearch(dna, k, t=-1, n=1):

    def searchFun(dna, k, t):
        bestMotifs = ['' for _ in dna]
        for i, m in enumerate(dna):
            d = random.randint(0, len(dna[0])-k)
            bestMotifs[i] = m[d:d+k]
        motifs = bestMotifs[:]
        bestCur = score(countMatrix(bestMotifs))
        while True:
            prof = profile(countMatrix(motifs))
            motifs = [profileMostProbable(x, k, prof) for x in dna]
            thisS = score(countMatrix(motifs))
            if thisS < bestCur:
                bestMotifs = motifs[:]
                bestCur = thisS
            else:
                return bestMotifs, bestCur

    if t == -1:
        t = len(dna)
    if t == 0:
        return []
    bestRun = 0
    bestOut = []
    for _ in xrange(n):
        thisOut, thisRun = searchFun(dna, k, t)
        if thisRun < bestRun or bestOut == []:
            bestOut = thisOut[:]
            bestRun = thisRun
    return bestOut

def gibbsSampler(dna, k, t=-1, n=1):

    def randomKmer(s, k, p):
        motifs = [s[x:x+k] for x in xrange(len(s)-k+1)]
        prob = [(profileProbability(j, p), j) for j in motifs]
        ps = sum([x for (x, _) in prob])
        prob = sorted([(x/ps, y) for (x, y) in prob], reverse=True)
        die = []
        d = 1
        for p in prob:
            d -= p[0]
            die.append((d, p[1]))
        z = random.random()
        for c in die:
            if c[0] < z:
                return c[1]
        return die[-1][1]

    def sampler(dna, k, t, n):
        bestMotifs = ['' for _ in dna]
        for i, m in enumerate(dna):
            d = random.randint(0, len(dna[0])-k)
            bestMotifs[i] = m[d:d+k]
        motifs = bestMotifs[:]
        bestCur = score(countMatrix(bestMotifs))
        for j in xrange(n):
            i = random.randint(0, t-1)
            prof = profile(countMatrix([m for z, m in enumerate(motifs) if z != i]))
            motifs[i] = randomKmer(dna[i], k, prof)
            thisS = score(countMatrix(motifs))
            if thisS < bestCur:
                bestMotifs = motifs[:]
                bestCur = thisS
        return bestMotifs, bestCur

    if t == -1:
        t = len(dna)
    if t == 0:
        return []
    bestRun = 0
    bestOut = []
    for _ in xrange(20):
        print '.',
        thisOut, thisRun = sampler(dna, k, t, n)
        if thisRun < bestRun or bestOut == []:
            bestOut = thisOut[:]
            bestRun = thisRun
    print bestRun
    return bestOut

# out = entropy(profile(countMatrix(aReadT('dna.txt'))))
# out = motifEnumeration(aReadT('dna.txt'), 5, 2)
# dna = ['ttaccttAAc', 'gAtAtctgtc', 'Acggcgttcg', 'ccctAAAgag', 'cgtcAgAggt']
# dna = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTACGGGACAG']
# st = medianString([x.upper() for x in dna], 3)
# print medianString(aReadT('dna.txt'), 6)
# print neighbours('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 3)
# out = profileMostProbable(fRead('str.txt'), 7, aReadF('dna.txt'))
# out = aReadF('dna.txt')
out = greedyMotifSearch(aReadT('dna.txt'), 12)
# out = randomisedMotifSearch(aReadT('dna.txt'), 15, n=1000)
# out = gibbsSampler(aReadT('dna.txt'), 15, n=1000)
# print out
for x in out:
    print x
