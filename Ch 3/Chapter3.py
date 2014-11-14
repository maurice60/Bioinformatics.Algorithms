__author__ = 'maurice'
from scipy.special import binom
import sys
import os
import math

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

def profile(matrix):
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

out = entropy(profile(countMatrix(aReadT('dna.txt'))))
# out = motifEnumeration(aReadT('dna.txt'), 5, 2)
print out
# for x in out:
#     print x
# dna = ['ttaccttAAc', 'gAtAtctgtc', 'Acggcgttcg', 'ccctAAAgag', 'cgtcAgAggt']
# dna = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTACGGGACAG']
# st = medianString([x.upper() for x in dna], 3)
# print medianString(aReadT('dna.txt'), 6)
# print neighbours('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 3)
