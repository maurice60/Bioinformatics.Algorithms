__author__ = 'maurice'
from scipy.special import binom
import sys
import os
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

# out = motifEnumeration(lReadT('dna.txt'), 5, 2)
# print len(out)
# for x in out:
#     print x,
# dna = ['ttaccttAAc', 'gAtAtctgtc', 'Acggcgttcg', 'ccctAAAgag', 'cgtcAgAggt']
# dna = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTACGGGACAG']
# st = medianString([x.upper() for x in dna], 3)
# print medianString(lReadT('dna.txt'), 6)
print neighbours('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA', 3)