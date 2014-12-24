__author__ = 'maurice'

import re

ADIC = {'A':0, 'C':1, 'G':2, 'T':3}

def fRead(fileName):
    with open(fileName, 'r') as x:
        txt = x.read()
    return txt

def tRead(fileName):
    fd = open(fileName)
    st = ''.join([line.rstrip('\n').rstrip('\r') for line in fd if not line.startswith('>g')])
    fd.close()
    return st

def lReadT(fileName):
    x = open(fileName, 'r')
    txt = x.read()
    x.close()
    out = []
    for i in txt.split(" "):
        out.append(i)
    return out

def aReadT(fileName):
    with open(fileName, 'r') as x:
        return [line.rstrip('\n').rstrip('\r') for line in x]

def aReadF(fileName):
    with open(fileName, 'r') as x:
        out = []
        raw = [l for l in [line.rstrip('\n').rstrip('\r').split() for line in x]]
        for l in raw:
            out.append([float(a) for a in l])
    return out

def gRead(fileName):
    with open(fileName, 'r') as x:
        raw = [l for l in [line.rstrip('\n').rstrip('\r').split() for line in x]]
        out = [[int(a) for a in l] for l in raw if l != ['-']]
    return out

def lReadB(fileName):
    out = []
    with open(fileName, 'r') as x:
        for line in x:
            raw = line.lstrip('(').rstrip(')').rstrip('\n').rstrip('\r').split()
            out = [int(a) for a in raw]
    return out

def lReadBA(fileName):
    out = []
    raw = ''
    with open(fileName, 'r') as x:
        for line in x:
            raw = [l for l in line.lstrip('(').rstrip(')').rstrip('\n').rstrip('\r').split(')(')]
    for s in raw:
        out.append([int(a) for a in s.split(' ')])
    return out

def hammingDistance(p, q):
    d = 0
    for i, c in enumerate(p):
        if q[i] != c:
            d += 1
    return d

def neighbours(pattern, d):
    def neighboursA(pattern, d):
        if len(pattern) == 1:
            return ['A', 'C', 'G', 'T']
        neighbourhood = set()
        suffixNeighbours = neighboursA(pattern[1:], d)
        for text in suffixNeighbours:
            if hammingDistance(pattern[1:], text)  < d:
                for x in ['A', 'C', 'G', 'T']:
                    neighbourhood.add(x+text)
            else:
                neighbourhood.add(pattern[0]+text)
        return neighbourhood

    if d == 0:
        return pattern
    return neighboursA(list(pattern), d)

def composition(dna, k, sort=False, reverse=False):
    if sort:
        return sorted([dna[i:i+k] for i in xrange(len(dna) - k + 1)], reverse=reverse)
    else:
        return [dna[i:i+k] for i in xrange(len(dna) - k + 1)]
    