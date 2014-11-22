__author__ = 'maurice'
import sys
import os
import re
import copy
# import math
# import random

sys.path.append(os.path.relpath('../'))
from common import *

def genomePathWalk(dna):
    genome = dna[0]
    for x in dna[1:]:
        genome += x[-1]
    return genome

def overlapGraph(dna):
    out = []
    for s in dna:
        for h in [x for x in dna if x.startswith(s[1:])]:
            out.append((s, h))
    return sorted(out)

def deBruijnPrint(patterns):

    def deBruijn(patterns):
        out = sorted([(i[:-1], i[1:]) for i in patterns])
        return out

    lines = []
    last = 0
    lastUsed = ''
    for x in deBruijn(patterns):
        if x[0] == lastUsed:
            lines[last-1] += ',{0}'.format(x[1])
        else:
            lines.append('{0} -> {1}'.format(*x))
            last += 1
            lastUsed = x[0]
    return lines

def toList(graph):
    nodes = {}
    l = 0
    for g in graph:
        m = re.findall(r'(\w+)', g)
        nodes[m[0]] = ([i for i in m[1:]])
        l += len(m) - 1
    return l, nodes

def eulerianCycle(graph):

    def stitch(xpth, path, n):
        if xpth == []:
            if len(path[0]) == 1:
                return path[:] + list(path[0])
            else:
                return path[:] + [path[0]]
        return xpth[:n] + path[:] + xpth[n:]

    pathS = 0
    length, nodes = toList(graph)
    src = sorted(nodes.keys())[0]
    path = []
    xpth = []
    pos = 0
    while pathS < length:
        d = nodes[src]
        if len(d) != 0:
            pathS += 1
            path.append(src)
            trg = d[0]
            nodes[src].remove(trg)
            src = trg
        else:
            xpth = stitch(xpth, path, pos)
            for n, b in enumerate(xpth[pos:]):
                if len(nodes[b]) != 0:
                    src = b
                    pos += n
                    break
            path = []
    return stitch(xpth, path, pos)

def eulerianPath(graph):

    def findUnbalanced(nodes):
        dicN = {}
        for n in nodes:
            for v in nodes[n]:
                try:
                    dicN[v].append(n)
                except KeyError:
                    dicN[v] = [n]
        start = '-1'
        end = '-1'
        for n in nodes:
            try:
                if len(nodes[n]) > len(dicN[n]):
                    start = n
                elif len(nodes[n]) < len(dicN[n]):
                    end = n
            except KeyError:
                start = n
        return start, end

    def stitch(xpth, path, n):
        try:
            if path[0] != xpth[n]:
                print(path[0], xpth[n-1])
        except IndexError:
            pass
        if xpth == []:
            return path[:]
        out = xpth[:n] + path[:] + xpth[n:]
        return out

    pathS = 0
    length, nodes = toList(graph)
    src, end = findUnbalanced(nodes)
    path = []
    xpth = []
    pos = 0
    while True:
        try:
            d = nodes[src]
            if len(d) != 0:
                pathS += 1
                path.append(src)
                trg = d[0]
                nodes[src].remove(trg)
                src = trg
            else:
                if src == end:
                    pathS += 1
                    path.append(src)
                if path == []:
                    break
                xpth = stitch(xpth, path, pos)
                for n, b in enumerate(xpth[pos:]):
                    try:
                        if len(nodes[b]) != 0:
                            src = b
                            pos += n
                            break
                    except KeyError:
                        continue
                path = []
        except KeyError:
            pathS += 1
            path.append(src)
            xpth = stitch(xpth, path, pos)
            if len(xpth) > length:
                break
            for n, b in enumerate(xpth[pos:]):
                if len(nodes[b]) != 0:
                    src = b
                    pos += n
                    break
            path = []
    return xpth

def pairedComposition(text, k, d, sorted=True):
    if sorted:
        return sorted([(text[i:i+k], text[i+k+d:i+2*k+d]) for i in xrange(len(text) - 2*k - d + 1)])
    return [(text[i:i+k], text[i+k+d:i+2*k+d]) for i in xrange(len(text) - 2*k - d + 1)]

def stringSpelledByGappedPatterns(d, pairs):

    prefix = []
    suffix = []

    def toListPairs():
        for p in pairs:
            m = re.findall(r'(\w+)', p)
            prefix.append(m[0])
            suffix.append(m[1])

    toListPairs()
    s1 = genomePathWalk(prefix)
    s2 = genomePathWalk(suffix)
    k = len(prefix[0])
    ans = ''
    for i in range(k+d, len(s1)-d+1):
        q1 = s1[i:]
        if s2.startswith(s1[i:]):
            ans = s1[:i] + s2
            break
    return ans

# out = composition('0000111100101101000', 4, sort=True)
# out = genomePathWalk(aReadT('dna.txt'))
# out = overlapGraph(aReadT('dna.txt'))
# k = 4
# bstr = ['{0:04b}'.format(d) for d in xrange(2**k)]
# print bstr
# print deBruijnPrint(bstr)
# print eulerianCycle(deBruijnPrint(bstr))
# out = genomePathWalk(eulerianCycle(deBruijnPrint(bstr))) #[:2**k]
# out = deBruijnPrint(bstr)
# out = genomePathWalk(eulerianPath(deBruijnPrint(aReadT('dna.txt'))))    # (composition(fRead('dna.txt'), 4))
# out = pairedComposition('AGCAGCTGCTGCA', 2, 1, sorted=False)
out = stringSpelledByGappedPatterns(200, aReadT('dna.txt'))
print out
# print '->'.join(out)
# out = deBruijnPrint(out)
# for x in out:
#     print '({0}|{1})'.format(*x),
# with open('outf.txt', 'w') as fil:
#     fil.writelines(out)

